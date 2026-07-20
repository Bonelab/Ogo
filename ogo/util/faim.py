"""Small FAIM/N88 command adapter used by Ogo FE CLIs."""

import csv
import os
from pathlib import Path
import re
import subprocess
from typing import Dict, Iterable, List, Optional, Sequence


_AXIS_INDEX = {"x": 0, "y": 1, "z": 2}


def _optional_path(value):
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    return Path(text).expanduser()


def resolve_faim_bin_dir(install_root=None, bin_dir=None):
    """Return the directory that likely contains the FAIM/N88 command-line tools."""
    explicit = _optional_path(bin_dir)
    if explicit is not None:
        return explicit

    root = _optional_path(install_root)
    if root is None:
        return None

    candidates = [
        root / "Contents" / "MacOS",
        root / "bin" / "64bit",
        root / "bin",
        root,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return root


def resolve_faim_command(command_name, install_root=None, bin_dir=None, explicit_path=None):
    """Resolve one FAIM/N88 command from an explicit path, install root, or PATH."""
    explicit = _optional_path(explicit_path)
    if explicit is not None:
        return str(explicit)

    resolved_bin_dir = resolve_faim_bin_dir(install_root=install_root, bin_dir=bin_dir)
    if resolved_bin_dir is not None:
        candidate = resolved_bin_dir / command_name
        if candidate.exists():
            return str(candidate)

    return command_name


def resolve_faim_solver_command(install_root=None, bin_dir=None, explicit_path=None):
    """Resolve the solver command, accepting either ``faim`` or the newer solver names."""
    explicit = _optional_path(explicit_path)
    if explicit is not None:
        return str(explicit)

    resolved_bin_dir = resolve_faim_bin_dir(install_root=install_root, bin_dir=bin_dir)
    if resolved_bin_dir is not None:
        for name in ("faim", "n88solver_slt", "n88solver_spt", "n88solver_sla"):
            candidate = resolved_bin_dir / name
            if candidate.exists():
                return str(candidate)

    return "faim"


def build_faim_env(bin_dir=None, license_dir=None, env_overrides=None):
    """Build the environment used when launching FAIM/N88 tools."""
    env = os.environ.copy()
    if env_overrides:
        env.update({str(k): str(v) for k, v in env_overrides.items()})

    resolved_bin_dir = _optional_path(bin_dir)
    if resolved_bin_dir is not None:
        existing_path = env.get("PATH", "")
        env["PATH"] = os.pathsep.join(
            [str(resolved_bin_dir)] + ([existing_path] if existing_path else [])
        )

    resolved_license_dir = _optional_path(license_dir)
    if resolved_license_dir is not None:
        env["NUMERICS88_LICENSE_DIR"] = str(resolved_license_dir)

    return env


def wrap_conda_command(cmd, conda_env=None, conda_executable="conda"):
    """Wrap a command in ``conda run`` when the user requested a named environment."""
    if conda_env and str(conda_env).strip() and str(conda_env).strip() != "current":
        return [conda_executable, "run", "-n", str(conda_env)] + list(cmd)
    return list(cmd)


def run_command(
    cmd,
    conda_env=None,
    conda_executable="conda",
    env=None,
    dry_run=False,
    capture_output=False,
    check=True,
):
    """Run one external command and optionally return stdout."""
    full_cmd = wrap_conda_command(cmd, conda_env=conda_env, conda_executable=conda_executable)
    if dry_run:
        print("[dry-run] " + " ".join(str(part) for part in full_cmd))
        return ""

    stdout = subprocess.PIPE if capture_output else None
    stderr = subprocess.PIPE if capture_output else None
    completed = subprocess.run(
        [str(part) for part in full_cmd],
        stdout=stdout,
        stderr=stderr,
        universal_newlines=True,
        env=env,
        check=check,
    )
    return completed.stdout if capture_output and completed.stdout is not None else ""


def _read_single_row_csv(path):
    if path is None or not Path(path).exists():
        return {}
    with open(str(path), newline="") as handle:
        rows = list(csv.DictReader(handle))
    return rows[0] if rows else {}


def _safe_float(value):
    if value in (None, ""):
        return ""
    try:
        return float(value)
    except (TypeError, ValueError):
        return value


def parse_pistoia_text(path):
    """Extract the main scalar values from a text Pistoia report."""
    if path is None or not Path(path).exists():
        return {}

    text = Path(path).read_text()
    out = {}
    patterns = {
        "critical_volume_pct": r"Critical volume \(%\):\s+([-\d.E+]+)",
        "critical_ees": r"Critical EES:\s+([-\d.E+]+)",
        "ees_at_crit_vol": r"EES at vol_crit:\s+([-\d.E+]+)",
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            out[key] = _safe_float(match.group(1))

    match = re.search(
        r"Failure load \(RF \* factor\):\s+([-\d.E+]+)\s+([-\d.E+]+)\s+([-\d.E+]+)",
        text,
        flags=re.IGNORECASE,
    )
    if match:
        out["failure_load_x_N"] = _safe_float(match.group(1))
        out["failure_load_y_N"] = _safe_float(match.group(2))
        out["failure_load_z_N"] = _safe_float(match.group(3))

    match = re.search(
        r"Axial stiffness:\s+([-\d.E+]+)\s+([-\d.E+]+)\s+([-\d.E+]+)",
        text,
        flags=re.IGNORECASE,
    )
    if match:
        out["stiffness_x_N_per_mm"] = _safe_float(match.group(1))
        out["stiffness_y_N_per_mm"] = _safe_float(match.group(2))
        out["stiffness_z_N_per_mm"] = _safe_float(match.group(3))

    return out


def _read_n88_node_coordinates_and_sets(model_file):
    """Read node coordinates and named node sets from an n88model file."""
    try:
        from netCDF4 import Dataset
    except ImportError as exc:
        raise RuntimeError(
            "netCDF4 is required to convert profile target displacements from percent to mm."
        ) from exc

    with Dataset(str(model_file), "r") as root:
        coords = root.groups["Parts"].groups["Part1"].variables["NodeCoordinates"][:]
        node_sets = {}
        sets_group = root.groups.get("Sets")
        if sets_group is not None and "NodeSets" in sets_group.groups:
            node_sets_group = sets_group.groups["NodeSets"]
            node_sets.update(
                {
                    name: group.variables["NodeNumber"][:]
                    for name, group in node_sets_group.groups.items()
                    if "NodeNumber" in group.variables
                }
            )
        constraints_group = root.groups.get("Constraints")
        if constraints_group is not None:
            node_sets.update(
                {
                    name: group.variables["NodeNumber"][:]
                    for name, group in constraints_group.groups.items()
                    if "NodeNumber" in group.variables
                }
            )
    return coords, node_sets


def _first_available_node_set_centroid(coords, node_sets, set_names):
    """Return the centroid for the first available node-set or constraint name."""
    for set_name in set_names:
        if set_name in node_sets:
            return _node_set_centroid(coords, node_sets, set_name)
    raise ValueError(
        "None of the expected node sets were found in the n88model: {}".format(
            ", ".join(set_names)
        )
    )


def _node_set_centroid(coords, node_sets, set_name):
    if set_name not in node_sets:
        raise ValueError(f"Node set '{set_name}' was not found in the n88model.")

    node_numbers = [int(node_number) for node_number in node_sets[set_name]]
    if not node_numbers:
        raise ValueError(f"Node set '{set_name}' is empty.")

    totals = [0.0, 0.0, 0.0]
    for node_number in node_numbers:
        coord = coords[node_number - 1]
        for axis_index in range(3):
            totals[axis_index] += float(coord[axis_index])
    return [total / len(node_numbers) for total in totals]


def _coordinate_axis_span(coords, axis_index):
    """Return the full model coordinate span along one physical axis."""
    values = [float(coord[axis_index]) for coord in coords]
    if not values:
        raise ValueError("Cannot infer characteristic length from a model with no nodes.")
    length = max(values) - min(values)
    if length <= 0:
        raise ValueError("Could not infer a positive model coordinate span.")
    return length


def infer_profile_characteristic_length_mm(model_file, report_profile, failure_axis):
    """Infer the physical length used to convert percent strain endpoints to mm."""
    profile = str(report_profile or "generic").strip().lower()
    axis_index = _AXIS_INDEX[str(failure_axis).lower()]
    coords, node_sets = _read_n88_node_coordinates_and_sets(model_file)

    if profile == "spine":
        first = _first_available_node_set_centroid(
            coords, node_sets, ["body_top", "top_displacement", "convergence_set"]
        )
        second = _first_available_node_set_centroid(
            coords, node_sets, ["body_bottom", "bottom_fixed_z"]
        )
    elif profile == "femur":
        return _coordinate_axis_span(coords, axis_index)
    else:
        return ""

    length = abs(first[axis_index] - second[axis_index])
    if length <= 0:
        raise ValueError(
            f"Could not infer a positive {profile} characteristic length along {failure_axis}."
        )
    return length


def _percent_to_mm(percent, characteristic_length_mm):
    return abs(float(percent)) * float(characteristic_length_mm) / 100.0


def _first_available_constraint_group(root, names):
    constraints = root.groups.get("Constraints")
    if constraints is None:
        return None
    for name in names:
        group = constraints.groups.get(name)
        if group is not None and "Value" in group.variables:
            return group
    return None


def read_prescribed_displacement(model_file, report_profile="generic"):
    """Read the active prescribed displacement from the model constraints."""
    try:
        from netCDF4 import Dataset
    except ImportError as exc:
        raise RuntimeError("netCDF4 is required to read n88model constraints.") from exc

    profile = str(report_profile or "generic").strip().lower()
    names = ["top_displacement", "convergence_set"]
    with Dataset(str(model_file), "r") as root:
        group = _first_available_constraint_group(root, names)
        if group is None:
            return ""
        values = group.variables["Value"][:]
        return float(values[0]) if len(values) else ""


def set_prescribed_displacement_from_percent(
    model_file,
    *,
    report_profile,
    failure_axis,
    target_displacement_percent,
    displacement_sign=-1.0,
):
    """Set model displacement constraints to a profile percent endpoint."""
    try:
        from netCDF4 import Dataset
    except ImportError as exc:
        raise RuntimeError("netCDF4 is required to update n88model constraints.") from exc

    characteristic_length = infer_profile_characteristic_length_mm(
        model_file, report_profile, failure_axis
    )
    target_displacement_mm = _percent_to_mm(target_displacement_percent, characteristic_length)
    signed_displacement = float(displacement_sign) * target_displacement_mm

    with Dataset(str(model_file), "r+") as root:
        for name in ("top_displacement", "convergence_set"):
            group = _first_available_constraint_group(root, [name])
            if group is not None:
                values = group.variables["Value"]
                values[:] = signed_displacement

    return {
        "characteristic_length_mm": characteristic_length,
        "target_displacement_mm": target_displacement_mm,
        "applied_displacement_mm": signed_displacement,
    }


def write_results_csv(
    *,
    output_file,
    model_file,
    analysis_file,
    pistoia_file,
    loads_csv,
    pistoia_csv,
    analysis_var,
    pistoia_vars,
    failure_axis,
    applied_displacement=None,
    target_displacement=0.2,
    report_profile="generic",
    warnings=None,
):
    """Write a compact, user-facing CSV with loads and stiffness."""
    pistoia_meta = parse_pistoia_text(pistoia_file)
    loads = _read_single_row_csv(loads_csv)
    pistoia = _read_single_row_csv(pistoia_csv)

    pistoia_vars = list(pistoia_vars or [])
    failure_var = pistoia_vars[0] if pistoia_vars else ""
    stiffness_var = pistoia_vars[1] if len(pistoia_vars) > 1 else ""

    failure_key = "failure_load_{}_N".format(failure_axis)
    stiffness_key = "stiffness_{}_N_per_mm".format(failure_axis)
    reaction_force = _safe_float(loads.get(analysis_var))
    pistoia_failure_load = _safe_float(
        pistoia.get(failure_var) or pistoia_meta.get(failure_key)
    )
    stiffness = _safe_float(
        pistoia.get(stiffness_var) or pistoia_meta.get(stiffness_key)
    )
    pistoia_failure_load_magnitude = (
        abs(float(pistoia_failure_load)) if pistoia_failure_load != "" else ""
    )

    # If the linear model was solved at a different prescribed displacement,
    # scale the reaction force to the requested reporting endpoint.
    profile = str(report_profile or "generic").strip().lower()
    characteristic_length = ""
    target_displacement_mm = ""
    rescale_factor = ""
    load_at_target = ""
    computed_stiffness = ""
    try:
        applied = abs(float(applied_displacement))
        target = abs(float(target_displacement))
    except (TypeError, ValueError):
        applied = ""
        target = ""

    if applied != "" and target != "":
        if profile in ("femur", "spine"):
            characteristic_length = infer_profile_characteristic_length_mm(
                model_file, profile, failure_axis
            )
            target_displacement_mm = _percent_to_mm(target, characteristic_length)
        else:
            target_displacement_mm = target
        if applied > 0 and reaction_force != "":
            rescale_factor = target_displacement_mm / applied
            load_at_target = abs(float(reaction_force)) * rescale_factor
            computed_stiffness = abs(float(reaction_force)) / applied

    common = {
        "model_file": str(model_file),
        "analysis_file": str(analysis_file),
        "analysis_var": analysis_var,
        "applied_displacement": applied_displacement if applied_displacement is not None else "",
        "reaction_force_N": reaction_force,
        "stiffness_N_per_mm": computed_stiffness if computed_stiffness != "" else stiffness,
    }
    pistoia_common = {
        "pistoia_file": str(pistoia_file),
        "pistoia_failure_var": failure_var,
        "pistoia_stiffness_var": stiffness_var,
        "pistoia_failure_load_N": pistoia_failure_load,
        "pistoia_failure_load_kN": pistoia_failure_load_magnitude / 1000.0
        if pistoia_failure_load_magnitude != ""
        else "",
        "pistoia_stiffness_N_per_mm": stiffness,
        "critical_volume_pct": pistoia_meta.get("critical_volume_pct", ""),
        "critical_ees": pistoia_meta.get("critical_ees", ""),
        "ees_at_crit_vol": pistoia_meta.get("ees_at_crit_vol", ""),
    }
    if profile == "femur":
        row = {
            "model_file": common["model_file"],
            "analysis_file": common["analysis_file"],
            "analysis_var": common["analysis_var"],
            "applied_displacement": common["applied_displacement"],
            "reaction_force_N": common["reaction_force_N"],
            "stiffness_N_per_mm": common["stiffness_N_per_mm"],
            "characteristic_length_mm": characteristic_length,
            **pistoia_common,
        }
    elif profile == "spine":
        row = {
            "model_file": common["model_file"],
            "analysis_file": common["analysis_file"],
            "analysis_var": common["analysis_var"],
            "applied_displacement": common["applied_displacement"],
            "reaction_force_N": common["reaction_force_N"],
            "stiffness_N_per_mm": common["stiffness_N_per_mm"],
            "characteristic_length_mm": characteristic_length,
            **pistoia_common,
        }
    else:
        row = {
            **common,
            "loads_csv": str(loads_csv),
            "pistoia_file": str(pistoia_file),
            "pistoia_csv": str(pistoia_csv),
            "pistoia_failure_var": failure_var,
            "pistoia_stiffness_var": stiffness_var,
            "target_displacement_percent": target_displacement,
            "rescale_factor_to_target": rescale_factor,
            "load_at_target_displacement_N": load_at_target,
            "load_at_target_displacement_kN": load_at_target / 1000.0 if load_at_target != "" else "",
            "pistoia_failure_load_N": pistoia_failure_load,
            "critical_volume_pct": pistoia_meta.get("critical_volume_pct", ""),
            "critical_ees": pistoia_meta.get("critical_ees", ""),
            "ees_at_crit_vol": pistoia_meta.get("ees_at_crit_vol", ""),
            "warnings": " | ".join(warnings or []),
        }

    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(str(output_file), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)
    return output_file


def run_faim_pipeline(
    *,
    model_file,
    output_prefix=None,
    analysis_var,
    pistoia_vars,
    failure_axis,
    threads=4,
    conda_env=None,
    conda_executable="conda",
    install_root=None,
    bin_dir=None,
    license_dir=None,
    faim_command=None,
    n88modelinfo_command=None,
    n88derivedfields_command=None,
    n88postfaim_command=None,
    n88pistoia_command=None,
    n88tabulate_command=None,
    n88copymodel_command=None,
    critical_volume=2.0,
    critical_strain=0.007,
    exclude=5000,
    run_pistoia=True,
    applied_displacement=None,
    target_displacement=0.2,
    report_profile="generic",
    solve_displacement_percent=None,
    compress=True,
    require_pistoia=False,
    dry_run=False,
):
    """Solve one N88 model and collect standard FAIM/Pistoia output files."""
    model_file = Path(model_file)
    output_prefix = Path(output_prefix) if output_prefix is not None else model_file.with_suffix("")
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    resolved_bin_dir = resolve_faim_bin_dir(install_root=install_root, bin_dir=bin_dir)
    env = build_faim_env(bin_dir=resolved_bin_dir, license_dir=license_dir)

    solver = resolve_faim_solver_command(
        install_root=install_root,
        bin_dir=resolved_bin_dir,
        explicit_path=faim_command,
    )
    modelinfo = resolve_faim_command(
        "n88modelinfo",
        install_root=install_root,
        bin_dir=resolved_bin_dir,
        explicit_path=n88modelinfo_command,
    )
    derivedfields = resolve_faim_command(
        "n88derivedfields",
        install_root=install_root,
        bin_dir=resolved_bin_dir,
        explicit_path=n88derivedfields_command,
    )
    postfaim = resolve_faim_command(
        "n88postfaim",
        install_root=install_root,
        bin_dir=resolved_bin_dir,
        explicit_path=n88postfaim_command,
    )
    pistoia = resolve_faim_command(
        "n88pistoia",
        install_root=install_root,
        bin_dir=resolved_bin_dir,
        explicit_path=n88pistoia_command,
    )
    tabulate = resolve_faim_command(
        "n88tabulate",
        install_root=install_root,
        bin_dir=resolved_bin_dir,
        explicit_path=n88tabulate_command,
    )
    copymodel = resolve_faim_command(
        "n88copymodel",
        install_root=install_root,
        bin_dir=resolved_bin_dir,
        explicit_path=n88copymodel_command,
    )

    analysis_file = output_prefix.with_name(output_prefix.name + "_analysis.txt")
    pistoia_file = output_prefix.with_name(output_prefix.name + "_pistoia.txt")
    loads_csv = output_prefix.with_name(output_prefix.name + "_loads.csv")
    pistoia_csv = output_prefix.with_name(output_prefix.name + "_pistoia.csv")
    results_csv = output_prefix.with_name(output_prefix.name + "_results.csv")
    log_file = output_prefix.with_name(output_prefix.name + ".log")

    warnings = []
    # Check whether the model already contains a solution so re-runs are cheap.
    modelinfo_text = run_command(
        [modelinfo, "--solutions", str(model_file)],
        conda_env=conda_env,
        conda_executable=conda_executable,
        env=env,
        dry_run=dry_run,
        capture_output=True,
        check=False,
    )

    if solve_displacement_percent is not None and not dry_run:
        if "Name :" in modelinfo_text:
            warnings.append(
                "Existing solution found; prescribed displacement was not updated before solve."
            )
            current_displacement = read_prescribed_displacement(model_file, report_profile)
            if current_displacement != "":
                applied_displacement = current_displacement
        else:
            current_displacement = read_prescribed_displacement(model_file, report_profile)
            try:
                sign_source = (
                    current_displacement
                    if current_displacement != ""
                    else applied_displacement
                    if applied_displacement is not None
                    else -1.0
                )
                displacement_sign = -1.0 if float(sign_source) < 0 else 1.0
            except (TypeError, ValueError):
                displacement_sign = -1.0
            update = set_prescribed_displacement_from_percent(
                model_file,
                report_profile=report_profile,
                failure_axis=failure_axis,
                target_displacement_percent=solve_displacement_percent,
                displacement_sign=displacement_sign,
            )
            applied_displacement = update["applied_displacement_mm"]

    if "Name :" not in modelinfo_text:
        run_command(
            [solver, "--engine=mt", "--threads={}".format(int(threads)), str(model_file)],
            conda_env=conda_env,
            conda_executable=conda_executable,
            env=env,
            dry_run=dry_run,
        )

    # Derived fields add quantities such as strain energy density to the model.
    run_command(
        [derivedfields, str(model_file)],
        conda_env=conda_env,
        conda_executable=conda_executable,
        env=env,
        dry_run=dry_run,
    )
    run_command(
        [postfaim, "--output_file", str(analysis_file), str(model_file)],
        conda_env=conda_env,
        conda_executable=conda_executable,
        env=env,
        dry_run=dry_run,
    )

    # Pistoia is useful but can fail for nonlinear or unusual models; make that
    # behavior configurable so the main solve result can still be collected.
    if run_pistoia:
        try:
            run_command(
                [
                    pistoia,
                    str(model_file),
                    "--critical_volume",
                    str(critical_volume),
                    "--critical_strain",
                    str(critical_strain),
                    "--exclude",
                    str(exclude),
                    "--output_file",
                    str(pistoia_file),
                ],
                conda_env=conda_env,
                conda_executable=conda_executable,
                env=env,
                dry_run=dry_run,
            )
        except subprocess.CalledProcessError:
            if require_pistoia:
                raise
            warnings.append("Pistoia postprocessing failed.")

    # Tabulate the reaction force and Pistoia values into small CSV files.
    run_command(
        [tabulate, "-H", "-d", ",", "-V", str(analysis_var), str(analysis_file), "-o", str(loads_csv)],
        conda_env=conda_env,
        conda_executable=conda_executable,
        env=env,
        dry_run=dry_run,
    )

    if run_pistoia and pistoia_vars:
        run_command(
            [
                tabulate,
                "-H",
                "-d",
                ",",
                "-V",
                ",".join(str(v) for v in pistoia_vars),
                str(pistoia_file),
                "-o",
                str(pistoia_csv),
            ],
            conda_env=conda_env,
            conda_executable=conda_executable,
            env=env,
            dry_run=dry_run,
        )

    if compress:
        run_command(
            [copymodel, "--compress", str(model_file), str(model_file)],
            conda_env=conda_env,
            conda_executable=conda_executable,
            env=env,
            dry_run=dry_run,
        )

    log_text = run_command(
        [modelinfo, "--log", str(model_file)],
        conda_env=conda_env,
        conda_executable=conda_executable,
        env=env,
        dry_run=dry_run,
        capture_output=True,
        check=False,
    )
    if not dry_run:
        log_file.write_text(log_text)
        write_results_csv(
            output_file=results_csv,
            model_file=model_file,
            analysis_file=analysis_file,
            pistoia_file=pistoia_file,
            loads_csv=loads_csv,
            pistoia_csv=pistoia_csv,
            analysis_var=analysis_var,
            pistoia_vars=pistoia_vars,
            failure_axis=failure_axis,
            applied_displacement=applied_displacement,
            target_displacement=target_displacement,
            report_profile=report_profile,
            warnings=warnings,
        )
        if str(report_profile or "").strip().lower() in ("spine", "femur"):
            for intermediate_file in (loads_csv, pistoia_csv):
                try:
                    Path(intermediate_file).unlink()
                except FileNotFoundError:
                    pass

    return {
        "analysis_file": analysis_file,
        "pistoia_file": pistoia_file,
        "loads_csv": loads_csv,
        "pistoia_csv": pistoia_csv,
        "results_csv": results_csv,
        "log_file": log_file,
        "warnings": warnings,
    }
