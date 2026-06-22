"""Small FAIM/N88 command adapter used by Ogo FE CLIs."""

import csv
import os
from pathlib import Path
import re
import subprocess
from typing import Dict, Iterable, List, Optional, Sequence


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
    warnings=None,
):
    """Write a compact, user-facing CSV with loads, stiffness, and Pistoia values."""
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

    # The benchmark endpoint is reported at 0.2. If the model was solved at a
    # different prescribed displacement, scale the linear reaction force back to
    # that standard point.
    rescale_factor = ""
    load_at_target = ""
    try:
        applied = abs(float(applied_displacement))
        target = abs(float(target_displacement))
        if applied > 0 and reaction_force != "":
            rescale_factor = target / applied
            load_at_target = abs(float(reaction_force)) * rescale_factor
    except (TypeError, ValueError):
        pass

    row = {
        "model_file": str(model_file),
        "analysis_file": str(analysis_file),
        "pistoia_file": str(pistoia_file),
        "loads_csv": str(loads_csv),
        "pistoia_csv": str(pistoia_csv),
        "analysis_var": analysis_var,
        "pistoia_failure_var": failure_var,
        "pistoia_stiffness_var": stiffness_var,
        "applied_displacement": applied_displacement if applied_displacement is not None else "",
        "target_displacement_percent": target_displacement,
        "rescale_factor_to_target": rescale_factor,
        "load_at_0p2_percent_N": load_at_target,
        "load_at_0p2_percent_kN": load_at_target / 1000.0 if load_at_target != "" else "",
        "reaction_force_N": reaction_force,
        "pistoia_failure_load_N": pistoia_failure_load,
        "stiffness_N_per_mm": stiffness,
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
    applied_displacement=None,
    target_displacement=0.2,
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

    if pistoia_vars:
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
            warnings=warnings,
        )

    return {
        "analysis_file": analysis_file,
        "pistoia_file": pistoia_file,
        "loads_csv": loads_csv,
        "pistoia_csv": pistoia_csv,
        "results_csv": results_csv,
        "log_file": log_file,
        "warnings": warnings,
    }
