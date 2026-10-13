"""Generate and solve finite-element models for spine vertebrae and hip femurs.

This module is intentionally a thin orchestration layer. The scientific model
generation lives in ``ogo.fea``; this command provides the maintained
user-facing entry point for running one or many targets.
"""

import argparse
from collections import namedtuple
import importlib.util
from pathlib import Path
import sys
from typing import Callable, List, Optional, Sequence

from ogo.fea.spine import (
    BENCHMARK_LINEAR_FE_DISPLACEMENT_MM,
    BENCHMARK_NONLINEAR_FE_DISPLACEMENT_MM,
    DEFAULT_SPINE_REGISTRATION_BACKEND,
    DEFAULT_SPINE_REGISTRATION_ITERATIONS,
    DEFAULT_SPINE_REGISTRATION_LANDMARKS,
    DEFAULT_SPINE_ICP_TARGET,
    SPINE_ICP_TARGETS,
    solve_report_profile as spine_solve_report_profile,
)
from ogo.fea.femur import (
    solve_report_profile as femur_solve_report_profile,
)


SpineTarget = namedtuple("SpineTarget", ["level", "body_label", "process_label"])

SPINE_PRESETS = {
    "none": [],
    "benchmark-linear": [
        "--fe_displacement",
        str(BENCHMARK_LINEAR_FE_DISPLACEMENT_MM),
        "--elastic_E_func",
        "kopperdahl_trab_E",
        "--cort_elastic_E_func",
        "kopperdahl_trab_E",
    ],
    "benchmark-nonlinear": [
        "--fe_displacement",
        str(BENCHMARK_NONLINEAR_FE_DISPLACEMENT_MM),
        "--pmma_yield_compression",
        "70.0",
        "--pmma_yield_tension",
        "70.0",
        "--elastic_E_func",
        "kopperdahl_trab_E",
        "--yield_comp_func",
        "kopperdahl_trab_yc",
        "--yield_tens_func",
        "kopperdahl_trab_yc",
        "--cort_elastic_E_func",
        "kopperdahl_trab_E",
        "--cort_yield_comp_func",
        "kopperdahl_trab_yc",
        "--cort_yield_tens_func",
        "kopperdahl_trab_yc",
    ],
}


def remove_extension(path: Path) -> str:
    """Remove common medical-image extensions while keeping the useful stem."""
    name = path.name
    if name.endswith(".nii.gz"):
        return name[:-7]
    return path.stem


def parse_spine_target(value: str) -> SpineTarget:
    """Parse ``LEVEL:BODY_LABEL:PROCESS_LABEL`` target syntax."""
    parts = [part.strip() for part in value.replace(",", ":").split(":")]
    if len(parts) != 3 or not all(parts):
        raise argparse.ArgumentTypeError(
            "spine targets must use LEVEL:BODY_LABEL:PROCESS_LABEL, for example L1:2:1"
        )
    level, body_label, process_label = parts
    try:
        return SpineTarget(level=level, body_label=int(body_label), process_label=int(process_label))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("body and process labels must be integers") from exc


def expand_sides(sides: Sequence[str]) -> List[str]:
    """Expand the friendly hip side syntax into concrete left/right jobs."""
    expanded = []  # type: List[str]
    for side in sides:
        if side == "both":
            expanded.extend(["left", "right"])
        elif side not in expanded:
            expanded.append(side)
    return expanded


def spine_preset_args(name: str) -> List[str]:
    """Return a copy of the lower-level arguments for one spine preset."""
    return list(SPINE_PRESETS[name])


def spine_registration_args(args: argparse.Namespace) -> List[str]:
    """Return lower-level spine registration options owned by the wrapper."""
    return [
        "--registration_backend",
        str(args.registration_backend),
        "--registration_landmarks",
        str(args.registration_landmarks),
        "--registration_iterations",
        str(args.registration_iterations),
        "--spine_icp_target",
        str(args.spine_icp_target),
    ]


def build_spine_command(
    *,
    calibrated_image: Path,
    bone_mask: Path,
    target: SpineTarget,
    output_path: Optional[Path],
    extra_args: Sequence[str],
) -> List[str]:
    """Build the lower-level spine FE generator command for one vertebra."""
    cmd = [
        str(calibrated_image),
        str(bone_mask),
        "--mask_threshold",
        str(target.body_label),
        "--process_mask_threshold",
        str(target.process_label),
        "--appendix",
        target.level,
    ]
    if output_path is not None:
        cmd.extend(["--output_path", str(output_path)])
    cmd.extend(extra_args)
    return cmd


def build_femur_command(
    *,
    calibrated_image: Path,
    bone_mask: Path,
    side: str,
    output_path: Optional[Path],
    extra_args: Sequence[str],
) -> List[str]:
    """Build the lower-level femur FE generator command for one side."""
    femur_side = {"left": "1", "right": "2"}[side]
    cmd = [
        str(calibrated_image),
        str(bone_mask),
        "--femur_side",
        femur_side,
    ]
    if output_path is not None:
        cmd.extend(["--output_path", str(output_path)])
    cmd.extend(extra_args)
    return cmd


def expected_spine_model_path(calibrated_image: Path, output_path: Optional[Path], target: SpineTarget) -> Path:
    """Predict the spine model path written by the spine builder."""
    output_dir = output_path if output_path is not None else calibrated_image.parent
    return output_dir / "{}_{}.n88model".format(remove_extension(calibrated_image), target.level)


def expected_femur_model_path(calibrated_image: Path, output_path: Optional[Path], side: str) -> Path:
    """Predict the femur model path written by the femur builder."""
    output_dir = output_path if output_path is not None else calibrated_image.parent
    stem = "LF" if side == "left" else "RF"
    return output_dir / "{}_{}.n88model".format(remove_extension(calibrated_image), stem)


def pistoia_mask_output_path(model_path: Path) -> Path:
    """Return the model-space Pistoia ROI mask sidecar path for a model."""
    model_path = Path(model_path)
    return model_path.with_name(model_path.with_suffix("").name + "_pistoia_mask.nii.gz")


def pistoia_mask_source_path(args: argparse.Namespace) -> Optional[Path]:
    """Return the input-space source image used for masked Pistoia, if any."""
    if getattr(args, "pistoia_mask", None) is not None:
        return Path(args.pistoia_mask)
    if getattr(args, "pistoia_mask_label", None):
        return Path(args.bone_mask)
    return None


def model_space_pistoia_mask_path(model_path: Path, args: argparse.Namespace) -> Optional[Path]:
    """Return the transformed Pistoia ROI mask sidecar expected by the solver."""
    if pistoia_mask_source_path(args) is None:
        return None
    return pistoia_mask_output_path(model_path)


def ensure_output_directory(output_path: Optional[Path]) -> None:
    """Create an explicit model output directory before generation starts."""
    if output_path is not None:
        output_path.mkdir(parents=True, exist_ok=True)


def option_value(argv: Sequence[str], option: str, default: Optional[str] = None) -> Optional[str]:
    """Return the final value for an option from an argv-style list."""
    value = default
    for index, token in enumerate(argv):
        prefix = option + "="
        if token.startswith(prefix):
            value = token[len(prefix) :]
        if token == option and index + 1 < len(argv):
            value = argv[index + 1]
    return value


def option_float(argv: Sequence[str], option: str, default: float) -> float:
    """Return the final numeric value for an argv-style option."""
    return float(option_value(argv, option, str(default)))


def option_optional_float(argv: Sequence[str], option: str) -> Optional[float]:
    """Return the final numeric value for an optional argv-style option."""
    value = option_value(argv, option)
    return None if value is None else float(value)


def option_int(argv: Sequence[str], option: str, default: int) -> int:
    """Return the final integer value for an argv-style option."""
    return int(option_value(argv, option, str(default)))


def option_present(argv: Sequence[str], option: str) -> bool:
    """Return whether a flag-like option is present."""
    return option in argv


def option_path(argv: Sequence[str], option: str, default: Optional[Path] = None) -> Optional[str]:
    """Return the final path value for an option as a string."""
    value = option_value(argv, option, str(default) if default is not None else None)
    return None if value is None else str(value)


def option_n_values(argv: Sequence[str], option: str, count: int, default: Sequence) -> list:
    """Return a fixed-width option value list from an argv-style list."""
    values = list(default)
    for index, token in enumerate(argv):
        if token == option and index + count < len(argv):
            values = list(argv[index + 1 : index + 1 + count])
    parsed = []
    for value in values:
        token = str(value).strip().lower() if value is not None else "none"
        parsed.append(None if token in {"", "none", "null", "auto"} else value)
    return parsed


def target_displacement_percent(args: argparse.Namespace, model_type: str) -> float:
    """Return the maintained endpoint displacement percentage for this model."""
    if args.target_displacement is not None:
        return args.target_displacement
    return solve_report_profile(args, model_type)["target_displacement_percent"]


def solve_report_profile(args: argparse.Namespace, model_type: str) -> dict:
    """Return site-specific solve/report settings with CLI overrides applied."""
    if model_type == "spine":
        profile = spine_solve_report_profile(preset=getattr(args, "preset", None))
    else:
        profile = femur_solve_report_profile()
    profile = profile.copy()
    if args.target_displacement is not None:
        profile["target_displacement_percent"] = args.target_displacement
    return profile


def critical_volume_percent(args: argparse.Namespace, model_type: str) -> float:
    """Return the Pistoia critical-volume setting for this model."""
    if args.critical_volume is not None:
        return args.critical_volume
    return float(solve_report_profile(args, model_type).get("critical_volume", 2.0))


def critical_strain(args: argparse.Namespace, model_type: str) -> float:
    """Return the Pistoia critical-strain setting for this model."""
    if args.critical_strain is not None:
        return args.critical_strain
    return float(solve_report_profile(args, model_type).get("critical_strain", 0.007))




def solve_model(
    model_path: Path,
    args: argparse.Namespace,
    model_type: str,
    generator_argv: Sequence[str],
) -> None:
    """Run the generic FAIM adapter on one generated model."""
    try:
        from ogo.util.faim import run_faim_pipeline
    except ModuleNotFoundError:
        # Support direct script use: python ogo/cli/GenerateFEM.py ...
        module_path = Path(__file__).resolve().parents[1] / "util" / "faim.py"
        spec = importlib.util.spec_from_file_location("ogo_util_faim", str(module_path))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        run_faim_pipeline = module.run_faim_pipeline

    profile = solve_report_profile(args, model_type)
    default_applied_displacement = str(profile["default_applied_displacement"])
    applied_displacement = option_value(
        generator_argv,
        "--fe_displacement",
        default_applied_displacement,
    )

    solve_displacement_percent = None if args.use_absolute_fe_displacement else profile["target_displacement_percent"]
    target_displacement = abs(float(applied_displacement)) if args.use_absolute_fe_displacement else profile["target_displacement_percent"]
    run_pistoia = args.run_pistoia or args.require_pistoia or pistoia_mask_source_path(args) is not None

    run_faim_pipeline(
        model_file=model_path,
        output_prefix=model_path.with_suffix(""),
        analysis_var=profile["analysis_var"],
        pistoia_vars=profile["pistoia_vars"],
        failure_axis=profile["failure_axis"],
        threads=args.threads,
        conda_env=args.faim_env,
        conda_executable=args.conda_executable,
        install_root=args.faim_install_root,
        bin_dir=args.faim_bin_dir,
        license_dir=args.faim_license_dir,
        faim_command=args.faim_command,
        n88modelinfo_command=args.n88modelinfo_command,
        n88derivedfields_command=args.n88derivedfields_command,
        n88postfaim_command=args.n88postfaim_command,
        n88pistoia_command=args.n88pistoia_command,
        n88tabulate_command=args.n88tabulate_command,
        n88copymodel_command=args.n88copymodel_command,
        critical_volume=critical_volume_percent(args, model_type),
        masked_critical_volume=getattr(args, "masked_critical_volume", None),
        critical_strain=critical_strain(args, model_type),
        exclude=args.exclude,
        run_pistoia=run_pistoia,
        pistoia_mask_file=model_space_pistoia_mask_path(model_path, args),
        applied_displacement=applied_displacement,
        target_displacement=target_displacement,
        report_profile=profile["report_profile"],
        solve_displacement_percent=solve_displacement_percent,
        compress=not args.no_compress,
        require_pistoia=args.require_pistoia,
        dry_run=args.dry_run,
    )


def write_modeling_metadata(
    model_path: Path,
    model_type: str,
    generator_argv: Sequence[str],
    args: argparse.Namespace,
    bc_audit_summary: Optional[dict] = None,
) -> Optional[Path]:
    """Record resolved construction settings and the requested solve profile."""
    from ogo.fea.metadata import write_modeling_metadata as write_record

    reporting = {
        "critical_volume_pct": critical_volume_percent(args, model_type),
        "critical_strain": critical_strain(args, model_type),
        "target_displacement_percent": target_displacement_percent(args, model_type),
    }
    return write_record(model_path, model_type, generator_argv, args,
                        bc_audit_summary, reporting=reporting)


def audit_generated_model(model_path: Path, model_type: str, args: argparse.Namespace) -> Optional[dict]:
    """Return BC audit summary and optionally write a debug PNG."""
    if args.skip_bc_audit or args.dry_run:
        return None
    if not model_path.exists():
        return None

    from ogo.cli.CheckFEModelBC import audit_model

    audit_kind = "spine-compression" if model_type == "spine" else "femur-sideways"
    result = audit_model(
        model_path,
        model=audit_kind,
        flat_tolerance=args.bc_audit_flat_tolerance,
        write_json=False,
        write_csv_file=False,
        write_plot=args.debug,
    )
    if result["png_path"] is not None:
        print(f"Wrote {result['png_path']}")
    if not result["passed"]:
        failed = [check["name"] for check in result["checks"] if not check["passed"]]
        raise RuntimeError("BC audit failed for {}: {}".format(model_path, "; ".join(failed)))
    return result["summary"]


def _call_cli(main_func: Callable, argv: Sequence[str]) -> None:
    """Call an anatomy builder without changing process-global arguments."""
    try:
        main_func(list(argv))
    except SystemExit as exc:
        if exc.code not in (None, 0):
            raise


def run_spine_command(argv: Sequence[str]) -> None:
    from ogo.fea.spine import main as spine_main

    _call_cli(spine_main, argv)


def run_femur_command(argv: Sequence[str]) -> None:
    from ogo.fea.femur import main as femur_main

    _call_cli(femur_main, argv)


def print_dry_run(program: str, argv: Sequence[str]) -> None:
    print(" ".join([program, *argv]))


def _add_common_image_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("calibrated_image", type=Path, help="Calibrated density image.")
    parser.add_argument("bone_mask", type=Path, help="Bone mask or labelled bone mask.")
    parser.add_argument(
        "--output_path",
        type=Path,
        default=None,
        help="Directory for generated .n88model files. Defaults to the input image directory.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print generated lower-level commands without running model generation or solving.",
    )
    parser.add_argument(
        "--no-solve",
        action="store_true",
        help="Only generate the .n88model; skip FAIM solve and postprocessing.",
    )
    parser.add_argument("--threads", type=int, default=4, help="Per-job VTK/ITK generation and FAIM solver thread limit.")
    parser.add_argument("--faim_env", default=None, help="Optional conda environment for FAIM/N88 tools.")
    parser.add_argument("--conda_executable", default="conda", help="Conda executable for --faim_env.")
    parser.add_argument("--faim_install_root", default=None, help="Optional FAIM install root.")
    parser.add_argument("--faim_bin_dir", default=None, help="Optional directory containing FAIM/N88 commands.")
    parser.add_argument("--faim_license_dir", default=None, help="Optional Numerics88 license directory.")
    parser.add_argument("--faim_command", default=None, help="Override FAIM solver command.")
    parser.add_argument("--n88modelinfo_command", default=None, help="Override n88modelinfo command.")
    parser.add_argument("--n88derivedfields_command", default=None, help="Override n88derivedfields command.")
    parser.add_argument("--n88postfaim_command", default=None, help="Override n88postfaim command.")
    parser.add_argument("--n88pistoia_command", default=None, help="Override n88pistoia command.")
    parser.add_argument("--n88tabulate_command", default=None, help="Override n88tabulate command.")
    parser.add_argument("--n88copymodel_command", default=None, help="Override n88copymodel command.")
    parser.add_argument(
        "--critical_volume",
        type=float,
        default=None,
        help="Pistoia critical volume percentage. Defaults are model-profile specific.",
    )
    parser.add_argument(
        "--masked_critical_volume",
        type=float,
        default=None,
        help="Regional Pistoia critical volume percentage; defaults to --critical_volume.",
    )
    parser.add_argument(
        "--critical_strain",
        type=float,
        default=None,
        help="Pistoia critical EES strain. Defaults are model-profile specific.",
    )
    parser.add_argument("--exclude", type=int, default=5000, help="Material ID excluded by Pistoia.")
    parser.add_argument(
        "--target_displacement",
        type=float,
        default=None,
        help=(
            "Profile reporting endpoint as percent strain. Spine presets default to "
            "0.68 for benchmark-linear and 4.0 for benchmark-nonlinear; hip defaults "
            "to 4.0. Values are converted to mm from generated model geometry."
        ),
    )
    parser.add_argument(
        "--use_absolute_fe_displacement",
        action="store_true",
        help=(
            "Solve at the explicit --fe_displacement value in mm instead of "
            "converting the profile target displacement percent from model geometry."
        ),
    )
    parser.add_argument(
        "--require_pistoia",
        action="store_true",
        help="Fail if Pistoia postprocessing fails.",
    )
    parser.add_argument(
        "--run_pistoia",
        action="store_true",
        help="Run Pistoia postprocessing and include Pistoia metrics when available.",
    )
    parser.add_argument(
        "--pistoia_mask",
        type=Path,
        default=None,
        help=(
            "Optional ROI mask in input image space for an additional masked Pistoia "
            "calculation. The model builder writes the transformed model-space mask "
            "next to the .n88model."
        ),
    )
    parser.add_argument(
        "--pistoia_mask_label",
        type=int,
        action="append",
        default=None,
        help=(
            "Label to keep from the Pistoia ROI source before masked Pistoia. Repeat "
            "for a multi-label ROI. If labels are supplied without --pistoia_mask, "
            "the main bone mask is used as the ROI source."
        ),
    )
    parser.add_argument(
        "--no_compress",
        action="store_true",
        help="Skip n88copymodel --compress after solving.",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Write debug sidecars such as BC audit PNGs and spine QC quick looks.",
    )
    parser.add_argument(
        "--skip_bc_audit",
        action="store_true",
        help="Skip automatic boundary-condition audit files after model generation.",
    )
    parser.add_argument(
        "--bc_audit_flat_tolerance",
        type=float,
        default=1.0e-4,
        help="Flat-plane tolerance for automatic boundary-condition audit checks.",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ogoFEA",
        description=(
            "Generate N88 finite-element models for spine vertebrae and hip femurs. "
            "Unknown options are forwarded to the selected lower-level FE generator."
        ),
    )
    subparsers = parser.add_subparsers(dest="model_type")

    spine = subparsers.add_parser(
        "spine",
        help="Generate one or more vertebral compression FE models.",
        description=(
            "Generate spine FE models. Repeat --vertebra for each level to run. "
            "Each target uses LEVEL:BODY_LABEL:PROCESS_LABEL."
        ),
    )
    _add_common_image_args(spine)
    spine.add_argument(
        "--vertebra",
        action="append",
        type=parse_spine_target,
        required=True,
        metavar="LEVEL:BODY_LABEL:PROCESS_LABEL",
        help="Vertebral target, for example L1:2:1. Repeat for all levels to process.",
    )
    spine.add_argument(
        "--preset",
        choices=sorted(SPINE_PRESETS),
        default="benchmark-linear",
        help=(
            "Spine FE parameter preset. The default benchmark-linear preset uses the "
            "public spineFE-benchmark linear model settings; use none to pass only "
            "explicit lower-level options."
        ),
    )
    spine.add_argument(
        "--reference_path",
        type=Path,
        default=None,
        help=(
            "Optional spine ICP reference surface. If omitted, the bundled "
            "reference matching --spine_icp_target is used."
        ),
    )
    spine.add_argument(
        "--registration_backend",
        choices=("vtk", "numpy"),
        default=DEFAULT_SPINE_REGISTRATION_BACKEND,
        help=(
            "Spine reference alignment backend. vtk keeps the original VTK ICP "
            "solver; numpy uses the deterministic point-cloud helper. "
            "(default: %(default)s)"
        ),
    )
    spine.add_argument(
        "--registration_landmarks",
        type=int,
        default=DEFAULT_SPINE_REGISTRATION_LANDMARKS,
        help=(
            "Maximum spine ICP landmarks/sampled points. For vtk this maps to "
            "SetMaximumNumberOfLandmarks; for numpy this caps sampled surface "
            "points. (default: %(default)s)"
        ),
    )
    spine.add_argument(
        "--registration_iterations",
        type=int,
        default=DEFAULT_SPINE_REGISTRATION_ITERATIONS,
        help="Maximum spine ICP iterations. (default: %(default)s)",
    )
    spine.add_argument(
        "--spine_icp_target",
        choices=SPINE_ICP_TARGETS,
        default=DEFAULT_SPINE_ICP_TARGET,
        help=(
            "Segmentation surface used for spine ICP. body uses the vertebral body "
            "label only; vertebra uses body plus posterior-process labels and "
            "therefore requires a matching full-vertebra reference surface. "
            "(default: %(default)s)"
        ),
    )

    hip = subparsers.add_parser(
        "hip",
        help="Generate hip sideways-fall FE models for left, right, or both femurs.",
    )
    _add_common_image_args(hip)
    hip.add_argument(
        "--side",
        action="append",
        choices=["left", "right", "both"],
        default=None,
        help="Femur side to process. Repeatable; defaults to both.",
    )

    return parser


def limit_generation_threads(threads: int) -> None:
    """Keep VTK/ITK model generation within the requested per-job CPU budget."""
    if threads < 1:
        raise ValueError("threads must be positive.")
    import SimpleITK as sitk
    import vtk

    vtk.vtkMultiThreader.SetGlobalMaximumNumberOfThreads(threads)
    vtk.vtkMultiThreader.SetGlobalDefaultNumberOfThreads(threads)
    vtk.vtkSMPTools.Initialize(threads)
    sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(threads)


def main(argv: Optional[Sequence[str]] = None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)

    parser = build_parser()
    args, extra_args = parser.parse_known_args(argv)

    if args.model_type is None:
        parser.error("model type is required: choose spine or hip")

    if not args.dry_run:
        limit_generation_threads(args.threads)
        ensure_output_directory(args.output_path)

    pistoia_mask_source = pistoia_mask_source_path(args)
    pistoia_mask_args = []
    if pistoia_mask_source is not None:
        pistoia_mask_args.extend(["--pistoia_mask", str(pistoia_mask_source)])
    for label in args.pistoia_mask_label or []:
        pistoia_mask_args.extend(["--pistoia_mask_label", str(label)])

    if args.model_type == "spine":
        reference_args = []
        if args.reference_path is not None:
            reference_args.extend(["--reference_path", str(args.reference_path)])
        spine_extra_args = (
            spine_preset_args(args.preset)
            + spine_registration_args(args)
            + reference_args
            + pistoia_mask_args
            + list(extra_args)
            + [
            "--quality_control",
            str(bool(args.debug)),
        ])
        for target in args.vertebra:
            cmd = build_spine_command(
                calibrated_image=args.calibrated_image,
                bone_mask=args.bone_mask,
                target=target,
                output_path=args.output_path,
                extra_args=spine_extra_args,
            )
            if args.dry_run:
                print_dry_run("ogoFEA-spine-builder", cmd)
            else:
                run_spine_command(cmd)
                model_path = expected_spine_model_path(args.calibrated_image, args.output_path, target)
                bc_audit_summary = audit_generated_model(
                    model_path,
                    "spine",
                    args,
                )
                write_modeling_metadata(model_path, "spine", cmd, args, bc_audit_summary)
            if not args.no_solve:
                solve_model(
                    expected_spine_model_path(args.calibrated_image, args.output_path, target),
                    args,
                    "spine",
                    cmd,
                )
        return

    if args.model_type == "hip":
        femur_extra_args = pistoia_mask_args + list(extra_args)
        for side in expand_sides(args.side or ["both"]):
            cmd = build_femur_command(
                calibrated_image=args.calibrated_image,
                bone_mask=args.bone_mask,
                side=side,
                output_path=args.output_path,
                extra_args=femur_extra_args,
            )
            if args.dry_run:
                print_dry_run("ogoFEA-hip-builder", cmd)
            else:
                run_femur_command(cmd)
                model_path = expected_femur_model_path(args.calibrated_image, args.output_path, side)
                bc_audit_summary = audit_generated_model(
                    model_path,
                    "hip",
                    args,
                )
                write_modeling_metadata(model_path, "hip", cmd, args, bc_audit_summary)
            if not args.no_solve:
                solve_model(
                    expected_femur_model_path(args.calibrated_image, args.output_path, side),
                    args,
                    "hip",
                    cmd,
                )
        return

    parser.error(f"Unsupported model type: {args.model_type}")


if __name__ == "__main__":
    main()
