"""Generate finite-element models for spine vertebrae and hip femurs.

This module is intentionally a thin orchestration layer.  The scientific model
generation remains in ``ogo.cli.ref.SpineCompressionFe`` and
``ogo.cli.ref.SidewaysFallFe``; this command provides a clean user-facing entry
point for running one or many targets.
"""

import argparse
from collections import namedtuple
import importlib.util
from pathlib import Path
import sys
from typing import Callable, List, Optional, Sequence


SpineTarget = namedtuple("SpineTarget", ["level", "body_label", "process_label"])

SPINE_PRESETS = {
    "none": [],
    "benchmark-linear": [
        "--fe_displacement",
        "-0.2",
        "--elastic_E_func",
        "kopperdahl_trab_E",
        "--cort_elastic_E_func",
        "kopperdahl_trab_E",
    ],
    "benchmark-nonlinear": [
        "--fe_displacement",
        "-2.0",
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
    """Predict the spine model path written by ``OgoSpineCompressionFe``."""
    output_dir = output_path if output_path is not None else calibrated_image.parent
    return output_dir / "{}_vertebra_{}_{}.n88model".format(
        remove_extension(calibrated_image),
        target.body_label,
        target.level,
    )


def expected_femur_model_path(calibrated_image: Path, output_path: Optional[Path], side: str) -> Path:
    """Predict the femur model path written by ``OgoSidewaysFallFe``."""
    output_dir = output_path if output_path is not None else calibrated_image.parent
    suffix = "_LT_FEMUR_SF" if side == "left" else "_RT_FEMUR_SF"
    return output_dir / "{}{}.n88model".format(remove_extension(calibrated_image), suffix)


def option_value(argv: Sequence[str], option: str, default: Optional[str] = None) -> Optional[str]:
    """Return the final value for an option from an argv-style list."""
    value = default
    for index, token in enumerate(argv):
        if token == option and index + 1 < len(argv):
            value = argv[index + 1]
    return value


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

    # Spine compression uses z reaction force; sideways-fall hip uses y.
    if model_type == "spine":
        analysis_var = "fz_ns1"
        pistoia_vars = ["pis_fz_fail", "pis_stiffz"]
        failure_axis = "z"
    else:
        analysis_var = "fy_ns1"
        pistoia_vars = ["pis_fy_fail", "pis_stiffy"]
        failure_axis = "y"

    applied_displacement = option_value(generator_argv, "--fe_displacement", "-1.0")

    run_faim_pipeline(
        model_file=model_path,
        output_prefix=model_path.with_suffix(""),
        analysis_var=analysis_var,
        pistoia_vars=pistoia_vars,
        failure_axis=failure_axis,
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
        critical_volume=args.critical_volume,
        critical_strain=args.critical_strain,
        exclude=args.exclude,
        applied_displacement=applied_displacement,
        target_displacement=args.target_displacement,
        compress=not args.no_compress,
        require_pistoia=args.require_pistoia,
        dry_run=args.dry_run,
    )


def _call_cli(main_func: Callable[[], None], program: str, argv: Sequence[str]) -> None:
    old_argv = sys.argv
    try:
        sys.argv = [program, *argv]
        main_func()
    finally:
        sys.argv = old_argv


def run_spine_command(argv: Sequence[str]) -> None:
    from ogo.cli.ref.SpineCompressionFe import main as spine_compression_main

    _call_cli(spine_compression_main, "OgoSpineCompressionFe", argv)


def run_femur_command(argv: Sequence[str]) -> None:
    from ogo.cli.ref.SidewaysFallFe import main as sideways_fall_main

    _call_cli(sideways_fall_main, "OgoSidewaysFallFe", argv)


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
    parser.add_argument("--threads", type=int, default=4, help="FAIM solver thread count.")
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
        default=2.0,
        help="Pistoia critical volume percentage.",
    )
    parser.add_argument(
        "--critical_strain",
        type=float,
        default=0.007,
        help="Pistoia critical EES strain.",
    )
    parser.add_argument("--exclude", type=int, default=5000, help="Material ID excluded by Pistoia.")
    parser.add_argument(
        "--target_displacement",
        type=float,
        default=0.2,
        help="Standard displacement threshold used for the rescaled load endpoint.",
    )
    parser.add_argument(
        "--require_pistoia",
        action="store_true",
        help="Fail if Pistoia postprocessing fails.",
    )
    parser.add_argument(
        "--no_compress",
        action="store_true",
        help="Skip n88copymodel --compress after solving.",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ogoGenerateFEM",
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


def _run_legacy(argv: Sequence[str]) -> None:
    """Keep the historical ``--model_type`` dispatch working."""
    parser = argparse.ArgumentParser(prog="ogoGenerateFEM")
    parser.add_argument("--model_type", choices=["vertebra", "femur"], required=True)
    args, remaining = parser.parse_known_args(argv)
    if args.model_type == "vertebra":
        run_spine_command(remaining)
    else:
        run_femur_command(remaining)


def main(argv: Optional[Sequence[str]] = None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--model_type" in argv:
        _run_legacy(argv)
        return

    parser = build_parser()
    args, extra_args = parser.parse_known_args(argv)

    if args.model_type is None:
        parser.error("model type is required: choose spine or hip")

    if args.model_type == "spine":
        spine_extra_args = spine_preset_args(args.preset) + list(extra_args)
        for target in args.vertebra:
            cmd = build_spine_command(
                calibrated_image=args.calibrated_image,
                bone_mask=args.bone_mask,
                target=target,
                output_path=args.output_path,
                extra_args=spine_extra_args,
            )
            if args.dry_run:
                print_dry_run("OgoSpineCompressionFe", cmd)
            else:
                run_spine_command(cmd)
            if not args.no_solve:
                solve_model(
                    expected_spine_model_path(args.calibrated_image, args.output_path, target),
                    args,
                    "spine",
                    cmd,
                )
        return

    if args.model_type == "hip":
        for side in expand_sides(args.side or ["both"]):
            cmd = build_femur_command(
                calibrated_image=args.calibrated_image,
                bone_mask=args.bone_mask,
                side=side,
                output_path=args.output_path,
                extra_args=extra_args,
            )
            if args.dry_run:
                print_dry_run("OgoSidewaysFallFe", cmd)
            else:
                run_femur_command(cmd)
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
