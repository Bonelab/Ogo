"""Audit FE model boundary-condition sets and displacement directions."""

import argparse
import csv
import json
from pathlib import Path
import subprocess

import numpy as np
from netCDF4 import Dataset

from ogo.fea.femur import (
    DISTAL_FEMUR_NODE_SET,
    FEMORAL_HEAD_NODE_SET,
    GREATER_TROCHANTER_NODE_SET,
    SIDEWAYS_FALL_NODE_SETS,
)


SENSE_TO_AXIS = {1: "x", 2: "y", 3: "z"}


def _node_coordinates(root, node_numbers):
    coords = root.groups["Parts"].groups["Part1"].variables["NodeCoordinates"][:]
    return coords[np.asarray(node_numbers, dtype=np.int64) - 1]


def _bounds(points):
    mins = points.min(axis=0)
    maxs = points.max(axis=0)
    return {
        "x_min": float(mins[0]),
        "x_max": float(maxs[0]),
        "y_min": float(mins[1]),
        "y_max": float(maxs[1]),
        "z_min": float(mins[2]),
        "z_max": float(maxs[2]),
    }


def _centroid(points):
    center = points.mean(axis=0)
    return {"x": float(center[0]), "y": float(center[1]), "z": float(center[2])}


def read_node_sets(root):
    out = {}
    for name, group in root.groups["Sets"].groups["NodeSets"].groups.items():
        node_numbers = group.variables["NodeNumber"][:]
        points = _node_coordinates(root, node_numbers)
        out[name] = {
            "count": int(len(node_numbers)),
            "bounds": _bounds(points),
            "centroid": _centroid(points),
        }
    return out


def read_constraints(root):
    out = {}
    for name, group in root.groups["Constraints"].groups.items():
        node_numbers = group.variables["NodeNumber"][:]
        senses = np.asarray(group.variables["Sense"][:])
        values = np.asarray(group.variables["Value"][:], dtype=float)
        points = _node_coordinates(root, node_numbers)
        unique_senses = sorted(set(int(value) for value in senses.tolist()))
        unique_values = sorted(set(float(value) for value in values.tolist()))
        out[name] = {
            "count": int(len(node_numbers)),
            "axes": [SENSE_TO_AXIS.get(value, str(value)) for value in unique_senses],
            "values": unique_values,
            "bounds": _bounds(points),
            "centroid": _centroid(points),
        }
    return out


def _single_axis(constraints, name):
    axes = constraints.get(name, {}).get("axes", [])
    return axes[0] if len(axes) == 1 else None


def _single_value(constraints, name):
    values = constraints.get(name, {}).get("values", [])
    return values[0] if len(values) == 1 else None


def audit_femur_sideways(node_sets, constraints, tolerance):
    checks = []

    for name in SIDEWAYS_FALL_NODE_SETS:
        count = node_sets.get(name, {}).get("count", 0)
        checks.append({"name": f"{name} exists", "passed": count > 0, "detail": str(count)})

    fh = node_sets.get(FEMORAL_HEAD_NODE_SET)
    gt = node_sets.get(GREATER_TROCHANTER_NODE_SET)
    if fh and gt:
        y_gap = gt["centroid"]["y"] - fh["centroid"]["y"]
        applied = _single_value(constraints, "top_displacement")
        axis = _single_axis(constraints, "top_displacement")
        compressive = (
            axis == "y"
            and applied is not None
            and applied != 0
            and y_gap != 0
            and np.sign(applied) == np.sign(y_gap)
        )
        checks.append({
            "name": "femoral head displacement moves toward fixed trochanter",
            "passed": bool(compressive),
            "detail": f"axis={axis}, value={applied}, gt_y_minus_fh_y={y_gap:.6g}",
        })

    distal = node_sets.get(DISTAL_FEMUR_NODE_SET)
    if distal:
        z_span = distal["bounds"]["z_max"] - distal["bounds"]["z_min"]
        checks.append({
            "name": "distal shaft set is one flat z plane",
            "passed": z_span <= tolerance,
            "detail": f"z_span={z_span:.6g}",
        })

    expected_zero_constraints = {
        "bottom_fixed_x": "x",
        "bottom_fixed_z": "z",
        "bottom_fixed_y_PMMA": "y",
    }
    for name, axis in expected_zero_constraints.items():
        actual_axis = _single_axis(constraints, name)
        value = _single_value(constraints, name)
        checks.append({
            "name": f"{name} fixes {axis}",
            "passed": actual_axis == axis and value == 0.0,
            "detail": f"axis={actual_axis}, value={value}",
        })

    return checks


def audit_spine_compression(node_sets, constraints, tolerance):
    checks = []

    required_sets = ["body_top", "body_bottom"]
    for name in required_sets:
        count = node_sets.get(name, {}).get("count", 0)
        checks.append({"name": f"{name} exists", "passed": count > 0, "detail": str(count)})

    top = node_sets.get("body_top")
    bottom = node_sets.get("body_bottom")
    if top and bottom:
        z_gap = bottom["centroid"]["z"] - top["centroid"]["z"]
        applied = _single_value(constraints, "top_displacement")
        axis = _single_axis(constraints, "top_displacement")
        compressive = (
            axis == "z"
            and applied is not None
            and applied != 0
            and z_gap != 0
            and np.sign(applied) == np.sign(z_gap)
        )
        checks.append({
            "name": "top displacement moves toward bottom support",
            "passed": bool(compressive),
            "detail": f"axis={axis}, value={applied}, bottom_z_minus_top_z={z_gap:.6g}",
        })

        for set_name, set_data in (("body_top", top), ("body_bottom", bottom)):
            z_span = set_data["bounds"]["z_max"] - set_data["bounds"]["z_min"]
            checks.append({
                "name": f"{set_name} set is one flat z plane",
                "passed": z_span <= tolerance,
                "detail": f"z_span={z_span:.6g}",
            })

    axis = _single_axis(constraints, "bottom_fixed_z")
    value = _single_value(constraints, "bottom_fixed_z")
    checks.append({
        "name": "bottom_fixed_z fixes z",
        "passed": axis == "z" and value == 0.0,
        "detail": f"axis={axis}, value={value}",
    })
    return checks


def audit_checks(model, node_sets, constraints, tolerance):
    if model == "femur-sideways":
        return audit_femur_sideways(node_sets, constraints, tolerance)
    if model == "spine-compression":
        return audit_spine_compression(node_sets, constraints, tolerance)
    raise ValueError(f"Unsupported audit model: {model}")


def write_csv(path, node_sets, constraints, checks):
    rows = []
    for kind, source in (("node_set", node_sets), ("constraint", constraints)):
        for name, data in source.items():
            bounds = data["bounds"]
            center = data["centroid"]
            rows.append({
                "kind": kind,
                "name": name,
                "count": data["count"],
                "axes": ",".join(data.get("axes", [])),
                "values": ",".join(str(value) for value in data.get("values", [])),
                "centroid_x": center["x"],
                "centroid_y": center["y"],
                "centroid_z": center["z"],
                **bounds,
            })

    check_rows = [
        {
            "kind": "check",
            "name": check["name"],
            "count": "",
            "axes": "",
            "values": "PASS" if check["passed"] else "FAIL",
            "centroid_x": "",
            "centroid_y": "",
            "centroid_z": "",
            "x_min": "",
            "x_max": "",
            "y_min": "",
            "y_max": "",
            "z_min": "",
            "z_max": "",
        }
        for check in checks
    ]
    rows.extend(check_rows)

    with open(path, "w", newline="") as handle:
        fieldnames = list(rows[0].keys()) if rows else []
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_extract_sets(model_file, command):
    subprocess.run([str(command), str(model_file)], check=True)


def read_all_node_coordinates(root):
    return root.groups["Parts"].groups["Part1"].variables["NodeCoordinates"][:]


def read_node_set_points(root):
    out = {}
    for name, group in root.groups["Sets"].groups["NodeSets"].groups.items():
        out[name] = _node_coordinates(root, group.variables["NodeNumber"][:])
    return out


def _axis_range(coords, axis_index):
    return float(coords[:, axis_index].min()), float(coords[:, axis_index].max())


def _plot_projection(ax, coords, node_set_points, projection, displacement):
    axis_a, axis_b, label_a, label_b = projection
    colors = {
        FEMORAL_HEAD_NODE_SET: "#d62728",
        GREATER_TROCHANTER_NODE_SET: "#1f77b4",
        DISTAL_FEMUR_NODE_SET: "#2ca02c",
        "body_top": "#d62728",
        "body_bottom": "#1f77b4",
    }
    labels = {
        FEMORAL_HEAD_NODE_SET: "Head",
        GREATER_TROCHANTER_NODE_SET: "GT",
        DISTAL_FEMUR_NODE_SET: "Shaft",
        "body_top": "Top",
        "body_bottom": "Bottom",
    }

    stride = max(1, int(np.ceil(coords.shape[0] / 30000.0)))
    sample = coords[::stride]
    ax.scatter(sample[:, axis_a], sample[:, axis_b], s=0.3, c="0.45", alpha=0.38, linewidths=0)

    for name, points in node_set_points.items():
        color = colors.get(name, "black")
        label = labels.get(name, name.replace("_PMMA_Nodes", "").replace("_Nodes", "").replace("_", " "))
        ax.scatter(points[:, axis_a], points[:, axis_b], s=4, c=color, alpha=0.9, linewidths=0, label=label)

    if displacement is not None:
        center = displacement["center"]
        vector = displacement["vector"]
        scale = 0.18 * max(
            _axis_range(coords, axis_a)[1] - _axis_range(coords, axis_a)[0],
            _axis_range(coords, axis_b)[1] - _axis_range(coords, axis_b)[0],
        )
        vec2 = np.array([vector[axis_a], vector[axis_b]], dtype=float)
        norm = float(np.linalg.norm(vec2))
        if norm > 0:
            vec2 = vec2 / norm * scale
            ax.arrow(
                center[axis_a],
                center[axis_b],
                vec2[0],
                vec2[1],
                color="#d62728",
                width=0.35,
                head_width=2.5,
                length_includes_head=True,
                alpha=0.95,
            )

    amin, amax = _axis_range(coords, axis_a)
    bmin, bmax = _axis_range(coords, axis_b)
    apad = max((amax - amin) * 0.05, 1.0)
    bpad = max((bmax - bmin) * 0.05, 1.0)
    ax.set_xlim(amin - apad, amax + apad)
    ax.set_ylim(bmin - bpad, bmax + bpad)
    ax.set_xlabel(f"{label_a} (mm)")
    ax.set_ylabel(f"{label_b} (mm)")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="0.9", linewidth=0.4)
    ax.tick_params(axis="both", which="major", labelsize=8, length=2.5, pad=1.5)


def _displacement_vector(node_sets, constraints):
    top = constraints.get("top_displacement")
    moving_set = node_sets.get(FEMORAL_HEAD_NODE_SET) or node_sets.get("body_top")
    if not top or not moving_set:
        return None
    axis = _single_axis(constraints, "top_displacement")
    value = _single_value(constraints, "top_displacement")
    axis_index = {"x": 0, "y": 1, "z": 2}.get(axis)
    if axis_index is None or value in (None, 0):
        return None
    vector = np.zeros(3, dtype=float)
    vector[axis_index] = float(np.sign(value))
    center = np.array(
        [moving_set["centroid"]["x"], moving_set["centroid"]["y"], moving_set["centroid"]["z"]],
        dtype=float,
    )
    return {"center": center, "vector": vector}


def write_png(path, coords, node_set_points, node_sets, constraints, checks):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    plt.rcParams.update({
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "font.size": 8,
        "legend.fontsize": 7,
    })
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.4), constrained_layout=True)
    displacement = _displacement_vector(node_sets, constraints)

    projections = [
        (0, 1, "X", "Y"),
        (0, 2, "X", "Z"),
        (1, 2, "Y", "Z"),
    ]
    titles = ["X-Y", "X-Z", "Y-Z"]
    for ax, projection, title in zip(axes.flat[:3], projections, titles):
        _plot_projection(ax, coords, node_set_points, projection, displacement)
        ax.set_title(title, pad=2)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4))

    legend_handles, legend_labels = axes.flat[0].get_legend_handles_labels()
    if legend_handles:
        axes.flat[0].legend(
            legend_handles,
            legend_labels,
            loc="upper left",
            fontsize=7,
            frameon=True,
            borderpad=0.2,
            labelspacing=0.25,
            handlelength=0.8,
            handletextpad=0.3,
            markerscale=1.2,
        )

    text_ax = axes.flat[3]
    text_ax.axis("off")
    lines = []
    for check in checks:
        status = "OK" if check["passed"] else "FAIL"
        lines.append(f"{status}: {check['name']}")
    text_ax.text(
        0.0,
        1.0,
        "\n".join(lines),
        va="top",
        ha="left",
        family="monospace",
        fontsize=7,
        linespacing=1.25,
        transform=text_ax.transAxes,
    )

    fig.savefig(path, dpi=300)
    plt.close(fig)


def audit_model(
    model_file,
    *,
    model="femur-sideways",
    output_prefix=None,
    flat_tolerance=1.0e-4,
    extract_sets=False,
    n88extractsets_command="n88extractsets",
    write_json=True,
    write_csv_file=True,
    write_plot=True,
):
    model_file = Path(model_file)
    output_prefix = Path(output_prefix) if output_prefix is not None else model_file.with_suffix("")
    json_path = output_prefix.with_name(output_prefix.name + "_bc_audit.json")
    csv_path = output_prefix.with_name(output_prefix.name + "_bc_audit.csv")
    png_path = output_prefix.with_name(output_prefix.name + "_bc_audit.png")

    if extract_sets:
        run_extract_sets(model_file, n88extractsets_command)

    with Dataset(str(model_file)) as root:
        coords = read_all_node_coordinates(root)
        node_set_points = read_node_set_points(root)
        node_sets = read_node_sets(root)
        constraints = read_constraints(root)

    checks = audit_checks(model, node_sets, constraints, flat_tolerance)
    passed = all(check["passed"] for check in checks)
    summary = {
        "model_file": str(model_file),
        "model": model,
        "passed": passed,
        "node_sets": node_sets,
        "constraints": constraints,
        "checks": checks,
    }

    if write_json:
        json_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    if write_csv_file:
        write_csv(csv_path, node_sets, constraints, checks)
    if write_plot:
        write_png(png_path, coords, node_set_points, node_sets, constraints, checks)
    return {
        "passed": passed,
        "json_path": json_path if write_json else None,
        "csv_path": csv_path if write_csv_file else None,
        "png_path": png_path if write_plot else None,
        "checks": checks,
        "summary": summary,
    }


def main():
    parser = argparse.ArgumentParser(
        prog="ogoCheckFEModelBC",
        description="Audit n88model node sets, constraints, and loading direction.",
    )
    parser.add_argument("model_file", type=Path)
    parser.add_argument(
        "--model",
        choices=["femur-sideways", "spine-compression"],
        default="femur-sideways",
    )
    parser.add_argument("--output_prefix", type=Path, default=None)
    parser.add_argument("--flat_tolerance", type=float, default=1.0e-4)
    parser.add_argument("--extract_sets", action="store_true")
    parser.add_argument("--n88extractsets_command", default="n88extractsets")
    parser.add_argument("--no_png", action="store_true", help="Skip writing the PNG quick-look audit plot.")
    args = parser.parse_args()

    result = audit_model(
        args.model_file,
        model=args.model,
        output_prefix=args.output_prefix,
        flat_tolerance=args.flat_tolerance,
        extract_sets=args.extract_sets,
        n88extractsets_command=args.n88extractsets_command,
        write_plot=not args.no_png,
    )

    for check in result["checks"]:
        status = "PASS" if check["passed"] else "FAIL"
        print(f"{status}: {check['name']} ({check['detail']})")
    print(f"Wrote {result['json_path']}")
    print(f"Wrote {result['csv_path']}")
    if result["png_path"] is not None:
        print(f"Wrote {result['png_path']}")
    raise SystemExit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
