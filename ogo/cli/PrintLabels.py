# /------------------------------------------------------------------------------+
# | 2-AUG-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

"""Format Ogo label definitions and expose them through ``ogo print-labels``."""

import argparse

import ogo.dat.OgoMasterLabels as lb


def render_labels(echo="labels"):
    """Return a complete table with a final newline, without printing.

    ``echo`` selects ``labels``, ``colors16``, or ``colors128``.
    Unknown selections raise ValueError. Definitions and row order come from
    OgoMasterLabels; formatting matches the original exported tables.
    """
    if echo == "labels":
        header = lb.labels_hdr
        rows = "\n".join(
            '{:5d} {:4d} {:4d} {:4d}  {:.1f} {:4d} {:4d}  "{}"'.format(
                label, value["RGB"][0], value["RGB"][1], value["RGB"][2],
                value["A"], value["VIS"], value["IDX"], value["LABEL"],
            )
            for label, value in lb.labels_dict.items()
        )
    else:
        if echo == "colors16":
            header, definitions = lb.colors16_hdr, lb.colors16_dict
        elif echo == "colors128":
            header, definitions = lb.colors128_hdr, lb.colors128_dict
        else:
            raise ValueError(f"Unknown label table: {echo!r}")
        rows = "\n".join(
            "{:5d} {:4d}  {:4d}  {:4d}   {:s}".format(
                label, value["RGB"][0], value["RGB"][1], value["RGB"][2],
                value["DESC"],
            )
            for label, value in definitions.items()
        )
    return header + "\n" + rows + "\n"


def add_arguments(parser):
    """Define the table selection for this command."""
    parser.add_argument(
        "--echo", default="labels", choices=["labels", "colors16", "colors128"],
        help="Specify output (default: %(default)s)",
    )


def main(argv=None):
    """Parse command arguments, write the selected table, and return success."""
    parser = argparse.ArgumentParser(
        prog="ogo print-labels",
        description=(
            "Print standard Ogo labels for bones, tissues, and calibration rods, "
            "or a 16/128-color table. Definitions come from OgoMasterLabels.py. "
            "The label output can be imported in ITK-SNAP through "
            "Segmentation > Import Label Descriptions."
        ),
        epilog=(
            "Examples:\n"
            "  ogo print-labels --echo labels > ogo-master-labels.txt\n"
            "  ogo print-labels --echo colors16"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    add_arguments(parser)
    args = parser.parse_args(argv)
    print(render_labels(args.echo), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
