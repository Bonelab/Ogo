"""Select a migrated command without importing unrelated scientific workflows."""

import argparse
import importlib

from ogo import __version__


"""
Command-line interface for the Ogo modernization pilot. This defines a set of commands available in the new "ogo" entry point.

Required demonstration of the new CLI:
- Read: Demonstrate that the new workflow reads an image (i.e. NIfTI, aims, etc.) and prints its metadata
- Write: Demonstrate that the new workflow writes an image
- Apply filter: Demonstrate that the new workflow applies a filter to an image (e.g., Gaussian smoothing)
- Threshold: Demonstrate that the new workflow thresholds an image (e.g., binary thresholding)
---
Potential candidate commands for demonstration:
- ogo.cli.ReplaceLabels: [Reads] an image, modifies its labels, and [writes] the output image.
---
Migrated commands:
- ogo.cli.PrintLabels: Print standard label definitions and color tables
---
"""

COMMANDS = {
    "print-labels": ("ogo.cli.PrintLabels", "Print standard label definitions and color tables"),
}


def main(argv=None):
    """Dispatch the selected command and return its exit status."""
    command_help = "Commands:\n" + "\n".join(
        f"  {name:<16} {description}"
        for name, (_, description) in COMMANDS.items()
    )
    parser = argparse.ArgumentParser(
        prog="ogo",
        description="Ogo modernization pilot. <Placeholder description of the new ogo CLI>",
        epilog=command_help,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("command", nargs="?", choices=COMMANDS, help="Command to run")
    parser.add_argument("arguments", nargs=argparse.REMAINDER, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    # If no command is provided, print help and exit
    if args.command is None:
        parser.print_help()
        return 0

    # Dispatch to the selected command module
    module_name, _ = COMMANDS[args.command]
    # Import the module and call its main() function with the remaining arguments
    module = importlib.import_module(module_name)
    # Call the main function of the selected command module with the remaining arguments
    return module.main(args.arguments)


if __name__ == "__main__":
    raise SystemExit(main())
