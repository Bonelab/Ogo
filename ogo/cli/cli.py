"""Select a migrated command without importing unrelated scientific workflows."""

import argparse
from contextlib import contextmanager
import importlib
import logging
import sys
import time

from ogo import __version__


COMMANDS = {
    "print-labels": ("ogo.cli.PrintLabels", "Print standard label definitions and color tables"),
    "check-labels": ("ogo.cli.CheckLabels", "Inspect segmentation labels or check required IDs"),
    "resample-image": ("ogo.cli.IsotropicResampling", "Resample a CT image or mask"),
}


class _ElapsedFormatter(logging.Formatter):
    """Prefix diagnostics with monotonic elapsed seconds and align extra lines."""

    def __init__(self, started):
        super().__init__("%(message)s")
        self.started = started

    def format(self, record):
        message = super().format(record)
        prefix = f"{time.monotonic() - self.started:8.2f}s [{record.levelname}] "
        return prefix + message.replace("\n", "\n" + " " * len(prefix))


@contextmanager
def _command_logging(level, started):
    """Own one stderr handler for this invocation, restoring application state."""
    logger = logging.getLogger("ogo")
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(level)
    handler.setFormatter(_ElapsedFormatter(started))
    previous_level, previous_propagate = logger.level, logger.propagate
    logger.addHandler(handler)
    logger.setLevel(level)
    logger.propagate = False
    try:
        yield
    finally:
        logger.removeHandler(handler)
        handler.close()
        logger.setLevel(previous_level)
        logger.propagate = previous_propagate


def main(argv=None):
    """Dispatch the selected command and return its exit status."""
    started = time.monotonic()
    command_help = "Commands:\n" + "\n".join(
        f"  {name:<16} {description}"
        for name, (_, description) in COMMANDS.items()
    )
    parser = argparse.ArgumentParser(
        prog="ogo",
        description="Ogo modernization pilot: run migrated research commands.",
        epilog=command_help,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument(
        "--verbose", dest="log_level", action="store_const", const=logging.DEBUG,
        help="Show debug diagnostics (before the command)",
    )
    verbosity.add_argument(
        "--quiet", dest="log_level", action="store_const", const=logging.WARNING,
        help="Show only warnings and errors (before the command)",
    )
    parser.set_defaults(log_level=logging.INFO)
    parser.add_argument("command", nargs="?", choices=COMMANDS, help="Command to run")
    parser.add_argument("arguments", nargs=argparse.REMAINDER, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 0

    module_name, _ = COMMANDS[args.command]
    with _command_logging(args.log_level, started):
        module = importlib.import_module(module_name)
        return module.main(args.arguments)


if __name__ == "__main__":
    raise SystemExit(main())
