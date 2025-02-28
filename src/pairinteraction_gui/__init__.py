__version__ = "0.1.0"


def main() -> None:
    """Entry point for the application."""
    import sys

    from pairinteraction_gui.app import run_app

    sys.exit(run_app())
