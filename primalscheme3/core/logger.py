import logging

from rich.console import Console, ConsoleRenderable
from rich.logging import RichHandler


class Handler(RichHandler):
    """Also logs Rich renderables (e.g. a Table) as themselves."""

    def render_message(self, record, message):
        if isinstance(record.msg, ConsoleRenderable):
            return record.msg
        return super().render_message(record, message)


def setup_rich_logger(logfile: str | None = None):
    """
    A logger writing INFO to the console and, if ``logfile`` is given, DEBUG
    to that file. Handlers are attached to this logger rather than through
    logging.basicConfig, which only takes effect once per process, so every
    run (e.g. several in one Python session) writes its own log file.
    """
    log = logging.getLogger(f"rich-{hash(logfile)}")
    log.setLevel(logging.DEBUG)
    # Replace the handlers of an earlier run with the same logfile
    for handler in list(log.handlers):
        log.removeHandler(handler)
        handler.close()

    handlers = [Handler(level=logging.INFO, markup=True, show_path=False)]
    # Add an optional file handler
    if logfile is not None:
        handlers.append(
            Handler(
                level=logging.DEBUG,
                console=Console(file=open(logfile, "w"), force_terminal=False),
                markup=True,
                show_path=False,
                omit_repeated_times=False,
            )
        )
    for handler in handlers:
        handler.setFormatter(logging.Formatter("%(message)s", datefmt="[%X]"))
        log.addHandler(handler)

    return log
