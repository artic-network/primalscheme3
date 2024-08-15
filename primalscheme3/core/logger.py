# Sets up the logger for the application

import pathlib
import sys

from loguru import logger

log = logger.opt(colors=True)


def setup_logger(OUTPUT_DIR: pathlib.Path | None):
    """
    Sets up the logger for the application
    :param OUTPUT_DIR: The output directory or None for stdout only
    :return: The loguru logger
    """
    ## Set up the loggers
    log.remove()  # Remove default stderr logger

    if OUTPUT_DIR is not None:
        # Add the deep log file
        log.add(
            OUTPUT_DIR / "work/file.log",
            colorize=False,
            format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}",
            enqueue=True,
        )

    # Add the nice stdout progress
    log.add(
        sys.stdout,
        colorize=True,
        format="{time:HH:mm:ss}|{level}| {message}",
        level="INFO",
    )

    return log


def setup_rich_logger(logfile: str | None = None):
    import logging

    from rich.console import Console
    from rich.logging import RichHandler

    handlers = [RichHandler(level=logging.INFO, markup=True, show_path=False)]
    # Add an optional file handler
    if logfile is not None:
        handlers.append(
            RichHandler(
                level=logging.DEBUG,
                console=Console(file=open(logfile, "w")),
                markup=True,
                show_path=False,
                omit_repeated_times=False,
            )
        )
    logging.basicConfig(
        format="%(message)s",
        datefmt="[%X]",
        handlers=handlers,
        level=logging.DEBUG,
    )
    log = logging.getLogger("rich")
    return log
