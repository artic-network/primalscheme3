# Sets up the logger for the application

import pathlib
import sys

from loguru import logger

log = logger.opt(colors=True)


def setup_loger(OUTPUT_DIR: pathlib.Path):
    """
    Sets up the logger for the application
    :param OUTPUT_DIR: The output directory
    :return: The loguru logger
    """
    ## Set up the loggers
    log.remove()  # Remove default stderr logger
    # Add the deep log file
    log.add(
        OUTPUT_DIR / "work/file.log",
        colorize=False,
        format="{time:YYYY-MM-DD at HH:mm:ss} | {level} | {message}",
        enqueue=True,
    )
    # Add the nice stdout progress
    log.add(sys.stdout, colorize=True, format="{message}", level="INFO")

    return log
