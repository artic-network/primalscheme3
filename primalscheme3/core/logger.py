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
        level="NOTSET", format="%(message)s", datefmt="[%X]", handlers=handlers
    )
    log = logging.getLogger(f"rich-{hash(logfile)}")

    return log
