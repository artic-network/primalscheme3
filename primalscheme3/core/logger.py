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
    )
    log = logging.getLogger(f"rich-{hash(logfile)}")
    log.setLevel(logging.DEBUG)

    # Add the handlers
    for handler in handlers:
        log.addHandler(handler)

    return log
