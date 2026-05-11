import sys
from loguru import logger


def setup_logging(verbose: bool) -> None:
    logger.remove()
    if verbose:
        logger.add(
            sys.stderr,
            level="DEBUG",
            format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan> - <level>{message}</level>",
        )
    else:
        logger.add(sys.stderr, level="WARNING", format="<level>{message}</level>")
