"""The utils module contains utility funtions for Phydrus.

"""

import logging
from logging import handlers

logger = logging.getLogger(__name__)


def show_versions():
    """
    Method to print the version of dependencies.

    Examples
    --------
    >>> import phydrus as ps
    >>> ps.show_versions()

    Python version: 3.8.2 (default, Mar 25 2020, 11:22:43)
    [Clang 4.0.1 (tags/RELEASE_401/final)]
    Numpy version: 1.19.2
    Pandas version: 1.2.1
    Phydrus version: 0.1.0
    Matplotlib version: 3.3.2

    """
    from phydrus import __version__ as ps_version
    from pandas import __version__ as pd_version
    from numpy import __version__ as np_version
    from matplotlib import __version__ as mpl_version
    from sys import version as os_version

    msg = (
        f"Python version: {os_version}\n"
        f"Numpy version: {np_version}\n"
        f"Pandas version: {pd_version}\n"
        f"Phydrus version: {ps_version}\n"
        f"Matplotlib version: {mpl_version}"
    )

    return print(msg)


def _initialize_logger(logger=None, level=logging.INFO):
    """
    Internal method to create a logger instance to log program output.

    Parameters
    -------
    logger : logging.Logger
        A Logger-instance. Use ps.logger to initialise the Logging instance
        that handles all logging throughout Phydrus,  including all sub modules
        and packages.

    """
    if logger is None:
        logger = logging.getLogger('phydrus')

    logger.setLevel(level)
    remove_file_handlers(logger)
    set_console_handler(logger)


def set_console_handler(logger=None, level=logging.INFO,
                        fmt="%(levelname)s: %(message)s"):
    """
    Method to add a console handler to the logger of Phydrus.

    Parameters
    -------
    logger : logging.Logger
        A Logger-instance. Use ps.logger to initialise the Logging instance
        that handles all logging throughout Phydrus, including all sub modules
        and packages.

    """
    if logger is None:
        logger = logging.getLogger('phydrus')
    remove_console_handler(logger)
    ch = logging.StreamHandler()
    ch.setLevel(level)
    formatter = logging.Formatter(fmt=fmt)
    ch.setFormatter(formatter)
    logger.addHandler(ch)


def set_log_level(level):
    """
    Set the log-level for which to log Phydrus messages.

    Parameters
    ----------
    level: str
        String with the level to log messages to the screen for. Options
        are: "INFO", "WARNING", and "ERROR".

    Examples
    --------
    >>> import phydrus as ps
    >>> ps.set_log_level("ERROR")

    """
    set_console_handler(level=level)


def remove_console_handler(logger=None):
    """
    Method to remove the console handler to the logger of Phydrus.

    Parameters
    ----------
    logger : logging.Logger
        A Logger-instance. Use ps.logger to initialise the Logging instance
        that handles all logging throughout Phydrus, including all sub modules
        and packages.

    """
    if logger is None:
        logger = logging.getLogger('phydrus')

    for handler in logger.handlers:
        if isinstance(handler, logging.StreamHandler):
            logger.removeHandler(handler)


def add_file_handlers(logger=None, filenames=('info.log', 'errors.log'),
                      levels=(logging.INFO, logging.ERROR), maxBytes=10485760,
                      backupCount=20, encoding='utf8', datefmt='%d %H:%M',
                      fmt='%(asctime)s-%(name)s-%(levelname)s-%(message)s'):
    """
    Method to add file handlers in the logger of Phydrus.

    Parameters
    -------
    logger : logging.Logger
        A Logger-instance. Use ps.logger to initialise the Logging instance
        that handles all logging throughout Phydrus, including all sub modules
        and packages.

    """
    if logger is None:
        logger = logging.getLogger('phydrus')
    # create formatter
    formatter = logging.Formatter(fmt=fmt, datefmt=datefmt)

    # create file handlers, set the level & formatter, and add it to the logger
    for filename, level in zip(filenames, levels):
        fh = handlers.RotatingFileHandler(filename, maxBytes=maxBytes,
                                          backupCount=backupCount,
                                          encoding=encoding)
        fh.setLevel(level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)


def remove_file_handlers(logger=None):
    """
    Method to remove any file handlers in the logger of Phydrus.

    Parameters
    -------
    logger : logging.Logger
        A Logger-instance. Use ps.logger to initialise the Logging instance
        that handles all logging throughout Phydrus, including all sub modules
        and packages.
    """
    if logger is None:
        logger = logging.getLogger('phydrus')
    for handler in logger.handlers:
        if isinstance(handler, handlers.RotatingFileHandler):
            logger.removeHandler(handler)
