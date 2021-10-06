"""The utils module contains utility funtions for Phydrus.

"""

import logging
from logging import handlers
from numpy import exp, maximum, argmin, array, flip


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
            
def partitioning_grass(P, ET, a=0.45, ch=5, k=0.463, return_SCF=False):
    
    """
    Partitioning according to equation 2.75 in the Manual v4.0 and 
    Sutanto, Wenninger, Coenders and Uhlenbrook [2021]
    
    Parameters
    ----------
    P - precipitation [cm]
    ET - potential evapotranspiration (Penman-Monteith) [cm]
    a - constant [cm]
    ch - cropheight (5-15 for clipped grass) [cm]
    k - radiation extinction by canopy (rExtinct) (0.463) [-]
    
    Internal variables
    ----------
    LAI - Leaf Area Index (0.24*ch) [cm/cm]
    SCF - Soil Cover Fraction (b) [-]
   
    Standard Output
    ----------
    Pnet - Net Precipitation (P - I) [cm]
    I - Interception [cm]
    Et,p - Potential Transpiration (rRoot) [cm]
    Es,p - Potential Soil Evaporation (rSoil) [cm]
    """
    
    LAI = 0.24 * ch
    SCF = 1 - exp(-k * LAI)
    I = a * LAI * (1 - 1 / (1 + SCF * P / (a * LAI)))
    Pnet = maximum(P - I, 0)
    Ep = maximum(ET - I, 0)
    Etp = Ep * SCF
    Esp = Ep * (1 - SCF)
    if return_SCF==True:
        return Pnet, I, Etp, Esp, SCF
    else:
        return Pnet, I, Etp, Esp
    
def get_recharge(nodinf, recharge_depth):
    """
    Function to obtain the flux at a certain depth from the NOD_INF.OUT file.
    
    INPUT:
    nodinf - Dictionary of the NOD_INF.OUT file with the timesteps as keys. 
             The data in the dictionary contains a DataFrame for each timestep 
             with the column 'Depth' and 'Flux'. 
             Take a look at read.read_nodinf()
    recharge_depth - Depth at which the recharge is extracted.

    OUTPUT:
    recharge - Numpy array with the recharge at a defined depth.
    """
    # recharge depth has to be negative
    idx = argmin(abs(nodinf[0]['Depth'].values-recharge_depth))
    recharge = []
    for i in list(nodinf.keys()):
        recharge.append(nodinf[i]['Flux'].iloc[idx])
    return array(recharge)

def get_gwt(nodinf):
    """
    Extract the location of the groundwater table from the NOD_INF.OUT file.
    The location of the groundwater table is taken where the pressure head = 0 (bottom up approach)
    
    INPUT:
    nodinf - Dictionary of the NOD_INF.OUT file with the timesteps as keys. 
             The data in the dictionary contains a DataFrame for each timestep 
             with the column 'Depth' and 'Head'.
    
    OUTPUT:
    gwl - Numpy array with the location of the groundwater table.
    """
    gwl = []
    for i in list(nodinf.keys()):
        j = 0
        while flip(nodinf[i]['Head'].values)[j] > 0 and j < len(nodinf[i])-1:
            j += 1
        idx = len(nodinf[0]['Head']) - j
        if idx == 0:
            gwl_new = nodinf[i]['Head'][idx]
        elif idx == len(nodinf[i]['Head']):
            gwl_new = nodinf[i]['Depth'].iloc[-1]
        else:
            gwl_new = (0-nodinf[i]['Head'].iloc[idx-1])*(nodinf[i]['Depth'].iloc[idx]-nodinf[i]['Depth'].iloc[idx-1])/(nodinf[i]['Head'].iloc[idx]-nodinf[i]['Head'].iloc[idx-1]) + nodinf[i]['Depth'].iloc[idx-1]
        gwl.append(gwl_new)
    return array(gwl)