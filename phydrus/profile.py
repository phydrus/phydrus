"""The profile module contains utility functions to help create the profile
DataFrame.

"""

from os import path

from numpy import linspace
from pandas import read_csv, DataFrame


def create_profile(top=0, bot=-1, dx=0.1, h=0, lay=1, mat=1, beta=0, ah=1.0,
                   ak=1.0, ath=1.0, temp=20.0, conc=None, sconc=None):
    """
    Method to create a DataFrame describing the soil profile.

    Parameters
    ----------
    top: float, optional
        Top of the soil column.
    bot: float or list of float, optional
        Bottom of the soil column. If a list is provided, multiple layers are
        created and other arguments need to be of the same length (e.g. mat).
    dx: float, optional
        Size of each grid cell. Default 0.1 meter.
    h: float, optional
        Initial values of the pressure head.
    lay: int or list of int, optional
        Subregion number (for mass balance calculations).
    mat: int or list of int, optional
        Material number (for heterogeneity).
    beta: float or list of float, optional
        Value of the water uptake distribution, b(x) [L-1], in the soil root
        zone at node n. Set Beta(n) equal to zero if node n lies outside the
        root zone.
    ah: float or list of float, optional
        Scaling factor for the pressure head (Axz in profile.dat).
    ak: float or list of float, optional
        Scaling factor the hydraulic conductivity (Bxz in profile.dat).
    ath: float or list of float, optional
        Scaling factor the the water content (Dxz in profile.dat).
    temp: float, optional
        Initial value of the temperature at node n [oC] (do not specify if
        both lTemp or lChem are equal to .false.; if lTemp=.false. and
        lChem=.true. then set equal to 0 or any other initial value to be used
        later for temperature dependent water flow and solute transport).
    conc: float, optional
    sconc: float, optional

    """
    if not isinstance(bot, list):
        bot = [bot]
    steps = int(top - bot[-1] / dx) + 1
    grid = linspace(top, bot[-1], steps)
    cols = ["x", "h", "Mat", "Lay", "Beta", "Axz", "Bxz", "Dxz", "Temp",
            "Conc", "SConc"]
    data = DataFrame(columns=cols)
    data["x"] = grid
    variables = [h, mat, lay, beta, ah, ak, ath, temp, conc, sconc]

    if len(bot) == 1:
        data[cols[1:]] = variables
    else:
        # If there are multiple layers
        for i, arg in enumerate(variables):
            if isinstance(arg, int) or isinstance(arg, float) or arg is None:
                variables[i] = [arg] * len(bot)

        for i, b in enumerate(bot):
            layer = ((data.loc[:, "x"] <= top) & (data.loc[:, "x"] > (b - dx)))
            data.loc[layer, cols[1:]] = [var[i] for var in variables]
            top = b
    data = data.fillna("")
    data.index = data.index + 1
    return data


def profile_from_file(fname="PROFILE.DAT", ws=None):
    """
    Method to create a profile DataFrame from a profile.dat file

    Parameters
    ----------
    fname: str, optional
        String with the path to the PROFILE.DAT file.
    ws: str, optional
        String with the path to the workspace.

    Returns
    -------
    data: pandas.DataFrame
        Pandas DataFrame with the soil profile data.

    Examples
    --------

    >>> profile = ps.create_profile(h=0.342)

    """
    fname = path.join(ws, fname)
    path.exists(fname)

    with open(fname) as file:
        # Find the starting line to read the profile
        for start, line in enumerate(file.readlines(1000)):
            if "Beta" in line:
                break
        file.seek(0)  # Go back to start of file
        # Read the profile data into a Pandas DataFrame
        data = read_csv(file, skiprows=start, skipfooter=2, index_col=0,
                        skipinitialspace=True, usecols=range(0, 11),
                        delim_whitespace=True, engine='python', names=None)
        data.columns = ["x", "h", "Mat", "Lay", "Beta", "Axz", "Bxz", "Dxz",
                        "Temp", "Conc"]
        data.index.name = None

    return data
