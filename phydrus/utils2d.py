import numpy as np
from pandas import read_csv


def read_mesh_bin(path="MESHTRIA.000"):
    """Generate x, y and triangle from the MESHTRIA.000 file from HYDRUS-2D.

    Returns
    -------
    x, y, triangles (numpy array):
        x and y locations of the nodes and the trianglular mesh of the cells.

    """
    with open(path, "rb") as f:
        _ = f.read(17)
        _ = np.fromfile(f, np.int16, 1)

        nP = np.fromfile(f, np.int64, 1)  # nNodes
        nE = np.fromfile(f, np.int64, 1)  # nCells
        x = np.empty(nP)
        y = np.empty(nP)

        for i in np.arange(nP):
            x[i] = np.fromfile(f, np.float32, 1)
            y[i] = np.fromfile(f, np.float32, 1)
            _ = np.fromfile(f, np.int64, 1)  # rdummy, idummy?

        tri = np.reshape(np.fromfile(f, np.int32), [nE[0], 3]) - 1
        return x, y, tri


def read_mesh_txt(path="MESHTRIA.TXT"):
    """Generate x, y and triangle from the MESHTRIA.TXT file from HYDRUS-2D.

    Returns
    -------
    x, y, triangles (numpy array):
        x and y locations of the nodes and the trianglular mesh of the cells.

    """
    with open(path) as fo:
        line = fo.readline().split()
        nNodes = int(line[1])
        nCells = int(line[3])
        x = np.empty(nNodes)
        y = np.empty(nNodes)

        for i in range(nNodes):
            line = fo.readline().split()
            x[i] = float(line[1])
            y[i] = float(line[2])

        for _ in range(3):
            line = fo.readline()

        tri = np.empty([nCells, 3])
        for i in range(nCells):
            line = fo.readline().split()
            tri[i, 0] = int(line[1])
            tri[i, 1] = int(line[2])
            tri[i, 2] = int(line[3])
        tri = tri - 1  # Python zero based
    return x, y, tri


def read_mesh(path):
    """Generate x, y and triangle lists from the MESHTRIA.000 or MESHTIRA.TXT
    file from HYDRUS-2D

    Returns
    -------
    x, y, triangles (numpy array):
        x and y locations of the nodes and the trianglular mesh of the cells.

    """
    extension = path.split(".")[-1].lower()

    if extension == "000":
        x, y, tri = read_mesh_bin(path)

    elif extension == "txt":
        x, y, tri = read_mesh_txt(path)

    return x, y, tri


def read_obsnod_out(path="ObsNod.out", col="hNew"):
    """Read the ObsNod.out file of HYDRUS-2D

    Parameters
    ----------
    path : str, optional
        Path to the file, by default 'ObsNod.out'
    col : str, optional
        Column to obtain from ObsNod.out file, by default 'hNew'
        Choice between: 'hNew', 'theta' & 'Temp'

    Returns
    -------
    DataFrame with nodes as columns and timesteps as index (pandas DataFrame)
    """

    with open(path) as f:
        for _ in range(4):
            line = f.readline()
        nodes = [int(node[5:-1]) for node in line.split()]
        df = read_csv(
            f, delim_whitespace=True, index_col=0, skipfooter=1, engine="python"
        )
    dtf = df.loc[:, [string for string in df.columns if col in string]].copy()
    dtf.columns = nodes
    return dtf


def read_bin_out(path="h.out"):
    """Read binary h.out, th.out or v.out file from HYDRUS-2D

    Parameters
    ----------
    path : str, optional
        Path to the file, by default 'h.out'

    Returns
    ----------
        array with floats32
    """
    return np.fromfile(path, "float32")


def read_h_out(path="h.out"):
    """Read binary h.out file from HYDRUS-2D

    Parameters
    ----------
    path : str, optional
        Path to the file, by default 'h.out'

    Returns
    ----------
        array with floats32
    """
    return read_bin_out(path)


def read_th_out(path="th.out"):
    """Read binary th.out file from HYDRUS-2D

    Parameters
    ----------
    path : str, optional
        Path to the file, by default 'th.out'

    Returns
    ----------
        array with floats32
    """
    return read_bin_out(path)


def read_v_out(path="v.out"):
    """Read binary v.out file from HYDRUS-2D

    Parameters
    ----------
    path : str, optional
        Path to the file, by default 'v.out'

    Returns
    ----------
        array with floats32
    """
    return read_bin_out(path)
