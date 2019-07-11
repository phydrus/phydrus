
import os
import pandas as pd

class Profile:
    def __init__(self, data=None, top=0, bot=-1, dx=0.1, lay=1, mat=1,
                 beta=1, ah=1, ak=1, ath=1, temp=None, conc=None, sconc=None):
        """

        Parameters
        ----------
        data
        top: float, optional
            Top of the soil column
        bot: float or list of float, optional
            Bottom of the soil column. If a list is provided, multiple
            layers are created and other arguments need to be of the same
            length (e.g. mat).
        dx: float: optional
            Size of each grid cell. Default 0.1 meter.
        lay: int or list of int, optional
            subregion number (for mass balance calculations)
        mat: int or list of int, optional
            Material number (for heterogeneity)
        beta: float or list of float, optional

        ah: float or list of float, optional
            Scaling factor for the pressure head (Axz in profile.dat).
        ak: float or list of float, optional
            Scaling factor the hydraulic conductivity (Bxz in profile.dat).
        ath: float or list of float, optional
            Scaling factor the the water content (Dxz in profile.dat).
        temp: float, optional
        conc: float, optional
        sconc: float, optional

        """
        self.data = data
        self.observations = []

        # Store all parameter to create the profile
        self.top = top
        self.bot = bot
        self.dx = dx
        self.lay = lay
        self.mat = mat
        self.beta = beta
        self.ah = ah
        self.ak = ak
        self.ath = ath
        self.temp = temp
        self.conc = conc
        self.sconc = sconc

    def create_data(self):

        pass
        cols = ["x", "h", "Mat", "Lay", "Beta", "Axz", "Bxz", "Dxz", "Temp",
                "Conc", "SConc"]
        data = pd.DataFrame(columns=cols)

        return data

    def add_observation(self, point):
                pass

    def write_file(self, fname="PROFILE_TEST.DAT", ws=""):
        """Method to write the profile.dat file

        """
        path = os.path.join(ws, fname)

        if os.path.exists(path):
            os.remove(path)

        with open(path, mode="w") as file:
            # Write start


            # Write data


            # Write observation points
            file.write(str(len(self.observations)) + os.linesep)
            file.write("".join(["   {}".format(i) for i in self.observations]))

    @classmethod
    def read_profile(self, fname="PROFILE.DAT", ws=None):
        """Method to read a profile.dat file

        Parameters
        ----------
        fname
        ws

        Returns
        -------
        profile: pydrus.profile.SoilProfile

        """
        path = os.path.join(ws, fname)
        os.path.exists(path)

        with open(path) as file:
            # Find the starting line to read the profile
            for start, line in enumerate(file.readlines(1000)):
                if "Beta" in line:
                    break
            file.seek(0)  # Go back to start of file
            # Read the profile data into a Pandas DataFrame
            data = pd.read_csv(file, skiprows=start, skipfooter=2, index_col=0,
                               skipinitialspace=True, delim_whitespace=True)

        return data

