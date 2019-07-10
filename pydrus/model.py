"""
This file contains the model class.
"""

import pandas as pd
import os
import subprocess


class Model:
    def __init__(self, exe_name, ws_name, name="model", description=None,
                 length_unit="m", time_unit="days", mass_units="mmol"):
        """Basic Hydrus model container.

        Parameters
        ----------
        name: str, optional
            String with the name of the model.
        description: str, optional
            String with the description of the model.
        length_unit: str, optional
            length units to use in the simulation. Options are "mm", "cm",
            and "m". Defaults to "cm".
        time_unit: str, optional
            time unit to use in the simulation, options are "seconds",
            "minutes", "hours", "days, "years". Defaults to "days".
        mass_units: str, optional
            Mass units to use in the simulation, Options are "mmol".
            Defaults to "mmol". Only used when transport process is added.

        """
        self.exe_name = exe_name
        self.ws_name = ws_name

        self.name = name
        self.description = description
        self.settings = dict()
        self.units = {
            "LUnit": length_unit,
            "TUnit": time_unit,
            "MUnit": mass_units
        }

        # The main processes to describe the simulation
        self.water_flow = None
        self.soil_profile = None
        self.time_information = None

        # The following processes will be implemented at a later stage.
        self.solute_transport = None
        self.heat_transport = None
        self.rootwater_uptake = None
        self.root_growth = None

    def add_profile(self):
        """Method to add the soil profile to the model.

        """
        pass

    def simulate(self):
        """Method to call the Hydrus-1D executable.

        """
        # 1. Check model

        # 2. Write files
        self.write_files()

        # 3. Run Hydrus executable.
        cmd = [self.exe_name, self.ws_name, "-1"]
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)

        # 4. Read output and check for success.
        return result

    def write_selector(self):
        """Write the selector.in file.

        """
        # Write block A: BASIC INFORMATION

        # Write block B: WATER FLOW INFORMATION

        # Write BLOCK C: TIME INFORMATION

        pass

    def write_hydrus(self):
        """Write the hydrus1s.dat file needed for simulation.

        """
        pass

    def write_files(self):
        pass

    def read_output(self):
        pass

    def read_profile(self, fname="PROFILE.OUT"):
        """

        Parameters
        ----------
        fname: str, optional
            String with the name of the profile out file. default is
            "PROFILE.OUT".

        Returns
        -------
        data: pandas.DataFrame
            Pandas with the profile data

        """
        path = os.path.join(self.ws_name, fname)
        os.path.exists(path)

        with open(path) as file:
            # Find the starting line to read the profile
            for start, line in enumerate(file.readlines(1000)):
                if "depth" in line:
                    break
            file.seek(0) # Go back to start of file
            # Read the profile data into a Pandas DataFrame
            data = pd.read_csv(file, skiprows=start, skipfooter=2, index_col=0,
                               skipinitialspace=True, delim_whitespace=True)
        return data

    def read_tlevel(self, fname="T_LEVEL.OUT"):
        """

        Parameters
        ----------
        fname: str, optional
            String with the name of the t_level out file. default is
            "T_LEVEL.OUT".

        Returns
        -------
        data: pandas.DataFrame
            Pandas with the t_level data

        """
        path = os.path.join(self.ws_name, fname)
        os.path.exists(path)

        with open(path) as file:
            # Find the starting line to read the profile
            for start, line in enumerate(file.readlines(1000)):
                if "rTop" in line:
                    break
            file.seek(0) # Go back to start of file
            # Read the profile data into a Pandas DataFrame
            data = pd.read_csv(file, skiprows=start, skipfooter=2, index_col=0,
                               skipinitialspace=True, delim_whitespace=True)
            # Fix the header
            data.columns = [col+str(col1) for col, col1 in data.iloc[0].items()]
            data = data.drop(data.index[0]).apply(pd.to_numeric)
        return data
