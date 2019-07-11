"""
This file contains the model class.
"""

import pandas as pd
import os
import subprocess
from numpy import logspace, log10


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
        # Store the hydrus executable and the project workspace
        self.exe_name = exe_name
        self.ws_name = ws_name

        self.name = name
        self.description = description

        self.basic_information = {
            "iVer": "0.1",
            "Hed": "Heading",
            "LUnit": length_unit,
            "TUnit": time_unit,
            "MUnit": mass_units,
            "lWat": True,
            "lChem": False,
            "lTemp": False,
            "lSink": False,
            "lRoot": False,
            "lShort": False,
            "lWDep": False,
            "lScreen": True,
            "AtmInf": False,
            "lEquil": True,
            "lInverse": False,
            "lSnow": False,
            "lHP1": False,
            "lMeteo": False,
            "lVapor": False,
            "lActRSU": False,
            "lFlux": False,
            "NMat": 1,
            "NLay": 1,
            "CosAlfa": 1,
        }

        self.time_information = {
            "dt": 0.001,
            "dtMin": 0.0001,
            "dtMax": 0.5,
            "dMul": 1.3,
            "dMul2": 0.7,
            "ItMin": 3,
            "ItMax": 7,
            "MPL": 4,
            "tInit": 0,
            "tMax": 1,
            "lPrint": False,
            "nPrintSteps": 1,
            "tPrintInterval": 1,
            "lEnter": False,
            "TPrint(1)": 0,
            "TPrint(2)": 1,
            "TPrint(MPL)": 1,
        }

        # The main processes to describe the simulation
        self.water_flow = None
        self.soil_profile = None

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
        result = subprocess.run(cmd)

        # 4. Read output and check for success.
        self.read_output()

        return result

    def write_files(self):
        self.write_selector()
        # self.profile.write_file()

    def write_selector(self, fname="SELECTOR.IN"):
        """Write the selector.in file.

        """
        # Create Header string
        string = "***{:{fill}{align}{width}}\n"

        lines = ["Pcp_File_Version=4\n"]

        # Write block A: BASIC INFORMATION
        lines.append(string.format(" BLOCK A: BASIC INFORMATION ", fill="*",
                                   align="<", width=72))
        lines.append("{}\n".format(self.basic_information["Hed"]))
        lines.append("{}\n".format(self.description))
        lines.append("LUnit  TUnit  MUnit  (indicated units are obligatory "
                     "for all input data)\n")
        lines.append("{}\n".format(self.basic_information["LUnit"]))
        lines.append("{}\n".format(self.basic_information["TUnit"]))
        lines.append("{}\n".format(self.basic_information["MUnit"]))

        vars2 = [['lWat', 'lChem', 'lTemp', 'lSink', 'lRoot', 'lShort',
                  'lWDep', 'lScreen', 'AtmInf', 'lEquil', 'lInverse', "\n"],
                 ['lSnow', 'lHP1', 'lMeteo', 'lVapor', 'lActRSU', 'lFlux',
                  "\n"]]

        for vars in vars2:
            lines.append("  ".join(vars))
            vals = []
            for var in vars[:-1]:
                val = self.basic_information[var]
                if val is True:
                    vals.append("t")
                elif val is False:
                    vals.append("f")
                else:
                    vals.append(str(val))
            vals.append("\n")
            lines.append("     ".join(vals))

        vars = ['NMat', 'NLay', 'CosAlfa', "\n"]
        lines.append("  ".join(vars))
        lines.append("   ".join([str(self.basic_information[var]) for var in
                                 vars[:-1]]))
        lines.append("\n")

        # Write block B: WATER FLOW INFORMATION
        lines.append(string.format(" BLOCK B: WATER FLOW INFORMATION ",
                                   fill="*", align="<", width=72))
        lines.append("\n")

        # Write BLOCK C: TIME INFORMATION
        lines.append(string.format(" BLOCK C: TIME INFORMATION ", fill="*",
                                   align="<", width=72))

        vars3 = [['dt', 'dtMin', 'dtMax', 'dMul', 'dMul2', 'ItMin', 'ItMax',
                  'MPL', "\n"], ["tInit", "tMax", "\n"],
                 ['lPrint', 'nPrintSteps', 'tPrintInterval', 'lEnter', "\n"]]
        for vars in vars3:
            lines.append(" ".join(vars))
            vals = []
            for var in vars[:-1]:
                val = self.time_information[var]
                if val is True:
                    vals.append("t")
                elif val is False:
                    vals.append("f")
                else:
                    vals.append(str(val))
            vals.append("\n")
            lines.append(" ".join(vals))

        lines.append("".join(["TPrint(1),TPrint(2),...,TPrint(MPL)\n"]))
        times = logspace(self.time_information["TPrint(1)"],
                         log10(self.time_information["TPrint(MPL)"]),
                         self.time_information["MPL"])
        lines.append(" ".join([str(time) for time in times]))
        lines.append("\n")

        # Write END statement
        lines.append(string.format(" END OF INPUT FILE 'SELECTOR.IN' ",
                                   fill="*", align="<", width=72))

        # Write the actual file
        with open(fname, "w") as file:
            file.writelines(lines)
        print("Succesfully wrote SELECTOR.IN")

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
            file.seek(0)  # Go back to start of file
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
            file.seek(0)  # Go back to start of file
            # Read the profile data into a Pandas DataFrame
            data = pd.read_csv(file, skiprows=start, skipfooter=2, index_col=0,
                               skipinitialspace=True, delim_whitespace=True)
            # Fix the header
            data.columns = [col + str(col1) for col, col1 in
                            data.iloc[0].items()]
            data = data.drop(data.index[0]).apply(pd.to_numeric)
        return data
