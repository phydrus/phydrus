"""
This file contains the model class.
"""


class Model:
    def __init__(self, name="model", description=None, length_unit="m",
                 time_unit="days", mass_units="mmol"):
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
        # 1. Check if all files exist.

        # 2. Run Hydrus executable.

        # 3. Read output and check for success.

        pass

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
