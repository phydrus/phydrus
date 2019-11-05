"""
This file contains the model class.
"""

import os
import subprocess

import numpy as np
import pandas as pd

from .plot import Plots
from .version import __version__


class Model:
    def __init__(self, exe_name, ws_name, name="model", description=None,
                 length_unit="m", time_unit="days", mass_units="mmol",
                 print_screen=False):
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
        print_screen: bool, optional
            Print the results to the screen during code execution.

        """
        # Store the hydrus executable and the project workspace
        self.exe_name = exe_name

        if not os.path.exists(ws_name):
            os.mkdir(ws_name)
            print("Directorty {} created".format(ws_name))

        self.ws_name = ws_name

        self.name = name
        self.description = description

        self.basic_information = {
            "iVer": "4",
            "Hed": "Created with Pydrus version {}".format(__version__),
            "LUnit": length_unit,
            "TUnit": time_unit,
            "MUnit": mass_units,
            "lWat": False,
            "lChem": False,
            "lTemp": False,
            "lSink": False,
            "lRoot": False,
            "lShort": True,
            "lWDep": False,
            "lScreen": print_screen,
            "AtmInf": False,
            "lEquil": True,
            "lInverse": False,
            "lSnow": False,
            "lHP1": False,
            "lMeteo": False,
            "lVapor": False,
            "lActRSU": False,
            "lFlux": False,
            "lIrrig": False,
            "NMat": 0,
            "NLay": 1,
            "CosAlfa": 1,
        }

        self.time_information = {
            "dt": 0.1,
            "dtMin": 0.0001,
            "dtMax": 0.5,
            "dMul": 1.3,
            "dMul2": 0.7,
            "ItMin": 3,
            "ItMax": 7,
            "MPL": None,  # Calculate automatically
            "tInit": 0.1,
            "tMax": 1,
            "lPrint": True,
            "nPrintSteps": 10,
            "tPrintInterval": 10,
            "lEnter": False,  # This should not be changed!
            "TPrint(1)": None,
            "TPrint(MPL)": None,
        }

        # The main processes to describe the simulation
        self.profile = None
        self.materials = None
        self.observations = []
        self.drains = None

        self.water_flow = None
        self.solute_transport = None
        self.heat_transport = None
        self.rootwater_uptake = None
        self.root_growth = None
        self.atmosphere_information = None
        self.atmosphere = None
        self.CO2Transport = None
        self.dual_porosity = None

        self.plots = Plots(ml=self)

    def add_profile(self, profile):
        """Method to add the soil profile to the model.

        """
        self.profile = profile
        self.basic_information["NLay"] = len(profile.loc[:, "Lay"].unique())

    def add_material(self, material):
        """Method to add a material to the model.

        Parameters
        ----------
        material: pandas.DataFrame
            Pandas Dataframe with the parameter names as columns and the
        values for each material as one row.

        Examples
        --------
        m = pd.DataFrame({1: [0.08, 0.3421, 0.03, 5, 1, -0.5]},
                 columns=["thr", "ths", "Alfa", "n" "Ks", "l"],
                 index=[1])
        ml.add_material(m)

        """
        if self.materials is None:
            self.materials = material
            self.basic_information["NMat"] = material.index.size
        else:
            self.materials = self.materials.append(material)
            self.basic_information["NMat"] += 1

    def add_drains(self):
        """Method to add a drain to the model.

        Returns
        -------

        """
        return NotImplementedError

    def add_observations(self, observations):
        """Method to add observation points.

        Parameters
        ----------
        observations: list of ints
            List of integers denoting the nodes to add a observation point
            to.

        """
        for obs in observations:
            nodes = self.profile.iloc[
                (self.profile['x'] - obs).abs().argsort()[:1]]
            node = nodes.index.item()
            self.observations.append(node)

    def add_waterflow(self, model=0, maxit=20, tolth=1e-4, tolh=0.1, ha=1e-3,
                      hb=1e3, topinf=False, botinf=False, kodtop=-1, kodbot=1,
                      linitw=True, free_drainage=False, seepage_face=False,
                      qdrain=False, hseep=0, rtop=0, rbot=0, rroot=0,
                      qgwlf=False, gw_level=None, aqh=None, bqh=None,
                      hysteresis=0, ikappa=-1, wlayer=False):
        """Method to add a water_flow module to the model.

        Parameters
        ----------
        model: int, optional
            Soil hydraulic properties model:
            0 = van Genuchten"s [1980] model with 6 parameters.
            1 = modified van Genuchten"s model  with 10 parameters [Vogel
            and Císlerová, 1988].
            2 = Brooks and Corey"s [1964] model with 6 parameters.
            3 = van Genuchten"s [1980] model with air-entry value of -2 cm
            and with 6 parameters.
            4 = Kosugi’s [1996] model with 6 parameters.
            5 = dual porosity model of Durner [1994] with 9 parameters.
            6 = dual-porosity system with transfer proportional to the
            effective saturation (9 parameters).
            7 = dual-porosity system with transfer proportional to the
            pressure head (11 parameters).
            9 = dual-permeability system with transfer proportional to the
            pressure head (17 parameters)
            !!! model>3 options are not available with the major ion chemistry.
            module.
        maxit: int, optional
            Maximum number of iterations allowed during any time step.
        tolth: float, optional
            Absolute water content tolerance for nodes in the unsaturated part
            of the flow region [-]. TolTh represents the maximum desired
            absolute change in the value of the water content, θ, between
            two successive iterations during a particular time step.
        tolh: float, optional
            Absolute pressure head tolerance for nodes in the saturated part of
             the flow region [L] (its recommended value is 0.1 cm). TolH
             represents the maximum desired absolute change in the value of the
             pressure head, h, between two successive iterations during a
             particular time step.
        ha: float, optional
            Absolute value of the upper limit [L] of the pressure head interval
            below which a table of hydraulic properties will be generated
            internally for each material.
        hb: float, optional
            Absolute value of the lower limit [L] of the pressure head interval
             for which a table of hydraulic properties will be generated
             internally for each material.
        topinf: bool, optional
            Set to True if the surface boundary condition is time-dependent.
        botinf: bool, optional
            Set to True if the bottom boundary condition is time-dependent.
        kodtop: int, optional
            Code specifying type of boundary condition (BC) for water flow at
            the surface. 1 for Dirichlet BC and -1 for Neumann BC. Set to 0
            when a prescribed BC can change from Dirichlet BC to Neumann BC
            and vice versa.
        kodbot: int, optional
            Code specifying type of boundary condition for water flow at the
            bottom of the profile. Set to -1 for a Dirichlet BC and to 1 for
            a Neumann BC. In case of a seepage face or free drainage BC set
            KodBot=-1.
        linitw: bool, optional
            Set to True if the initial condition is given in terms of the
            water content. Set to False if given in terms of pressure head.
        free_drainage: bool, optional
            True if free drainage is to be considered as bottom BC.
        seepage_face: bool, optional
            True if seepage face is to be considered as the bottom BC.
        qdrain: bool, optional
            True if flow to horizontal drains is considered as bottom BC.
        hseep: float, optional
            Pressure head (i.e., 0) that initiates flow over the seepage face
            bottom boundary.
        rtop: float, optional
            Prescribed top flux [LT-1] (in case of a Dirichlet BC set this
            variable equal to zero).
        rbot: float, optional
            Prescribed bottom flux [LT-1] (in case of a Dirichlet BC set this
            variable equal to zero).
        rroot: float, optional
            Prescribed potential transpiration rate [LT-1] (if no transpiration
            occurs or if transpiration is variable in time set this variable
            equal to zero).
        qgwlf: bool, optional
            Set to True if the discharge-groundwater level relationship q(
            GWL) is applied as bottom boundary condition.
        gw_level: float, optional
            Reference position of the groundwater table (e.g., the
            x-coordinate of the soil surface).
        aqh: float, optional
            Value of the parameter Aqh [LT-1] in the q(GWL)-relationship.
        bqh: float, optional
            Value of the parameter Bqh [L-1] in the q(GWL)-relationship.
        hysteresis: int, optional
            Hysteresis in the soil hydraulic properties:
            0 = No hysteresis.
            1 = Hysteresis in the retention curve only.
            2 = Hysteresis in both the retention and hydraulic conductivity
            functions.
            3 = Hysteresis using Robert Lenhard’s model [Lenhard et al.,
            1991; Lenhard and Parker, 1992]. (Not available with major ion
            chemistry module.)
        ikappa: int, optional
            Set to -1 if the initial condition is to be calculated from the
            main drying branch. Set to 1 if the initial condition is to be
            calculated from the main wetting branch.
        wlayer: bool: optional
            Set to True if water can accumulate at the surface with zero
            surface runoff.

        """
        # If qgwlf is user as bottom boundary condition
        if qgwlf:
            for var in [gw_level, aqh, bqh]:
                if var is None:
                    raise TypeError("When the groundwater level is used as "
                                    "bottom boundary condition, the keyword "
                                    "{} needs to be provided".format(var))
        if self.water_flow is None:
            if topinf:
                # In thecase of 'Atmospheric BC' set KodTop=-1.
                kodtop = -1
            if free_drainage:
                # In case of a seepage face or free drainage BC set KodBot=-1.
                kodbot = -1
            if seepage_face:
                # In case of a seepage face or free drainage BC set KodBot=-1.
                kodbot = -1

            self.water_flow = {
                "MaxIt": maxit,  # Maximum No. of Iterations
                "TolTh": tolth,  # [-]
                "TolH": tolh,  # [L], default is 0.1 cm
                "TopInf": topinf,
                "WLayer": wlayer,
                "KodTop": kodtop,  # Depends on boundary condition
                "lInitW": linitw,
                "BotInf": botinf,
                "qGWLF": qgwlf,
                "FreeD": free_drainage,
                "SeepF": seepage_face,
                "KodBot": kodbot,  # Depends on boundary condition
                "qDrain": qdrain,
                "hSeep": hseep,  # [L]
                "rTop": rtop,  # [LT-1]
                "rBot": rbot,  # [LT-1]
                "rRoot": rroot,  # [LT-1]
                "GWL0L": gw_level,  # [L]
                "Aqh": aqh,  # [LT-1]
                "Bqh": bqh,  # [L-1]
                "ha": ha,  # [L]
                "hb": hb,  # [L]
                "iModel": model,
                "iHyst": hysteresis,
                "iKappa": ikappa,
            }
        else:
            raise Warning("Water flow was already provided. Please delete "
                          "the old information first.")

        self.basic_information["lWat"] = True

    def add_atmosphere(self, atmosphere, ldailyvar=False, lsinusvar=False,
                       llai=False, lbccycles=False, linterc=False,
                       rextinct=0.463, hcrits=1e30, tatm=0, prec=0, rsoil=0,
                       rroot=0, hcrita=1e5, rb=0, hb=0, ht=0, ttop=0, tbot=0,
                       ampl=0):
        """Method to add the atmospheric boundary condition to the model.

        Parameters
        ----------
        atmosphere: Pandas.DataFrame
            Pandas DataFrame with at least the following columns: tAtm,
            Prec, rSoil, rRoot, hCritA, rB, hB, hT, tTop, tBot, and Ampl.
        ldailyvar: bool, optional
            True if HYDRUS-1D is to generate daily variations in evaporation
            and transpiration. False otherwise.
        lsinusvar: bool, optional
            True if HYDRUS-1D is to generate sinusoidal variations in
            precipitation. False otherwise.
        llai: bool, optional
            Boolean indicating that potential evapotranspiration is
            to be divided into potential evaporation and potential
            transpiration using eq. (2.75) of the manual.
        lbccycles: bool, optional
            ???
        linterc: bool, optional
            ???
        rextinct: float, optional
            A constant for the radiation extinction by the canopy
            (rExtinct=0.463) [-]. only used when lLai is True.
        hcrits: float, optional
            Maximum allowed pressure head at the soil surface [L]. Default is
            1e+30.

        Notes
        -----
        The index of the atmosphere DataFrame should be a RangeIndex with
        integers.

        """
        if self.atmosphere_information is None:
            self.atmosphere_information = {
                "lDailyVar": ldailyvar,
                "lSinusVar": lsinusvar,
                "lLai": llai,
                "lBCCycles": lbccycles,
                "lInterc": linterc,
                "lExtinct": rextinct,
                "hCritS": hcrits,
            }
        else:
            raise Warning("Atmospheric information was already provided. "
                          "Please delete the old information first through "
                          "ml.del_atmosphere().")

        data = {"tAtm": tatm, "Prec": prec, "rSoil": rsoil, "rRoot": rroot,
                "hCritA": hcrita, "rB": rb, "hB": hb, "ht": ht, "tTop": ttop,
                "tBot": tbot, "Ampl": ampl}

        self.atmosphere = pd.DataFrame(data=data, index=atmosphere.index)
        self.atmosphere.update(atmosphere)

        # Enable atmosphere module
        self.basic_information["AtmInf"] = True
        self.water_flow["TopInf"] = True
        self.water_flow["KodTop"] = -1

    def add_rootwater_uptake(self, model=0, crootmax=0, omegac=0.5, p0=-10,
                             p2h=-200, p2l=-800, p3=-8000, r2h=0.5, r2l=0.1,
                             poptm=None):
        """Method to add rootwater update modeule to the model.

        Parameters
        ----------
        model: int, optional
            Type of root water uptake stress response function. 0 = Feddes
            et al. [1978]. 1 = S-shaped, van Genuchten [1987]
        crootmax: float, optional
            Maximum allowed concentration in the root solute uptake term for
            the first solute [ML-3]. When the nodal concentration is lower than
            cRootMax, all solute is taken up. When the nodal concentration
            is higher than cRootMax, additional solute stays behind.
        omegac: float, optional
            Maximum allowed concentration in the root solute uptake term for
            the last solute [ML-3].
        p0: float, optional
            Only used if model=0. Value of the pressure head, h1, below
            which roots start to extract water from the soil.
        p2h: float, optional
            Only used if model=0. Value of the limiting pressure head, h3,
            below which the roots cannot extract water at the maximum rate
            (assuming a potential transpiration rate of r2H).
        p2l: float, optional
            Only used if model=0. As above, but for a potential transpiration
            rate of r2L.
        p3: float, optional
            Only used if model=0. Value of the pressure head, h4, below which
            root water uptake ceases (usually equal to the wilting point).
        r2h: float, optional
            Only used if model=0. Potential transpiration rate [LT-1]
            (currently set at 0.5 cm/day).
        r2l: float, optional
            Only used if model=0. Potential transpiration rate [LT-1]
            (currently set at 0.1 cm/day).
        poptm: list, optional
            Value of the pressure head, h2, below which roots start to extract
            water at the maximum possible rate. The length of poptm should
            equal the No. of materials.

        """
        if model == 1:
            raise NotImplementedError("Sorry, the S-shaped model has not been "
                                      "implemented yet!")

            # Number of pressure heads should equal the number of materials.
        if poptm:
            if len(poptm) != self.basic_information["NMat"]:
                raise KeyError("Length of pressure heads poptm does not "
                               "equal the number of materials!")

        if self.rootwater_uptake is None:
            self.rootwater_uptake = {
                "iMoSink": model,
                "cRootMax": crootmax,
                "OmegaC": omegac,
                "POptm": poptm,
                "P0": p0,
                "P2H": p2h,
                "P2L": p2l,
                "P3": p3,
                "r2H": r2h,
                "r2L": r2l,
            }
            self.basic_information["lSink"] = True
        else:
            msg = "Root water uptake model is already present in the model." \
                  " Remove the old root water uptake model first using " \
                  "ml.del_root_water_uptake()"
            raise InterruptedError(msg)

    def add_root_growth(self, irootin=0, ngrowth=None, tgrowth=None,
                        rootdepth=None, irfak=None, trmin=None, trmed=None,
                        trharv=None, xrmin=None, xrmed=None, xrmax=None,
                        trperiod=None):
        """Method to add root growth to the model.

        Parameters
        ----------
        irootin: int
            0 = (default) The root depth is specified together with other
            time-variable boundary condition, such as atmospheric fluxes.
            1 = the root depth is given in a table
            2 = the root depth is calculated using the growth function.
        ngrowth: int
            Number of data points in the root depth table. Only used when
            irootin = 1.
        tgrowth: float, optional
            Days. Only used when irootin = 1.
        rootdepth: list of float, optional
            Rooting depth [L]. List has a length of ngrowth. Only used when
            irootin = 1.
        irfak: int, optional
            Method to calculate the root growth factor, r. Only used when
            irootin = 2.
            0= the root growth factor is calculated from given data [xRMed,
            tRMed].
            1 = the root growth factor is calculated based on the assumption
            that 50% of the rooting depth, (xRMax+xRMin)/2., is reached at
            the midpoint of the growing season, (tRMin+tRHarv)/2.
        trmin: float, optional
            Initial time of the root growth period [T]. Only used when
            irootin = 2.
        trmed: float, optional
            Time of known rooting depth (set equal to zero if iRFak=1) [T].
            Only used when irootin = 2.
        trharv: float, optional
            Time at the end of the root water uptake period [T]. Only used when
            irootin = 2.
        xrmin: float, optional
            Initial value of the rooting depth at the beginning of the
            growth period (recommended value = 1 cm) [L]. Only used when
            irootin = 2.
        xrmed: float, optional
            Value of known rooting depth (set equal to zero if iRFak=1) [L].
            Only used when irootin = 2.
        xrmax: float, optional
            Maximum rooting depth, which may be reached at infinite time [L].
            Only used when irootin = 2.
        trperiod: float, optional
            Time period at which the growth function repeats itself. Only
            used when irootin = 2.

        """

        # Store the root growth information depending on the model.
        if irootin is 0:
            root_growth = {
                "iRootIn": irootin
            }
        elif irootin is 1:
            root_growth = {
                "iRootIn": irootin,
                "nGrowht": ngrowth,
                "tGrwoth": tgrowth,
                "RootDepth": rootdepth
            }
        elif irootin is 2:
            root_growth = {
                "iRootIn": irootin,
                "iRFak": irfak,
                "tRMin": trmin,
                "tRMed": trmed,
                "tRHarv": trharv,
                "xRMin": xrmin,
                "xRMed": xrmed,
                "xRMax": xrmax,
                "tRPeriod": trperiod
            }
        else:
            raise Warning("Option %s for irootin is not support in Hydrus."
                          % irootin)

        if self.root_growth is None:
            self.root_growth = root_growth
            self.basic_information["lRoot"] = True
        else:
            raise Warning("Root growth model already exists. Please delete "
                          "the old root growth model first using "
                          "ml.del_root_growth().")

    def simulate(self):
        """Method to call the Hydrus-1D executable.

        """
        # Run Hydrus executable.
        cmd = [self.exe_name, self.ws_name, "-1"]
        result = subprocess.run(cmd)

        return result

    def get_print_times(self):
        """Method to get the print times for the simulation.

        Returns
        -------

        """
        if self.time_information["TPrint(1)"] is None:
            tmin = self.time_information["tInit"] + 1
        else:
            tmin = self.time_information["TPrint(1)"]

        if self.time_information["TPrint(MPL)"] is None:
            tmax = self.time_information["tMax"]
        else:
            tmax = self.time_information["TPrint(MPL)"]

        times = range(tmin, tmax)

        return times

    def write_files(self):
        self.write_selector()

        # Write Block Block I - Atmospheric information
        if self.basic_information["AtmInf"]:
            self.write_atmosphere()

        self.write_profile()

    def write_selector(self, fname="SELECTOR.IN"):
        """Write the selector.in file.

        """
        # Create Header string
        string = "*** {:{fill}{align}{width}}\n"

        lines = ["Pcp_File_Version={}\n".format(
            self.basic_information["iVer"]),
            string.format("BLOCK A: BASIC INFORMATION ", fill="*",
                          align="<", width=72),
            "{}\n".format(self.basic_information["Hed"]),
            "{}\n".format(self.description),
            "LUnit  TUnit  MUnit  (indicated units are obligatory "
            "for all input data)\n",
            "{}\n".format(self.basic_information["LUnit"]),
            "{}\n".format(self.basic_information["TUnit"]),
            "{}\n".format(self.basic_information["MUnit"])]

        # Write block A: BASIC INFORMATION

        vars_list = [["lWat", "lChem", "lTemp", "lSink", "lRoot", "lShort",
                      "lWDep", "lScreen", "AtmInf", "lEquil", "lInverse",
                      "\n"],
                     ["lSnow", "lHP1", "lMeteo", "lVapor", "lActRSU", "lFlux",
                      "lIrrig", "\n"]]

        for variables in vars_list:
            lines.append("  ".join(variables))
            values = []
            for var in variables[:-1]:
                val = self.basic_information[var]
                if val is True:
                    values.append("t")
                elif val is False:
                    values.append("f")
                else:
                    values.append(str(val))
            values.append("\n")
            lines.append("     ".join(values))

        variables = ["NMat", "NLay", "CosAlfa", "\n"]
        lines.append("  ".join(variables))
        lines.append("   ".join([str(self.basic_information[var]) for var in
                                 variables[:-1]]))
        lines.append("\n")

        # Write block B: WATER FLOW INFORMATION
        lines.append(string.format("BLOCK B: WATER FLOW INFORMATION ",
                                   fill="*", align="<", width=72))
        lines.append("MaxIt  TolTh  TolH   (maximum number of iterations and "
                     "tolerances)\n")
        variables = ["MaxIt", "TolTh", "TolH"]
        lines.append(
            "   ".join([str(self.water_flow[var]) for var in variables]))
        lines.append("\n")

        vars_list = [["TopInf", "WLayer", "KodTop", "lInitW", "\n"],
                     ["BotInf", "qGWLF", "FreeD", "SeepF", "KodBot", "qDrain",
                      "hSeep", "\n"]]

        #        if (self.water_flow["KodTop"] >= 0) or \
        #                (self.water_flow["KodBot"] >= 0):
        #            vars_list.append(["rTop", "rBot", "rRoot", "\n"])
        # Records 8a and 9a are provided only when lower or upper boundary
        # conditions are independent of time and at least one of them is a 
        # Neumann BC.

        if self.water_flow["qGWLF"]:
            vars_list.append(["GWL0L", "Aqh", "Bqh", "\n"])

        vars_list.append(["ha", "hb", "\n"])
        vars_list.append(["iModel", "iHyst", "\n"])

        if self.water_flow["iHyst"] > 0:
            vars_list.append(["iKappa", "\n"])

        for variables in vars_list:
            lines.append("  ".join(variables))
            values = []
            for var in variables[:-1]:
                val = self.water_flow[var]
                if val is True:
                    values.append("t")
                elif val is False:
                    values.append("f")
                else:
                    values.append(str(val))
            values.append("\n")
            lines.append("     ".join(values))

        if self.drains:
            raise NotImplementedError

        # Write the material parameters
        lines.append(self.materials.to_string(index=False))
        lines.append("\n")

        # Write BLOCK C: TIME INFORMATION
        lines.append(string.format("BLOCK C: TIME INFORMATION ", fill="*",
                                   align="<", width=72))

        times = self.get_print_times()

        #        self.time_information["MPL"] = len(times)

        vars_list = [
            ["dt", "dtMin", "dtMax", "dMul", "dMul2", "ItMin", "ItMax",
             "MPL", "\n"], ["tInit", "tMax", "\n"],
            ["lPrint", "nPrintSteps", "tPrintInterval", "lEnter", "\n"]]
        for variables in vars_list:
            lines.append(" ".join(variables))
            values = []
            for var in variables[:-1]:
                val = self.time_information[var]
                if val is True:
                    values.append("t")
                elif val is False:
                    values.append("f")
                else:
                    values.append(str(val))
            values.append("\n")
            lines.append(" ".join(values))

        lines.append("TPrint(1),TPrint(2),...,TPrint(MPL)\n")
        #        for i in range(int(len(times) / 6) + 1):
        #            lines.append(
        #                " ".join([str(time) for time in times[i * 6:i * 6 + 6]]))
        #            lines.append("\n")
        printtimes = np.linspace(start=self.time_information["tInit"],
                                 stop=self.time_information["tMax"],
                                 num=self.time_information["MPL"] + 1)[1:]
        roundtimes = printtimes.round(6)
        strtimes = "  ".join(str(e) for e in roundtimes)
        self.strtimes = strtimes
        lines.append(strtimes)
        lines.append("\n")

        # Write BLOCK D: Root Growth Information
        if self.basic_information["lRoot"]:
            raise NotImplementedError

        # Write Block E - Heat transport information
        if self.basic_information["lTemp"]:
            raise NotImplementedError

        # Write Block F - Solute transport information
        if self.basic_information["lChem"]:
            raise NotImplementedError

        # Write Block G - Root water uptake information
        if self.basic_information["lSink"]:
            lines.append(
                string.format("BLOCK G: ROOT WATER UPTAKE INFORMATION ",
                              fill="*", align="<", width=72))
            vars_list = [["iMoSink", "cRootMax", "OmegaC", "\n"]]

            if self.rootwater_uptake["iMoSink"] is 0:
                vars_list.append(
                    ["P0", "P2H", "P2L", "P3", "r2H", "r2L", "\n"])

            for variables in vars_list:
                lines.append(" ".join(variables))
                values = []
                for var in variables[:-1]:
                    val = self.rootwater_uptake[var]
                    if var:
                        if val is True:
                            values.append("t")
                        elif val is False:
                            values.append("f")
                        else:
                            values.append(str(val))
                lines.append(" ".join(values))
                lines.append("\n")

            lines.append("POptm(1),POptm(2),...,POptm(NMat)\n")
            lines.append(" ".join(str(p) for p in self.rootwater_uptake[
                "POptm"]))
            lines.append("\n")

        # Write Block H - Nodal information

        # Write Block J - Inverse solution information
        if self.basic_information["lInverse"]:
            raise NotImplementedError("The inverse modeling module from "
                                      "Hydrus-1D will not be supported. "
                                      "Python packages are used for this.")

        # Write Block K – Carbon dioxide transport information

        # Write Block L – Major ion chemistry information
        if self.basic_information["lChem"]:
            raise NotImplementedError

        # Write Block M – Meteorological information
        if self.basic_information["lMeteo"]:
            raise NotImplementedError

        # Write END statement
        lines.append(string.format(" END OF INPUT FILE SELECTOR.IN ",
                                   fill="*", align="<", width=72))

        # Write the actual file
        fname = os.path.join(self.ws_name, fname)
        with open(fname, "w") as file:
            file.writelines(lines)
        print("Succesfully wrote {}".format(fname))

    def write_atmosphere(self, fname="ATMOSPH.IN"):
        """Method to write the ATMOSPH.IN file

        """
        # 1 Write Header information
        lines = ["Pcp_File_Version={}\n".format(
            self.basic_information["iVer"]),
            "*** BLOCK I: ATMOSPHERIC INFORMATION  "
            "**********************************\nMaxAL "
            "(MaxAL = number of atmospheric data-records)\n"]

        # Print some values
        nrow = self.atmosphere.index.size
        lines.append("{}\n".format(nrow))

        vars5 = ["lDailyVar", "lSinusVar", "lLai", "lBCCycles", "lInterc",
                 "\n"]

        lines.append(" ".join(vars5))
        vals = []
        for var in vars5[:-1]:
            val = self.atmosphere_information[var]
            if var:
                if val is True:
                    vals.append("t")
                elif val is False:
                    vals.append("f")
                else:
                    vals.append(str(val))
        lines.append(" ".join(vals))
        lines.append("\n")

        lines.append("hCritS (max. allowed pressure head at the soil surface)")
        lines.append("\n{}\n".format(self.atmosphere_information["hCritS"]))

        lines.append(self.atmosphere.to_string(index=False))
        lines.append("\n")
        lines.append("end*** END OF INPUT FILE ATMOSPH.IN "
                     "**********************************\n")
        # Write the actual file
        fname = os.path.join(self.ws_name, fname)
        with open(fname, "w") as file:
            file.writelines(lines)
        print("Succesfully wrote {}".format(fname))

    def write_profile(self, fname="PROFILE.DAT", ws=""):
        """Method to write the profile.dat file

        """
        # 1 Write Header information
        lines = ["Pcp_File_Version=4\n"]
        lines.append("0\n")  # TODO Figure out what these lines do.

        # Print some values
        nrow = self.profile.index.size
        ii = 0  # TODO Figure this out
        ns = 0  # TODO Number of solutes
        lines.append("{} {} {} ".format(nrow, ii, ns))

        # 2. Write the profile data
        lines.append(self.profile.to_string())
        lines.append("\n")

        # 3. Write observation points
        lines.append(str(len(self.observations)) + os.linesep)
        lines.append("".join(["   {}".format(i) for i in self.observations]))

        # Write the actual file
        fname = os.path.join(self.ws_name, fname)
        with open(fname, "w") as file:
            file.writelines(lines)
        print("Succesfully wrote {}".format(fname))

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
        if not os.path.exists(path):
            raise FileNotFoundError(
                "File {} has not been found.".format(path))

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

    def read_tlevel(self, fname="T_LEVEL.OUT", use_cols=None):
        """Method to read the T_LEVEL.OUT output file.

        Parameters
        ----------
        fname: str, optional
            String with the name of the t_level out file. default is
            "T_LEVEL.OUT".
        use_cols: list of str optional
            List with the names of the columns to import. By default
            only the real fluxes are imported and not the cumulative
            fluxes. Options are: "rTop", "rRoot", "vTop", "vRoot", "vBot",
            "sum(rTop)", "sum(rRoot)", "sum(vTop)", "sum(vRoot)", "sum(vBot)",
            "hTop", "hRoot", "hBot", "RunOff", "sum(RunOff)", "Volume",
            "sum(Infil)", "sum(Evap)", "TLevel", "Cum(WTrans)", "SnowLayer".

        Returns
        -------
        data: pandas.DataFrame
            Pandas with the t_level data

        """
        path = os.path.join(self.ws_name, fname)
        if not os.path.exists(path):
            raise FileNotFoundError(
                "File {} has not been found.".format(path))

        if use_cols is None:
            use_cols = ["Time", "rTop", "rRoot", "vTop", "vRoot",
                        "vBot", "sum(rTop)", "sum(rRoot)", "sum(vTop)",
                        "sum(vRoot)", "sum(vBot)", "hTop", "hRoot", "hBot",
                        "RunOff", "Volume", ]

            if self.water_flow["iModel"] > 4:
                use_cols.append("Cum(WTrans)")
            if self.basic_information["lSnow"]:
                use_cols.append("SnowLayer")

        with open(path) as file:
            # Find the starting line to read the profile
            for start, line in enumerate(file.readlines(1000)):
                if "rTop" in line:
                    break
            file.seek(0)  # Go back to start of file

            skiprows = list(range(start))
            skiprows.append(start + 1)

            # Read the profile data into a Pandas DataFrame
            data = pd.read_csv(file, skiprows=skiprows, skipfooter=2,
                               index_col="Time", skipinitialspace=True,
                               delim_whitespace=True, usecols=use_cols,
                               engine="python")

        return data

    def read_nod_inf(self, fname="NOD_INF.OUT", times=None):
        """Method to read the NOD_INF.OUT output file.

        Parameters
        ----------
        fname: str, optional
            String with the name of the NOD_INF out file. default is
            "NOD_INF.OUT".
        times = int, optional
            Create a DataFrame with nodal values of the pressure head, 
            the water content, the solution and sorbed concentrations, and 
            temperature, etc, at the time "times". default is None.
        Returns
        -------
        data: dict
            Dictionary with the time as a key and a Pandas DataFrame as a
            value.

        """
        path = os.path.join(self.ws_name, fname)
        if not os.path.exists(path):
            raise FileNotFoundError(
                "File {} has not been found.".format(path))

        use_times = []
        start = []
        end = []
        data = {}

        with open(path) as file:
            # Find the starting times to read the information
            for i, line in enumerate(file.readlines()):
                if "Time" in line and not "Date" in line:
                    time = float(line.replace(" ", "").split(":")[
                                     1].replace("\n", ""))
                    use_times.append(time)
                elif "Node" in line:
                    start.append(i)
                elif "end" in line:
                    end.append(i)

            if times is None:
                times = use_times
                # Read the profile data into a Pandas DataFrame
                for start, end, time in zip(start, end, use_times):
                    if time in times:
                        file.seek(0)  # Go back to start of file
                        df = pd.read_csv(file, skiprows=start, index_col=0,
                                         skipinitialspace=True,
                                         delim_whitespace=True,
                                         nrows=end - start - 2)
                        #
                        # # Fix the header
                        df = df.drop(df.index[0]).apply(pd.to_numeric)
                        df.index = pd.to_numeric(df.index)

                        data[time] = df
            elif not isinstance(times, list):
                time = use_times.index(times)
                file.seek(0)  # Go back to start of file
                df = pd.read_csv(file, skiprows=start[time], index_col=0,
                                 skipinitialspace=True,
                                 delim_whitespace=True,
                                 nrows=end[time] - start[time] - 2)
                # # Fix the header
                df = df.drop(df.index[0]).apply(pd.to_numeric)
                df.index = pd.to_numeric(df.index)
                data = df
        if len(data) is 1:
            return next(iter(data.values()))
        else:
            return data

    def read_run_inf(self, fname="RUN_INF.OUT", use_cols=None):
        """Method to read the RUN_INF.OUT output file.

        Parameters
        ----------
        fname: str, optional
            String with the name of the run_inf out file. default is
            "RUN_INF.OUT".
        use_cols: list of str optional
            List with the names of the columns to import. By default:
			"TLevel", "Time", "dt", "Iter", "ItCum",
            "KodT", "KodB", "Convergency".

        Returns
        -------
        data: pandas.DataFrame
            Pandas with the run_inf data

        """
        path = os.path.join(self.ws_name, fname)
        if not os.path.exists(path):
            raise FileNotFoundError(
                "File {} has not been found.".format(path))

        if use_cols is None:
            use_cols = ["TLevel", "Time", "dt", "Iter", "ItCum", "KodT",
                        "KodB", "Convergency", ]
            if not self.solute_transport == None:
                use_cols.append("IterC")

        with open(path) as file:
            # Find the starting line to read the profile
            for start, line in enumerate(file.readlines(1000)):
                if "TLevel" in line:
                    break
            file.seek(0)  # Go back to start of file

            skiprows = list(range(start))
            skiprows.append(start + 1)

            # Read the profile data into a Pandas DataFrame
            data = pd.read_csv(file, skiprows=skiprows, skipfooter=1,
                               index_col="TLevel", skipinitialspace=True,
                               delim_whitespace=True, usecols=use_cols,
                               engine="python")

        return data

    def read_balance(self, fname="BALANCE.OUT", use_cols=None):
        """Method to read the BALANCE.OUT output file.

        Parameters
        ----------
        fname: str, optional
            String with the name of the run_inf out file. default is
            "BALANCE.OUT".
        use_cols: list of str optional
            List with the names of the columns to import. By default:
            "Area","W-volume","In-flow","h Mean","Top Flux", 
            "Bot Flux","WatBalT","WatBalR".

        Returns
        -------
        data: pandas.DataFrame
            Pandas with the balance data

        """
        path = os.path.join(self.ws_name, fname)
        if not os.path.exists(path):
            raise FileNotFoundError(
                "File {} has not been found.".format(path))

        if use_cols is None:
            use_cols = ["Area", "W-volume", "In-flow", "h Mean", "Top Flux",
                        "Bot Flux", "WatBalT", "WatBalR"]
            if not self.solute_transport == None:
                use_cols.append("ConcVol", "ConcVolIm", "cMean", "CncBalT",
                                "CncBalR")

            if not self.heat_transport == None:
                use_cols.append("TVol", "TMean")

            if not self.CO2Transport == None:
                use_cols.append("COVol", "COMean", "CO2BalT", "CncBalT")

            if not self.dual_porosity == None:
                use_cols.append("W-VolumeI", "cMeanIm")

        List = open(path).readlines()
        List_copy = List.copy()
        use_times = []
        start = []
        end = [16]
        for i, line in enumerate(List_copy):
            for x in use_cols:
                if x in line:
                    line1 = x
                    line2 = line.replace("\n", "").split(" ")[-1]
                    line3 = line.replace("  ", " ").split(" ")[-2]
                    List_copy[i] = [line1, line2, line3]

            if "Time" in line and not "Date" in line:
                time = float(
                    line.replace(" ", "").split("]")[1].replace("\n", ""))
                use_times.append(time)
            if "Area" in line:
                start.append(i)
            if "WatBalR" in line:
                end.append(i + 1)
            if "Sub-region" in line:
                subregions = \
                    line.replace("  ", " ").replace("\n", "").split(" ")[-1]
        data = {}
        for s, e, time in zip(start, end, use_times):
            List_copy1 = List_copy[s:e]
            df = pd.DataFrame(List_copy1).set_index(0).T
            indexc = {}
            for x in range(int(subregions) + 1):
                indexc[x + 1] = x
                df = df.rename(index=indexc)
            data[time] = df

        return data

    def read_obs_node(self, fname="OBS_NODE.OUT", times=None):
        """Method to read the OBS_NODE.OUT output file.

        Parameters
        ----------
        fname: str, optional
            String with the name of the OBS_NODE out file. default is
            "OBS_NODE.OUT".

        Returns
        -------
        data: dict
            Dictionary with the node as a key and a Pandas DataFrame as a
            value.

        """
        path = os.path.join(self.ws_name, fname)
        if not os.path.exists(path):
            raise FileNotFoundError(
                "File {} has not been found.".format(path))

        data = {}
        with open(path) as file:
            # Find the starting times to read the information
            for i, line in enumerate(file.readlines()):
                if "time" in line:
                    start = i
                elif "end" in line:
                    end = i
            for i, nod in enumerate(self.observations):
                file.seek(0)
                if i == 0:
                    usecols = ["time", "h", "theta", "Temp"]
                else:
                    usecols = ["time", "h." + str(1), "theta." + str(1),
                               "Temp." + str(1)]
                df = pd.read_csv(file, skiprows=start, index_col=0,
                                 skipinitialspace=True,
                                 delim_whitespace=True,
                                 nrows=end - start - 1,
                                 usecols=usecols)
                df.columns = ["h", "theta", "Temp"]
                data[nod] = df

        return data

    def read_I_check(self, fname="I_CHECK.OUT", times=None):
        """Method to read the I_CHECK.OUT output file.

        Parameters
        ----------
        fname: str, optional
            String with the name of the I_Check out file. default is
            "I_Check.OUT".

        Returns
        -------
        data: dict
            Dictionary with the node as a key and a Pandas DataFrame as a
            value.

        """
        usecols = ["theta", "h", "log_h", "C", "K", "log_K", "S", "Kv"]
        path = os.path.join(self.ws_name, fname)
        if not os.path.exists(path):
            raise FileNotFoundError(
                "File {} has not been found.".format(path))

        with open(path) as file:
            # Find the starting times to read the information
            for i, line in enumerate(file.readlines()):
                if "theta" in line:
                    start = i
                elif "end" in line:
                    end = i

            # Read the profile data into a Pandas DataFrame
            file.seek(0)
            df = pd.read_csv(file, skiprows=start + 1, index_col=None,
                             header=None, names=usecols,
                             skipinitialspace=True,
                             delim_whitespace=True,
                             nrows=int(end) - int(start) - 2)
        return df
