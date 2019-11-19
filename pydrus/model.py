"""
This file contains the model class.
"""

import os
from subprocess import run

import numpy as np
import pandas as pd
from pandas import DataFrame

from .plot import Plots
from .read import read_profile, read_nod_inf, read_run_inf, read_tlevel, \
    read_balance, read_i_check, read_obs_node, read_solute, read_alevel
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
        if not os.path.exists(exe_name):
            raise Warning("Path to the Hydrus-1D executable seems incorrect, "
                          "please check the path to the executable.")
        else:
            self.exe_name = exe_name

        if not os.path.exists(ws_name):
            os.mkdir(ws_name)
            print("Directory {} created".format(ws_name))

        self.ws_name = ws_name

        self.name = name
        self.description = description

        # The main processes to describe the simulation
        self.profile = None
        self.materials = None
        self.observations = []
        self.drains = None
        self.times = None

        self.water_flow = None
        self.solute_transport = None
        self.solutes = None
        self.heat_transport = None

        self.root_uptake = None
        self.root_growth = None

        self.atmosphere_info = None
        self.atmosphere = None

        self.basic_info = {
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
            "CosAlfa": 1,
        }

        self.time_info = {
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

        self.plots = Plots(ml=self)

    @property
    def n_materials(self):
        if self.materials is None:
            return 0
        else:
            return self.materials.index.size

    @property
    def n_solutes(self):
        if self.solutes is None:
            return 0
        else:
            return self.solutes.index.size

    @property
    def n_layers(self):
        if self.profile is None:
            return 0
        else:
            return len(self.profile.loc[:, "Lay"].unique())

    def add_profile(self, profile):
        """Method to add the soil profile to the model.

        """
        self.profile = profile

    def add_material(self, material):
        """Method to add a material to the model.

        Parameters
        ----------
        material: pandas.DataFrame
            Pandas DataFrame with the parameter names as columns and the
            values for each material as one row. The index for each is the
            reference number for each material and must be unique. The
            number of columns depends on the water flow model that has been
            chosen.

        Examples
        --------
        m = pd.DataFrame({1: [0.08, 0.3421, 0.03, 5, 1, -0.5]}, index=[1],
                         columns=["thr", "ths", "Alfa", "n" "Ks", "l"])
        ml.add_material(m)

        """
        if self.materials is None:
            raise Warning("The water flow module has to be chosen before "
                          "adding the soil materials. Use ml.add_water_flow("
                          ") to add water flow to the model")
        elif material.columns.size != self.materials.columns.size:
            raise TypeError("the number of parameters (columns) describing "
                            "the material does not match the water flow "
                            "model. Please check the number of parameters.")
        else:
            self.materials = self.materials.append(material)

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
            List of floats denoting the depth of the nodes. The depth is
            defined in the same length unit as selected in ml.model function.
            The function defines the closest node to the desired depth.

        """
        for obs in observations:
            nodes = self.profile.iloc[
                (self.profile['x'] - obs).abs().argsort()[:1]]
            node = nodes.index.item()
            self.observations.append(node)

    def add_waterflow(self, model=0, maxit=20, tolth=1e-4, tolh=0.1, ha=1e-3,
                      hb=1e3, linitw=True, free_drainage=False,
                      seepage_face=False, qdrain=False, hseep=0, rtop=0,
                      rbot=0, rroot=0, qgwlf=False, gw_level=None, aqh=None,
                      bqh=None, hysteresis=0, ikappa=-1, wlayer=False):
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
            5 = dual-porosity model of Durner [1994] with 9 parameters.
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
            self.water_flow = {
                "MaxIt": maxit,  # Maximum No. of Iterations
                "TolTh": tolth,  # [-]
                "TolH": tolh,  # [L], default is 0.1 cm
                "TopInf": None,  # Depends on boundary condition
                "WLayer": wlayer,
                "KodTop": None,  # Depends on boundary condition
                "lInitW": linitw,
                "BotInf": None,  # Depends on boundary condition
                "qGWLF": qgwlf,
                "FreeD": free_drainage,
                "SeepF": seepage_face,
                "KodBot": None,  # Depends on boundary condition
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

            self.basic_info["lWat"] = True
            self.materials = self.get_empty_material_df()
        else:
            raise Warning("Water flow was already provided. Please delete "
                          "the old information first.")

    def add_atmospheric_bc(self, atmosphere, ldailyvar=False, lsinusvar=False,
                           llai=False, lbccycles=False, linterc=False,
                           rextinct=0.463, hcrits=1e30, tatm=0, prec=0,
                           rsoil=0, rroot=0, hcrita=1e5, rb=0, hb=0, ht=0,
                           ttop=0, tbot=0, ampl=0):
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
        if self.atmosphere_info is None:
            self.atmosphere_info = {
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

        self.atmosphere = DataFrame(data=data, index=atmosphere.index)
        self.atmosphere.update(atmosphere)

        # Enable atmosphere module
        self.basic_info["AtmInf"] = True

    def add_root_uptake(self, model=0, crootmax=0, omegac=0.5, p0=-10,
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
            if len(poptm) != self.n_materials:
                raise Warning("Length of pressure heads poptm does not "
                              "equal the number of materials!")

        if self.root_uptake is None:
            self.root_uptake = {
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
            self.basic_info["lSink"] = True
        else:
            msg = "Root water uptake model is already present in the model." \
                  " Remove the old root water uptake model first using " \
                  "ml.del_root_water_uptake()"
            raise InterruptedError(msg)

    def add_root_growth(self, irootin=0, ngrowth=None, tgrowth=None,
                        rootdepth=None, irfak=None, trmin=None, trmed=None,
                        trmax=None, xrmin=None, xrmed=None, xrmax=None,
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
        trmax: float, optional
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
                "tRMax": trmax,
                "xRMin": xrmin,
                "xRMed": xrmed,
                "xRMax": xrmax,
                "tRPeriod": trperiod
            }
            if irfak is 1:
                root_growth["tRMed"] = 0
                root_growth["xRMed"] = 0
        else:
            raise Warning("Option %s for irootin is not support in Hydrus."
                          % irootin)

        if self.root_growth is None:
            self.root_growth = root_growth
            self.basic_info["lRoot"] = True
        else:
            raise Warning("Root growth model already exists. Please delete "
                          "the old root growth model first using "
                          "ml.del_root_growth().")

    def add_solute_transport(self, model=0, epsi=0.5, lupw=False, lartd=False,
                             ltdep=False, ctola=0.0, ctolr=0.0, maxitc=20,
                             pecr=0.0, ltort=True, ibacter=0, lfiltr=False,
                             lwatdep=False, ktopch=None, kbotch=None,
                             dsurf=None, catm=None, tpulse=1):
        """Method to add solute transport to the model.

        Parameters
        ----------
        model: int, optional
            Code describing type of nonequilibrium for solute transport:
            0 = equilibrium solute transport (Default)
            1 = one-site sorption model (chemical nonequilibrium)
            2 = two-site sorption model (chemical nonequilibrium)
            3 = two kinetic sorption sites model (attachment/detachment;
            chemical nonequilibrium). This model is often used for particle
            (viruses, colloids, bacteria) transport.
            4 = two kinetic sorption sites model (attachment/detachment) (
            chemical nonequilibrium). Attachment coefficients are calculated
            using filtration theory.
            5 = dual-porosity model (mobile-immobile regions; physical
            nonequilibrium).
            6 = dual-porosity model (mobile-immobile regions) with two-site
            sorption in the mobile zone (physical and chemical nonequilibrium).
            7 = dual-permeability model (physical nonequilibrium).
            8 = dual-permeability model with either an immobile region in the
            matrix domain (physical nonequilibrium) or with two-site
            sorption in both domains (physical and chemical nonequilibrium).
        epsi: float, optional
            Temporal weighing coefficient. 0.0 for an explicit scheme (
            Default). 0.5 for a Crank-Nicholson implicit scheme. =1.0 for
            a fully implicit scheme.
        lupw: bool, optional
            True if upstream weighing formulation is to be used. False if
            the original Galerkin formulation is to be used.
        lartd: bool, optional
            True if artificial dispersion is to be added in order to fulfill
            the stability criterion PeCr (see Section 8.4.4), else False.
        ltdep: bool, optional
            True if at least one transport or reaction coefficient (ChPar) is
            temperature dependent, else False. If ltdep=True, then all
            values of ChPar(i,M) should be specified at a reference
            temperature Tr=20 degrees celsius.
        ctola: float, optional
            Absolute concentration tolerance [ML-3], the value is dependent
            on the units used (set equal to zero if nonlinear adsorption is
            not considered).
        ctolr: float, optional
            Relative concentration tolerance [-] (set equal to zero if
            nonlinear adsorption is not considered).
        maxitc: int, optional
            Maximum number of iterations allowed during any time step for
            solute transport - usually 20 (set equal to zero if nonlinear
            adsorption is not considered).
        pecr: float optional
            Stability criteria (see Section 8.4.4). Set equal to zero when
            lUpW is equal to True.
        ltort: bool, optional
            True if the tortuosity factor [Millington and Quirk, 1961] is to be
            used. False if the tortuosity factor is assumed to be equal to one.
        ibacter: int, optional
            Set equal to 1 if attachment/detachment approach is to be used
            to calculate nonequilibrium transport of viruses, colloids,
            or bacteria. Set equal to 0 if original formulations, i.e.,
            physical nonequilibrium or two-site sorption is to be used to
            describe nonequilibrium solute transport.
        lfiltr: bool, optional
            Set this logical variable equal to True. if the attachment
            coefficient is to be evaluated using the filtration theory.
        lwatdep: bool, optional
            True if at least one degradation coefficient (ChPar) is water
            content dependent.
        ktopch: int, optional
            Code which specifies the type of upper boundary condition
            1 = Dirichlet boundary condition,
            -1 = Cauchy boundary condition.
            -2 = a special type of boundary condition for volatile solutes
            as described by equation (3.46).
        kbotch: int, optional
            Code which specifies the type of lower boundary condition:
            1 = Dirichlet boundary condition,
            0 = continuous concentration profile,
            -1 = Cauchy boundary condition.
        dsurf: float, optional
            Thickness of the stagnant boundary layer, d [L]. Only when
            kTopCh=-2.
        catm: float, optional
            Concentration above the stagnant boundary layer, g_atm [ML-3].
            Only when kTopCh=-2.
        tpulse: float, optional
            Time duration of the concentration pulse [T].

        """
        if self.solute_transport is None:
            self.solute_transport = {
                "Epsi": epsi,
                "lUpW": lupw,
                "lArtD": lartd,
                "ltDep": ltdep,
                "cTolA": ctola,
                "cTolR": ctolr,
                "MaxItC": maxitc,
                "PeCr": pecr,
                "lTort": ltort,
                "iBacter": ibacter,
                "lFiltr": lfiltr,
                "iNonEqual": model,
                "lWatDep": lwatdep,
                "lDualEq": True if model == 6 else False,
                "kTopCh": ktopch,
                "kBotCh": kbotch,
                "dSurf": dsurf,
                "cAtm": catm,
                "tPulse": tpulse
            }
            self.basic_info["lChem"] = True
        else:
            raise Warning("Solute transport model already exists. Please "
                          "delete the old solute transport model first using "
                          "ml.del_solute_transport().")

    def add_solutes(self, data, ):
        if self.solutes is None:
            self.solutes = data

    def add_heat_transport(self, ampl, top_bc, top_temp, bot_bc, bot_temp,
                           tperiod=1, icampbell=1, snowmelt=0.40):
        """Method to add heat transport to the model.

        Parameters
        ----------
        ampl: float, optional
            Temperature amplitude at the soil surface [K].
        top_bc: int, optional
            Code which specifies the type of upper boundary condition:
            1 = Dirichlet boundary condition,
            -1 = Cauchy boundary condition.
        top_temp: float, optional
            Temperature of the upper boundary, or temperature of the
            incoming fluid [degree Celsius].
        bot_bc: int, optional
            Code which specifies the type of lower boundary condition:
            1 = Dirichlet boundary condition,
            0 = continuous temperature profile, zero gradient,
            -1 = Cauchy boundary condition.
        bot_temp: float, optional
            Temperature of lower boundary, or temperature of the incoming
            fluid [degree Celsius].
        tperiod: float, optional
            Time interval for completion of one temperature cycle (usually 1
            day, default) [T].
        icampbell: int, optional
            Set equal to 1 if Campbell [1985] formula is to be used to
            calculate the thermal conductivity (Default_. Set equal to 0,
            when Chung and Horton [1987] formula is to be used.
        snowmelt: float, optional
            Amount of snow that will melt during one day for each degree
            Celsius (e.g., 0.40cm default).

        """
        if self.heat_transport is None:
            self.heat_transport = {
                "Ampl": ampl,
                "tPeriod": tperiod,
                "iCampbell": icampbell,
                "SnowMF": snowmelt,
                "kTopT": top_bc,
                "tTop": top_temp,
                "kBotT": bot_bc,
                "tBot": bot_temp
            }
            self.basic_info["lTemp"] = True
        else:
            raise Warning("Heat transport model already exists. Please "
                          "delete the old heat transport model first using "
                          "ml.del_solute_transport().")

    def add_time_info(self, tinit=0, tmax=1, dt=0.1,
                      dtmin=0.0001, dtmax=0.5, print_times=False,
                      printinit=None, printmax=None, dtprint=None,
                      nsteps=None, from_atmo=False):
        """Method to produce time information.

        Parameters
        ----------
        tinit: int, optional
            Initial time of the simulation [T].
        tmax: int, optional
            Final time of the simulation [T].
        print_times: boolean, optional
            Set to True. if information of pressure head, water contents,
            temperatures, and concentrations in observation nodes, and 
            the water and solute fluxes is to be printed at a constant 
            time interval of 1 time unit.
        printinit: int, optional
            First specified print-time [T].
        printmax:int, optional
            Last specified print-time [T].
        dtprint: int, optional
            Specified time increment for print times [T].
        nsteps: str, optional
            Number of required time steps between the first specified 
            print-time (printinit) and the final specified 
            print-time (printmax)".
        from_atmo: boolean, optional.
            Set to True. If time information is determined based oon the 
            atmospheric boundary condition input.
        dt: int, optional
            Initial time increment [T].
        dtmin: int, optional 
            Minimum permitted time increment [T].        
        dtmax: int, optional
            Maximum permittedtime increment [T].
        """
        self.time_info["tInit"] = tinit
        self.time_info["tMax"] = tmax
        self.time_info["dt"] = dt
        self.time_info["dtMax"] = dtmax
        self.time_info["dtMin"] = dtmin
        if from_atmo:
            if self.atmosphere is None:
                raise Warning("Atmospheric information not provided. Please "
                              "provide atmosheric information through: "
                              "ml.add_atmospheric_bc().")
            if isinstance(self.atmosphere.index, pd.DatetimeIndex):
                times = self.atmosphere.index.dayofyear
            else:
                times = self.atmosphere.index
            self.time_info["tInit"] = times[0]
            self.time_info["tMax"] = times[-1]
            self.times = times[1:-1]
        else:
            if print_times:
                if printinit is None:
                    printinit = tinit
                if printmax is None:
                    printmax = tmax
                if nsteps is None:
                    times = np.arange(printinit, printmax,
                                      step=dtprint)
                else:
                    times = np.linspace(printinit, printmax,
                                        num=nsteps + 1)
                if printinit == tinit:
                    self.times = times[1:]
                else:
                    self.times = times
                self.time_info["MPL"] = len(self.times)
            else:
                self.time_info["MPL"] = 0
        return self.times

    def simulate(self):
        """Method to call the Hydrus-1D executable.

        """
        # Run Hydrus executable.
        cmd = [self.exe_name, self.ws_name, "-1"]
        result = run(cmd)

        return result

    def write_input(self):
        """Method to write the input files for the HYDRUS-1D simulation."""
        # 1. Write SELECTOR.IN
        self.write_selector()
        # 2. Write PROFILE.DAT
        self.write_profile()

        # 3. Write ATMOSPH.IN
        if self.basic_info["AtmInf"]:
            self.write_atmosphere()
        # 4. Write METEO.IN
        if self.basic_info["lMeteo"]:
            self.write_meteo()
        # 5. Write FIT.IN
        if self.basic_info["lInverse"]:
            self.write_fit()

    def write_selector(self, fname="SELECTOR.IN"):
        """Write the selector.in file.

        """
        self._set_bc_settings()

        # Create Header string
        string = "*** BLOCK {:{}{}{}}\n"

        # Write block A: BASIC INFORMATION
        lines = [
            "Pcp_File_Version={}\n"
            "{}{}\n"
            "{}\n"
            "LUnit TUnit MUnit\n"
            "{}\n"
            "{}\n"
            "{}\n".format(self.basic_info["iVer"],
                          string.format("A: BASIC INFORMATION ", "*", "<", 72),
                          self.basic_info["Hed"], self.description,
                          self.basic_info["LUnit"], self.basic_info["TUnit"],
                          self.basic_info["MUnit"])
        ]

        vars_list = [["lWat", "lChem", "lTemp", "lSink", "lRoot", "lShort",
                      "lWDep", "lScreen", "AtmInf", "lEquil", "lInverse",
                      "\n"],
                     ["lSnow", "lHP1", "lMeteo", "lVapor", "lActRSU", "lFlux",
                      "lIrrig", "\n"]]

        for variables in vars_list:
            lines.append("  ".join(variables))
            lines.append("  ".join("t" if self.basic_info[var] else "f" for
                                   var in variables[:-1]))
            lines.append("\n")

        lines.append("NMat NLay CosAlfa \n"
                     "{} {} {}\n".format(self.n_materials, self.n_layers,
                                         self.basic_info["CosAlfa"]))

        # Write block B: WATER FLOW INFORMATION
        lines.append(string.format("B: WATER FLOW INFORMATION ", "*", "<", 72))
        lines.append("MaxIt  TolTh  TolH   (maximum number of iterations and "
                     "tolerances)\n")
        variables = ["MaxIt", "TolTh", "TolH"]
        lines.append(
            "   ".join([str(self.water_flow[var]) for var in variables]))
        lines.append("\n")

        vars_list = [["TopInf", "WLayer", "KodTop", "lInitW", "\n"],
                     ["BotInf", "qGWLF", "FreeD", "SeepF", "KodBot", "qDrain",
                      "hSeep", "\n"]]

        upper_condition = (self.water_flow["KodTop"] < 0
                           and not self.water_flow["TopInf"])

        lower_condition = ((self.water_flow["KodBot"] < 0)
                           and not self.water_flow["BotInf"]
                           and not self.water_flow["qGWLF"]
                           and not self.water_flow["FreeD"]
                           and not self.water_flow["SeepF"]
                           )

        if upper_condition or lower_condition:
            vars_list.append(["rTop", "rBot", "rRoot", "\n"])

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
        lines.append(string.format("C: TIME INFORMATION ", "*", "<", 72))
        vars_list = [
            ["dt", "dtMin", "dtMax", "dMul", "dMul2", "ItMin", "ItMax",
             "MPL", "\n"], ["tInit", "tMax", "\n"],
            ["lPrint", "nPrintSteps", "tPrintInterval", "lEnter", "\n"]]
        for variables in vars_list:
            lines.append(" ".join(variables))
            values = []
            for var in variables[:-1]:
                val = self.time_info[var]
                if val is True:
                    values.append("t")
                elif val is False:
                    values.append("f")
                else:
                    values.append(str(val))
            values.append("\n")
            lines.append(" ".join(values))

        lines.append("TPrint(1),TPrint(2),...,TPrint(MPL)\n")
        for i in range(int(len(self.times) / 6) + 1):
            lines.append(
                " ".join([str(time) for time in self.times[i * 6:i * 6 + 6]]))
            lines.append("\n")

        # Write BLOCK D: Root Growth Information
        if self.basic_info["lRoot"]:
            lines.append(
                string.format("D: ROOT GROWTH INFORMATION ", "*", "<", 72))
            lines.append("iRootDepthEntry\n{}\n".format(self.root_growth[
                                                            "iRootIn"]))
            d = self.root_growth.copy()
            d.pop("iRootIn")
            d["\n"] = "\n"
            lines.append("    ".join(d.keys()))
            lines.append("    ".join(str(p) for p in d.values()))

        # Write Block E - Heat transport information
        if self.basic_info["lTemp"]:
            lines.append(string.format("E: HEAT TRANSPORT INFORMATION ",
                                       "*", "<", 72))
            lines.append(self.heat_parameters.to_string(index=False))
            lines.append(
                "tAmpl tPeriod Campbell SnowMF lDummy lDummy lDummy "
                "lDummy lDummy\n"
                "{} {} {} {} f f f f f\n"
                "kTopT TTop kBotT TBot\n"
                "{} {} {} {}\n".format(self.heat_transport["Ampl"],
                                       self.heat_transport["tPeriod"],
                                       self.heat_transport["iCampbell"],
                                       self.heat_transport["SnowMF"],
                                       self.heat_transport["kTopT"],
                                       self.heat_transport["tTop"],
                                       self.heat_transport["kBotT"],
                                       self.heat_transport["tBot"]))

        # Write Block F - Solute transport information
        if self.basic_info["lChem"]:
            raise NotImplementedError

        # Write Block G - Root water uptake information
        if self.basic_info["lSink"]:
            lines.append(string.format("G: ROOT WATER UPTAKE INFORMATION ",
                                       "*", "<", 72))
            vars_list = [["iMoSink", "cRootMax", "OmegaC", "\n"]]

            if self.root_uptake["iMoSink"] is 0:
                vars_list.append(
                    ["P0", "P2H", "P2L", "P3", "r2H", "r2L", "\n"])

            for variables in vars_list:
                lines.append(" ".join(variables))
                lines.append("    ".join(str(self.root_uptake[var]) for var in
                                         variables[:-1]))
                lines.append("\n")

            lines.append("POptm(1),POptm(2),...,POptm(NMat)\n")
            lines.append("    ".join(str(p) for p in self.root_uptake[
                "POptm"]))
            lines.append("\n")

        # Write Block J - Inverse solution information
        if self.basic_info["lInverse"]:
            raise NotImplementedError("The inverse modeling module from "
                                      "Hydrus-1D will not be supported. "
                                      "Python packages are used for this.")

        # Write Block K – Carbon dioxide transport information

        # Write Block L – Major ion chemistry information
        if self.basic_info["lChem"]:
            raise NotImplementedError

        # Write Block M – Meteorological information
        if self.basic_info["lMeteo"]:
            raise NotImplementedError

        # Write END statement
        lines.append(string.format("END OF INPUT FILE SELECTOR.IN ",
                                   "*", "<", 72))

        # Write the actual file
        fname = os.path.join(self.ws_name, fname)
        with open(fname, "w") as file:
            file.writelines(lines)
        print("Successfully wrote {}".format(fname))

    def write_atmosphere(self, fname="ATMOSPH.IN"):
        """Method to write the ATMOSPH.IN file

        """
        # 1 Write Header information
        lines = ["Pcp_File_Version={}\n".format(
            self.basic_info["iVer"]),
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
            val = self.atmosphere_info[var]
            if var:
                if val is True:
                    vals.append("t")
                elif val is False:
                    vals.append("f")
                else:
                    vals.append(str(val))
        lines.append(" ".join(vals))

        lines.append("\nhCritS (max. allowed pressure head at the soil "
                     "surface)\n{}\n".format(self.atmosphere_info["hCritS"]))

        lines.append(self.atmosphere.to_string(index=False))
        lines.append("\nend*** END OF INPUT FILE ATMOSPH.IN "
                     "**********************************\n")
        # Write the actual file
        fname = os.path.join(self.ws_name, fname)
        with open(fname, "w") as file:
            file.writelines(lines)
        print("Successfully wrote {}".format(fname))

    def write_profile(self, fname="PROFILE.DAT"):
        """Method to write the profile.dat file.

        """
        # 1 Write Header information
        lines = ["Pcp_File_Version={}\n"
                 "0\n"
                 "{} {} {} {}".format(self.basic_info["iVer"],
                                      self.profile.index.size,
                                      self.n_solutes, 0, 0),
                 # 2. Write the profile data
                 self.profile.to_string(),
                 # 3. Write observation points
                 "\n{}\n".format(len(self.observations)),
                 "".join(["   {}".format(i) for i in self.observations])]

        # Write the actual file
        fname = os.path.join(self.ws_name, fname)
        with open(fname, "w") as file:
            file.writelines(lines)
        print("Successfully wrote {}".format(fname))

    def write_fit(self, fname="FIT.IN"):
        raise NotImplementedError

    def write_meteo(self, fname="METEO.IN"):
        raise NotImplementedError

    def read_output(self):
        raise NotImplementedError

    def read_profile(self, fname="PROFILE.OUT"):
        path = os.path.join(self.ws_name, fname)
        data = read_profile(path=path)
        return data

    def read_nod_inf(self, fname="NOD_INF.OUT", times=None):
        path = os.path.join(self.ws_name, fname)
        data = read_nod_inf(path=path, times=times)
        return data

    def read_run_inf(self, fname="RUN_INF.OUT", usecols=None):
        path = os.path.join(self.ws_name, fname)

        if usecols is None:
            usecols = ["TLevel", "Time", "dt", "Iter", "ItCum", "KodT",
                       "KodB", "Convergency", ]
            if self.solute_transport is not None:
                usecols.append("IterC")

        data = read_run_inf(path, usecols=usecols)
        return data

    def read_balance(self, fname="BALANCE.OUT", usecols=None):
        path = os.path.join(self.ws_name, fname)
        if not os.path.exists(path):
            raise FileNotFoundError(
                "File {} has not been found.".format(path))

        if usecols is None:
            usecols = ["Area", "W-volume", "In-flow", "h Mean", "Top Flux",
                       "Bot Flux", "WatBalT", "WatBalR"]
            if self.solute_transport is not None:
                usecols.append("ConcVol", "ConcVolIm", "cMean", "CncBalT",
                               "CncBalR")

            if self.heat_transport is not None:
                usecols.append("TVol", "TMean")

            if self.CO2Transport is not None:
                usecols.append("COVol", "COMean", "CO2BalT", "CncBalT")

            if self.water_flow["iModel"] in [5, 6, 7]:
                usecols.append("W-VolumeI", "cMeanIm")

        data = read_balance(path=path, usecols=usecols)

        return data

    def read_obs_node(self, fname="OBS_NODE.OUT", nodes=None, times=None):
        path = os.path.join(self.ws_name, fname)
        if nodes is None:
            nodes = self.observations

        data = read_obs_node(path=path, nodes=nodes, times=times)
        return data

    def read_i_check(self, fname="I_CHECK.OUT", times=None):
        raise NotImplementedError

    def read_tlevel(self, fname="T_LEVEL.OUT", usecols=None):
        path = os.path.join(self.ws_name, fname)

        if usecols is None:
            usecols = ["Time", "rTop", "rRoot", "vTop", "vRoot",
                       "vBot", "sum(rTop)", "sum(rRoot)", "sum(vTop)",
                       "sum(vRoot)", "sum(vBot)", "hTop", "hRoot", "hBot",
                       "RunOff", "Volume", ]

            if self.water_flow["iModel"] > 4:
                usecols.append("Cum(WTrans)")
            if self.basic_info["lSnow"]:
                usecols.append("SnowLayer")

        data = read_tlevel(path=path, usecols=usecols)

        return data

    def read_alevel(self, fname="A_LEVEL.OUT", usecols=None):
        path = os.path.join(self.ws_name, fname)
        data = read_alevel(path=path, usecols=usecols)
        return data

    def read_solutes(self, fname="SOLUTE.OUT"):
        path = os.path.join(self.ws_name, fname)
        data = read_solute(path=path)
        return data

    def get_empty_material_df(self):
        """Returns an empty dataframe with the soil parameters as columns.

        return
        ----------
        pandas.DataFrame
            Pandas DataFrame with the soil parameters as columns.

        """
        models = {
            0: ["thr", "ths", "Alfa", "n", "Ks", "l"],
            1: ["thr", "ths", "Alfa", "n", "Ks", "l", "thm", "tha", "thk",
                "Kk"],
            2: ["thr", "ths", "Alfa", "n", "Ks", "l"],
            3: ["thr", "ths", "Alfa", "n", "Ks", "l"],
            4: ["thr", "ths", "Alfa", "n", "Ks", "l"],
            5: ["thr", "ths", "Alfa", "n", "Ks", "l", "w", "Alfa2", "n2"],
            6: ["thr", "ths", "Alfa", "n", "Ks", "l", "thr_im", "ths_im",
                "omega"],
            7: ["thr", "ths", "Alfa", "n", "Ks", "l", "thr_im", "ths_im",
                "Alfa_im", "n_im", "Ka"],
            9: list(range(17))
        }

        cols = models[self.water_flow["iModel"]]

        if self.solute_transport is not None:
            models = {
                1: [],
                2: [],
                3: [],
                4: [],
                5: [],
                6: [],
                7: [],
                8: []
            }
            cols.extend(models[self.solute_transport["lNonEqual"]])

        return DataFrame(columns=cols)

    def get_empty_heat_df(self):
        """Returns an empty DataFrame to fill in the heat parameters.

        """
        columns = ["thn", "tho", "lambda", "b1", "b2", "b3", "Cn", "C0", "Cw"]
        return DataFrame(columns=columns, index=self.materials.index)

    def get_empty_solute_df(self):
        """Returns an empty DataFrame with the solute parameters as columns.
        """
        models = {
            1: "",

        }

        df = DataFrame(columns=models[self.solute_transport["lNonEqual"]])
        return df

    def _set_bc_settings(self):
        """Internal method to set the boundary condition settings.

        Returns
        -------

        Notes
        -----
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

        """
        kodtop = -1
        topinf = False
        kodbot = 1
        botinf = False

        # In case of a seepage face or free drainage BC set KodBot=-1.
        if self.water_flow["SeepF"] or self.water_flow["FreeD"]:
            kodbot = -1

        if self.basic_info["AtmInf"]:
            topinf = True
            kodtop = -1

        self.water_flow["KodTop"] = kodtop
        self.water_flow["TopInf"] = topinf
        self.water_flow["KodBot"] = kodbot
        self.water_flow["BotInf"] = botinf


# Copy all the docstrings from the read methods
Model.read_profile.__doc__ = "{}".format(read_profile.__doc__)
Model.read_nod_inf.__doc__ = "{}".format(read_nod_inf.__doc__)
Model.read_run_inf.__doc__ = "{}".format(read_run_inf.__doc__)
Model.read_balance.__doc__ = "{}".format(read_balance.__doc__)
Model.read_obs_node.__doc__ = "{}".format(read_obs_node.__doc__)
Model.read_i_check.__doc__ = "{}".format(read_i_check.__doc__)
Model.read_tlevel.__doc__ = "{}".format(read_tlevel.__doc__)
Model.read_alevel.__doc__ = "{}".format(read_alevel.__doc__)
Model.read_solutes.__doc__ = "{}".format(read_solute.__doc__)
