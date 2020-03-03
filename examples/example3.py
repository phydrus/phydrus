"""
Example script of how to setup a basic HYDRUS-1D model using Pydrus.

Author: R.A. Collenteur, University of Graz, 2019

"""

import os

import pandas as pd

import phydrus as ps

ws = "example3"
exe = os.path.join(os.getcwd(), "hydrus")

# Create the basic model
desc = "Root uptake with meteorological data"

ml = ps.Model(exe_name=exe, ws_name=ws, name="model", description=desc,
              mass_units="mmol", time_unit="days", length_unit="cm")

# Time info
times = ml.add_time_info(tmax=213, print_times=True)

# Water flow info
ml.add_waterflow(maxit=20, tolh=1, linitw=False, top_bc=3, bot_bc=6, hseep=-60)

# Add materials
m = ml.get_empty_material_df(n=1)
m.loc[[1, 1]] = [[0.095, 0.41, 0.019, 1.31, 3.4, 0.5]]
ml.add_material(m)

# Profile
profile = ps.create_profile(bot=-100, dx=10, h=-70.5614)
ml.add_profile(profile)

# Atmospheric data
atm = pd.read_csv("data/ex3.csv", index_col=0)
atm = atm.drop("RootDepth", 1)
ml.add_atmospheric_bc(atm)

# Root uptake
ml.add_root_uptake(model=0, poptm=[-25])

ml.write_input()
ml.simulate()

df = ml.read_tlevel()
df['vBot'].plot()
