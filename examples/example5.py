"""
Example script of how to setup a basic HYDRUS-1D model using Pydrus.

Author: R.A. Collenteur, University of Graz, 2019

"""

import os

import pandas as pd

import pydrus as ps

ws = "example5"
exe = os.path.join(os.getcwd(), "hydrus")

# Create the basic model
desc = "Example 5 - Grass Field Problem (Hupselse Beek 1982)"

ml = ps.Model(exe_name=exe, ws_name=ws, name="model", description=desc,
              mass_units="-", time_unit="days", length_unit="m",
              print_screen=True)

times = ml.add_time_info(tinit=90, tmax=273, print_times=True)

# Water inflow parameters
ml.add_waterflow(linitw=False, top_bc=3, bot_bc=4, ha=1e-6, hb=1e4)
ml.add_solute_transport(tpulse=1)

m = ml.get_empty_material_df(n=2)
m.loc[1:2] = [[0.01, 0.399, 0.0174, 1.3757, 209.75, 0.5, 1.9, 130, 1, 0],
              [0.01, 0.339, 0.0139, 1.6024, 45.34, 0.5, 1.9, 100, 1, 0]]
ml.add_material(m)

profile = ps.create_profile(0, [-50, -100], h=-200, dx=10, mat=m.index,
                            conc=0, sconc=0)
ml.add_profile(profile)
ml.add_obs_nodes([10, 20])

atm = pd.read_csv("data/ex2.csv", index_col=0)
atm["cTop"] = 0
atm["cBot"] = 0
atm.loc[0, "cTop"] = 1
ml.add_atmospheric_bc(atm)

sol1 = ml.get_empty_solute_df()
sol1["beta"] = 1
ml.add_solute(sol1, difw=0, difg=0)

ml.write_input()
rs = ml.simulate()

df = ml.read_solutes()
df.loc[:, ["cBot"]].plot(subplots=True)
