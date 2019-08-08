"""
Example script of how to setup a basic HYDRUS-1D model using Pydrus.
"""

import os
import pydrus as ps
import pandas as pd

ws = "example2"
exe = os.path.join(os.getcwd(), "hydrus")

# Create the basic model
desc = "Example 2 - Grass Field Problem (Hupselse Beek 1982)"

ml = ps.Model(exe_name=exe, ws_name=ws, name="model", description=desc,
              mass_units="-", time_unit="days", length_unit="cm")

ml.time_information["tInit"] = 90
ml.time_information["tMax"] = 273

# Water inflow parameters
ml.add_waterflow(kodbot=-1, linitw=False, free_drainage=True, ha=1e-6, hb=1e4)

m = pd.DataFrame(columns=["thr", "ths", "Alfa", "n", "Ks", "l"],
                 data=[[0.0001, 0.399, 0.0174, 1.3757, 29.75, 0.5],
                       [0.0001, 0.339, 0.0139, 1.6024, 45.34, 0.5]],
                 index=[1, 2])

ml.add_material(m)

profile = ps.create_profile(0, -230, h=-200, dx=10)
profile.loc[11:, ["Mat", "Lay"]] = 2
ml.add_profile(profile)
ml.add_observations([10, 20])


atm = pd.read_csv("data/ex2.csv", index_col=0)
ml.add_atmosphere(atm)

ml.add_rootwater_uptake(model=0, poptm=[-25, -25])

ml.write_files()
rs = ml.simulate()

df = ml.read_tlevel()
df['vBot[L/T]'].plot()
