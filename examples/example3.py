"""
Example script of how to setup a basic HYDRUS-1D model using Pydrus.
"""

import os
import pydrus as ps
import pandas as pd

ws = "example3"
exe = os.path.join(os.getcwd(), "hydrus")

# Create the basic model
desc = "Root uptake with meteorological data"

ml = ps.Model(exe_name=exe, ws_name=ws, name="model", description=desc,
              mass_units="mmol", time_unit="days", length_unit="cm")

# Time info
ml.time_information["tInit"] = 0
ml.time_information["tMax"] = 214
ml.time_information["dtMax"] = 1
ml.time_information["MPL"] = 1

# Add materials
m = pd.DataFrame(columns=["thr", "ths", "Alfa", "n", "Ks", "l"], index=[1],
                 data=[[0.095, 0.41, 0.019, 1.31, 3.4, 0.5]])
ml.add_material(m)

# Profile
profile = ps.create_profile(bot=-100, dx=1, h=-70.5614)
ml.add_profile(profile)

# Water flow info
ml.add_waterflow(maxit=20, tolh=1, linitw=False, kodbot=-1,
                 seepage_face=True, hseep=-60)

# Atmospheric data
atm = pd.read_csv("data/ex3.csv")
ml.add_atmosphere(atm)

# Root uptake: does not work yet!
# ml.add_rootwater_uptake(model=0, poptm=[-25] * m.index.size)

ml.write_files()
ml.simulate()

df = ml.read_tlevel()
df['vBot[L/T]'].plot()
