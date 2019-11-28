"""
This example implements a basic example of water flow in a layered profile,
following the example from the following link:

https://www.pc-progress.com/Downloads/Tutorials/Tutorial_H1D_4.pdf

A steady state simulation is performed and the pressure heads are used as the
initial condition for a transient model.

Author: R.A. Collenteur, University of Graz, 2019

"""

import os

import pandas as pd
import pydrus as ps

ws = "example4"
exe = os.path.join(os.getcwd(), "hydrus")

# Create the basic model
desc = "Steady state water flow in a layered soil profile"

ml = ps.Model(exe_name=exe, ws_name=ws, name="model", description=desc,
              mass_units="-", time_unit="days", length_unit="cm")

# Time info
times = ml.add_time_info(tmax=100, print_times=True, dt=0.001,
                         dtmin=0.000001)
# Water flow info
ml.add_waterflow(linitw=False, top_bc=1, bot_bc=4, ha=1e-6, hb=1e4, rtop=-0.12)

# Add materials
m = pd.DataFrame(columns=["thr", "ths", "Alfa", "n", "Ks", "l"],
                 data=[[0.065, 0.42, 0.016, 1.94, 95.040, 0.5],
                       [0.035, 0.42, 0.015, 3.21, 311.04, 0.5],
                       [0.042, 0.40, 0.016, 1.52, 38.880, 0.5],
                       [0.044, 0.40, 0.028, 2.01, 864.00, 0.5],
                       [0.039, 0.40, 0.023, 2.99, 1209.6, 0.5],
                       [0.030, 0.42, 0.021, 2.99, 1209.6, 0.5],
                       [0.021, 0.39, 0.021, 2.99, 1209.6, 0.5]],
                 index=range(1, 8))
ml.add_material(m)

# Profile
profile = ps.create_profile(bot=[-7, -19, -24, -28, -50, -75, -100], dx=1,
                            h=-100, mat=m.index)
ml.add_profile(profile)
ml.add_obs_nodes([50, 100])

# run steady state simulation
ml.write_input()
ml.simulate()

# Update the initial head for the transient simulation with the steady state
df = ml.read_nod_inf()
ml.profile["h"] = df["Head"].values

# Atmospheric data
atm = pd.read_csv("data/ex4.csv", decimal=",", sep=";")
ml.add_atmospheric_bc(atm)

times1 = ml.add_time_info(tmax=360, print_times=True, dt=0.001,
                          dtmin=0.000001)

# Root uptake
ml.add_root_uptake(model=0, poptm=[-125] * m.index.size)

ml.write_input()
ml.simulate()

df = ml.read_tlevel()
df['vBot'].plot()
