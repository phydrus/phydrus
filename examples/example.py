"""
Example script of how to setup a basic HYDRUS-1D model using Pydrus.
"""

import os
import pydrus as ps
import pandas as pd

ws = "example"
exe = os.path.join(os.getcwd(), "hydrus")

# Create the basic model
desc = "Infiltration and drainage in a large caisson"

ml = ps.Model(exe_name=exe, ws_name=ws, name="model", description=desc,
              mass_units="mmol", time_unit="min", length_unit="cm")

ml.time_information["tInit"] = 90
ml.time_information["tMax"] = 273
ml.time_information["dt"] = 0.1
ml.time_information["dtMax"] = 0.5

ml.time_information["TPrint(1)"] = 120
ml.time_information["TPrint(MPL)"] = 273

ml.add_waterflow()

m = pd.DataFrame(data=[[0.08, 0.3421, 0.03, 5, 1, -0.5],
                       [0.08, 0.3421, 0.03, 5, 0.1, -0.5]],
                 columns=["thr", "ths", "Alfa", "n", "Ks", "l"], index=[1, 2])
ml.add_material(m)

profile = ps.create_profile(h=0.342)
profile.loc[5:11, "Mat"] = 2
ml.add_profile(profile)

# atm = pd.read_csv("seep/orig/ATMOSPH.IN", skiprows=5, skipfooter=1,
#                   skipinitialspace=True, delim_whitespace=True)
# ml.add_atmosphere(atm)

ml.write_files()
rs = ml.simulate()
ml.plots.profile()

# df = ml.read_tlevel()
# df['vBot[L/T]'].plot()
