"""
Example script of how to setup a basic HYDRUS-1D model using Pydrus.
"""

import os
import pydrus as ps
import pandas as pd

ws = "TEST"
exe = os.path.join(os.getcwd(), "hydrus")

# Create the basic model
desc = "Infiltration and drainage in a large caisson"

ml = ps.Model(exe_name=exe, ws_name=ws, name="model", description=desc,
              mass_units="mmol", time_unit="min", length_unit="cm")
ml.time_information["tMax"] = 500
ml.time_information["dtMax"] = 500
ml.time_information["TPrint(MPL)"] = 500

m = pd.DataFrame(data=[[0.08, 0.3421, 0.03, 5, 1, -0.5]],
                 columns=["thr", "ths", "Alfa", "n", "Ks", "l"])
ml.add_material(m)


profile = ps.create_profile(h=0.342)

ml.add_profile(profile)

ml.write_files()
rs = ml.simulate()
