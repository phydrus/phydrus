"""
Example script of how to setup a basic HYDRUS-1D model using Pydrus.
"""

import os
import pydrus as ps
import pandas as pd

ws = "DRAINAGE"
exe = os.path.join(os.getcwd(), "hydrus")

# Create the basic model
desc = "Infiltration and drainage in a large caisson"

ml = ps.Model(exe_name=exe, ws_name=ws, name="model", description=desc,
              mass_units="mmol", time_unit="min", length_unit="cm")
ml.time_information["tMax"] = 500
ml.time_information["dtMax"] = 500
ml.time_information["TPrint(MPL)"] = 500

# Create a Profile
# data = ps.Profile.read_profile(ws=ws)
# p = ps.Profile(data=data)
# p.observations = [1,2,3]
#
# ml.add_profile(p)
#
# # p.write_file()


rs = ml.simulate()
