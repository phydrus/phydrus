"""
Example script of how to setup a basic HYDRUS-1D model using Pydrus.

Author: R.A. Collenteur, University of Graz, 2019

"""

import os

import pydrus as ps

ws = "example"
exe = os.path.join(os.getcwd(), "hydrus")

# Create the basic model
desc = "Infiltration and drainage in a large caisson"

ml = ps.Model(exe_name=exe, ws_name=ws, name="model", description=desc,
              mass_units="mmol", time_unit="min", length_unit="cm")

times = ml.add_time_info(tinit=90, tmax=273, print_times=True, dt=0.1,
                         dtmax=0.5, printinit=120)

ml.add_waterflow(top_bc=1, bot_bc=0, rtop=0)

m = ml.get_empty_material_df(n=2)
m.loc[1:3] = [[0.08, 0.3421, 0.03, 5, 1, -0.5],
              [0.08, 0.3421, 0.03, 5, 0.1, -0.5]]

ml.add_material(m)

profile = ps.create_profile(h=0.342)
profile.loc[5:11, "Mat"] = 2
ml.add_profile(profile)

# atm = pd.read_csv("seep/orig/ATMOSPH.IN", skiprows=5, skipfooter=1,
#                   skipinitialspace=True, delim_whitespace=True)
# ml.add_atmosphere(atm)

ml.write_input()
rs = ml.simulate()
ml.plots.profile()

# df = ml.read_tlevel()
# df['vBot[L/T]'].plot()
