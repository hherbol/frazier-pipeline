# Code for orca's dft simulation

# Python Imports
import os
import numpy as np
import cPickle as pickle

# fpl Imports
import fpl_constants, fpl_utils

# Clancelot Imports
import files, units, orca

def job(run_name, prev_run_name, system, solvent_name, path=os.getcwd()+"/", debug=False,
	    route = "! OPT B97-D3 SV GCP(DFT/TZ) ECP{def2-TZVP} Grid7 SlowConv LooseOpt",
	    extra_section = "%basis aux auto NewECP Pb \"def2-SD\" \"def2-TZVP\" end NewECP Cs \"def2-SD\" \"def2-TZVP\" end NewGTO S \"def2-TZVP\" end end",
 		implicit_solvent = True,
 		dft_params = {"queue":"long", "mem":5000, "procs":1, "charge":0, "xhost":None} ):
	
	# Step 0 - Ensure proper variables
	if not path.endswith("/"): path += "/"
	if type(route) is int:
		route = ["! OPT B97-D3 SV GCP(DFT/TZ) ECP{def2-TZVP} Grid7 SlowConv LooseOpt",
				 "! OPT B97-D3 SV GCP(DFT/TZ) ECP{def2-TZVP} Grid7 SlowConv LooseOpt",
				 "! OPT B97-D3 SV GCP(DFT/TZ) ECP{def2-TZVP} Grid7 SlowConv LooseOpt"][route]
		raise Exception("Currently all 0, 1, and 2 are the same (lvl 0, lowest)")

	# Step 1 - Generate the structure for our orca simulation
	## Grab last frame of this smaller, newly optimized system and store
	elements_by_index = fpl_utils.get_xyz_elems(prev_run_name)
	xyz = files.read_xyz(system.name+'.xyz')
	## Convert Ba to Pb (we had used Ba in LAMMPs, but in reality we want Pb)
	for a,b in zip(system.atoms,xyz[-1]):
		a.x, a.y, a.z, a.element = b.x, b.y, b.z, elements_by_index[b.element] if elements_by_index[b.element] != "Ba" else "Pb"
		if any( [np.isnan(x) for x in (a.x,a.y,a.z)] ):
			return None
	## Generate outputs
	files.write_xyz(system.atoms, run_name+'.cluster')	
	pickle.dump(system, open(run_name+'.pickle', 'w'))

	os.chdir(path)

	if implicit_solvent:
		route += " COSMO"
		extra_section = "%cosmo SMD true epsilon " + str(fpl_constants.solvent[solvent_name]["dielectric"]) + " end " + extra_section

	if not debug:
		orca.job(run_name, route, system.atoms, queue=dft_params["queue"], mem=dft_params["mem"], procs=dft_params["procs"], extra_section=extra_section, charge=dft_params["charge"], xhost=dft_params["xhost"])

	return system
