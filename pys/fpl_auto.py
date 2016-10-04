import os

import fpl
import fpl_calc
from fpl_lmp_large import job as lmp_large_job
from fpl_lmp_small import job as lmp_small_job
from fpl_orca import job as orca_job

def get_enthalpy_solvation(solute, solvent):
	########################################
	#solute = "pb2+"
	#solvent = "THTO"
	run_name = "%s_%s" % (solvent, solute)

	# Generate initial object
	fpl_obj = fpl.fpl_job(run_name, solvent, solute)
	fpl_obj.cml_dir="/fs/home/hch54/frazier-pipeline/cml/"

	# Set parameters
	fpl_obj.num_solvents=1

	# Generate system
	fpl_obj.generate_system()

	########################################
	# Add the tasks here

	## Task 1 - Large Lammps Simulation
	### PARAMETERS
	task1 = run_name + "_large_lammps"
	fpl_obj.queue=None
	fpl_obj.procs=1
	fpl_obj.lmp_run_len=10000
	fpl_obj.trj_file=None
	### ADD TASK
	task = lmp_large_job(fpl_obj, task1)
	fpl_obj.add_task(task)

	## Task 2 - Small Lammps Simulation
	### PARAMETERS
	task2 = run_name + "_small_lammps"
	fpl_obj.queue=None
	fpl_obj.procs=1
	fpl_obj.lmp_run_len=10000
	### ADD TASK
		# Note, you can always overwrite the callback function
		# tsk = lmp_small_job(fpl_obj, task2)
		# tsk.callback = None
		# fpl_obj.add_task( task2, tsk )
	task = lmp_small_job(fpl_obj, task2)
	fpl_obj.add_task(task)

	## Task 3 - Orca Simulation
	### PARAMETERS
	task3 = run_name + "_orca"
	fpl_obj.queue="batch"
	fpl_obj.procs=4
	fpl_obj.route = "! OPT B97-D3 SV GCP(DFT/TZ) ECP{def2-TZVP} Grid7 SlowConv LooseOpt"
	### ADD TASK
	task = orca_job(fpl_obj, task3)
	fpl_obj.add_task(task)

	########################################

	# Run the simulation here
	fpl_obj.start(save=False)

	########################################

	## Task 4 - Calculate Enthalpy of Solvation
	### PARAMETERS
	task4 = run_name + "_Hsolv"
	fpl_obj.queue = "batch"
	fpl_obj.procs = 4

	fpl_obj.route = "! B97-D3 SV GCP(DFT/TZ) ECP{def2-TZVP} Grid7 SlowConv"
	fpl_obj.extra_section = "%basis aux auto NewECP Pb \"def2-SD\" \"def2-TZVP\" end NewECP Cs \"def2-SD\" \"def2-TZVP\" end NewGTO S \"def2-TZVP\" end end" 
	fpl_obj.charge_and_multiplicity = "0 1"

	fpl_obj.route_solute = "! B97-D3 SV GCP(DFT/TZ) ECP{def2-TZVP} Grid7 SlowConv"
	fpl_obj.extra_section_solute = "%basis aux auto NewECP Pb \"def2-SD\" \"def2-TZVP\" end NewECP Cs \"def2-SD\" \"def2-TZVP\" end end" 
	fpl_obj.charge_and_multiplicity_solute = "0 1"

	fpl_obj.route_solvent = "! B97-D3 SV GCP(DFT/TZ) ECP{def2-TZVP} Grid7 SlowConv"
	fpl_obj.extra_section_solvent = "%basis aux auto NewGTO S \"def2-TZVP\" end end" 
	fpl_obj.charge_and_multiplicity_solvent = "0 1"

	### ADD TASK
	tasks = fpl_calc.enthalpy_solvation(fpl_obj, task4)
	fpl_obj.add_task(tasks, parallel=True)

	fpl_obj.start(save=False)

	H_solv = fpl_calc.post_enthalpy_solvation(fpl_obj)

	return H_solv