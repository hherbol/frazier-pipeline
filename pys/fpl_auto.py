import os
import types

import jobs

import fpl
import fpl_calc
from fpl_lmp_large import job as lmp_large_job
from fpl_lmp_small import job as lmp_small_job
from fpl_orca import job as orca_job

def get_enthalpy_solvation(solute, solvent, num_solvents=1, on_queue=False,
	queue="batch", nprocs="1", xhost=None, unit="kT_300",
	charge_and_multiplicity="0 1", charge_and_multiplicity_solute="0 1",
	charge_and_multiplicity_solvent="0 1", name_append="",
	route_lvls=[1,1,1,1]):
	if num_solvents < 1:
		raise Exception("You must have a number of solvents greater than or\
		equal to 1")

	if on_queue:
		pysub_str = """import fpl_auto
e_solv = fpl_auto.get_enthalpy_solvation("$SOLUTE","$SOLVENT", num_solvents=$NUM_SOLVENTS, unit="$UNIT", charge_and_multiplicity="$CHARGE_AND_MULTIPLICITY", name_append="$NAME_APPEND", route_lvls=$ROUTE_LVLS)
print e_solv
"""
		pysub_str = pysub_str.replace("$SOLUTE",solute)
		pysub_str = pysub_str.replace("$SOLVENT",solvent)
		pysub_str = pysub_str.replace("$NUM_SOLVENTS",str(num_solvents))
		pysub_str = pysub_str.replace("$UNIT",unit)
		pysub_str = pysub_str.replace("$CHARGE_AND_MULTIPLICITY",charge_and_multiplicity)
		pysub_str = pysub_str.replace("$NAME_APPEND",name_append)
		pysub_str = pysub_str.replace("$ROUTE_LVLS",str(route_lvls))

		job_name = "%s_%s%s.py" % (solute,solvent,name_append)
		fptr = open(job_name, "w")
		fptr.write(pysub_str)
		fptr.close()
		running_job = jobs.pysub(job_name, nprocs=nprocs, queue=queue, xhost=xhost, path=os.getcwd(), remove_sub_script=True)
		def enthalpy(self):
			if not self.is_finished():
				return None
			else:
				H = float(open(self.name+".log", "r").read().strip().split("\n")[-1])
				return H
		running_job.enthalpy = types.MethodType( enthalpy, running_job)

		return running_job

	########################################
	run_name = "%s_%s" % (solvent, solute)
	run_name += name_append

	# Generate initial object
	fpl_obj = fpl.fpl_job(run_name, solvent, solute)
	fpl_obj.cml_dir="/fs/home/hch54/frazier-pipeline/cml/"

	# Set parameters
	fpl_obj.num_solvents=num_solvents

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
	fpl_obj.route = route_lvls[0]
	fpl_obj.charge_and_multiplicity = charge_and_multiplicity
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

	fpl_obj.route = route_lvls[1]
	fpl_obj.extra_section = "%basis aux auto NewECP Pb \"def2-SD\" \"def2-TZVP\" end NewECP Cs \"def2-SD\" \"def2-TZVP\" end NewGTO S \"def2-TZVP\" end end" 
	fpl_obj.charge_and_multiplicity = charge_and_multiplicity

	fpl_obj.route_solute = route_lvls[2]
	fpl_obj.extra_section_solute = "%basis aux auto NewECP Pb \"def2-SD\" \"def2-TZVP\" end NewECP Cs \"def2-SD\" \"def2-TZVP\" end NewGTO S \"def2-TZVP\" end end" 
	fpl_obj.charge_and_multiplicity_solute = charge_and_multiplicity_solute

	fpl_obj.route_solvent = route_lvls[3]
	fpl_obj.extra_section_solvent = "%basis aux auto NewECP Pb \"def2-SD\" \"def2-TZVP\" end NewECP Cs \"def2-SD\" \"def2-TZVP\" end NewGTO S \"def2-TZVP\" end end" 
	fpl_obj.charge_and_multiplicity_solvent = charge_and_multiplicity_solvent

	### ADD TASK
	tasks = fpl_calc.enthalpy_solvation(fpl_obj, task4)
	fpl_obj.add_task(tasks, parallel=True)

	fpl_obj.start(save=False)

	H_solv = fpl_calc.post_enthalpy_solvation(fpl_obj, unit=unit)

	return H_solv
