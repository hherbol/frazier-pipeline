"""
Automated calculations for the Frazier Pipeline are as follows:

- :func:`get_MBO`
- :func:`get_UMBO`

The following is still prone to bugs.  It should work for num_solvents=1;
however, any more and it is prone to blowing up (ie, irrationally large
enthalpy of solvations).  For now, it has been hidden to avoid incorrect
calculations.

- :func:`get_enthalpy_solvation`

For a quick start guide, in theory the following will suffice:

.. code-block:: python

	import fpl_auto

	halides = [ "Cl", "Br", "I" ]
	cations = ["MA","FA","Cs"]

	# Due to not having all the solvent densities and dielectrics in
	# fpl_constants, only the following are currently available.
	solvents = ["acetone",
                "gbl",
                "DMF",
                "dmso",
                "nmp",
                "THTO"]

	# Note, you can also do mixed halides by passing a list of halides
	# such as h = ["Cl","Cl","I"].  In this example though, we just
	# have h = "Cl" which is equivalent to h = ["Cl","Cl","Cl"]
	h, c, s = halides[0], cations[0], solvents[0]

	# Get the UMBO of the Oxygen-Carbon bond in a single Acetone molecule 
	# adsorbed onto a solute of PbCl3MA
	umbo = fpl_auto.get_UMBO(h,c,s)

------------

"""

import os
import types

import jobs
import files
import geometry

import fpl
import fpl_calc
import fpl_constants
import fpl_utils
from fpl_lmp_large import job as lmp_large_job
from fpl_lmp_small import job as lmp_small_job
from fpl_orca import job as orca_job

def _get_enthalpy_solvation(solute, solvent, num_solvents=1, on_queue=False,
	queue="batch", nprocs=1, xhost=None, unit="kT_300",
	charge_and_multiplicity="0 1", charge_and_multiplicity_solute="0 1",
	charge_and_multiplicity_solvent="0 1", name_append="",
	route_lvls=[1,1,1,1]):
	"""
	This automates the process of calculating the enthalpy of solvation for
	a solvent bonded to a given solute.

	NOTE - DO NOT USE THIS METHOD FOR NOW!

	**Parameters**

		solute : *str*
			A solute for which to be solvated.

		solvent : *str*
			The solvent of the system.

		num_solvents : *int, optional*
			The number of solvents to explicitly simulate.

		on_queue : *bool, optional*
			Whether to run this simulation on the queue.

		queue : *str, optional*
			Which queue to run the simulation on.

		nprocs : *int, optional*
			How many processors to use.

		xhost : *list, str or str, optional*
			A list of processors, or a single processor for which to submit
			the simulations (on the queue).

		unit : *str, optional*
			What units the energy should be returned in.

		charge_and_multiplicity : *str, optional*
			The charge and multiplicity of the full system.

		charge_and_multiplicity_solute : *str, optional*
			The charge and multiplicity of the solute.

		charge_and_multiplicity_solvent : *str, optional*
			The charge and multiplicity of the solvent.

		name_append : *str, optional*
			What to append to the name for these simulations.

		route_lvls : *list, int, optional*
			What level of theory (DFT) to run these calculations at.

	**Return**

		H_solv : *float*
			The enthalpy of solvation.
	"""
	if num_solvents < 1:
		raise Exception("You must have a number of solvents greater than or\
		equal to 1")

	if on_queue:
		pysub_str = """import fpl_auto
e_solv = fpl_auto.get_enthalpy_solvation("$SOLUTE","$SOLVENT", num_solvents=$NUM_SOLVENTS, 
	unit="$UNIT", charge_and_multiplicity="$CHARGE_AND_MULTIPLICITY",
	charge_and_multiplicity_solute="$CHARGE_AND_MULTIPLICITY_SOLUTE",
	charge_and_multiplicity_solvent="$CHARGE_AND_MULTIPLICITY_SOLVENT",
	name_append="$NAME_APPEND", route_lvls=$ROUTE_LVLS)
print e_solv
"""
		pysub_str = pysub_str.replace("$SOLUTE",solute)
		pysub_str = pysub_str.replace("$SOLVENT",solvent)
		pysub_str = pysub_str.replace("$NUM_SOLVENTS",str(num_solvents))
		pysub_str = pysub_str.replace("$UNIT",unit)
		pysub_str = pysub_str.replace("$CHARGE_AND_MULTIPLICITY",charge_and_multiplicity)
		pysub_str = pysub_str.replace("$CHARGE_AND_MULTIPLICITY_SOLUTE",charge_and_multiplicity_solute)
		pysub_str = pysub_str.replace("$CHARGE_AND_MULTIPLICITY_SOLVENT",charge_and_multiplicity_solvent)
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

	# Generate a new system here
	fpl_obj_solvent = fpl.fpl_job(run_name, solvent, None)
	fpl_obj_solvent.cml_dir="/fs/home/hch54/frazier-pipeline/cml/"
	# Set parameters
	fpl_obj_solvent.num_solvents=num_solvents
	# Generate system
	fpl_obj_solvent.generate_system()

	## Task 4 - Large Lammps Simulation of only solvent
	### PARAMETERS
	task1 = run_name + "_large_lammps_solv"
	fpl_obj_solvent.queue=None
	fpl_obj_solvent.procs=1
	fpl_obj_solvent.lmp_run_len=10000
	fpl_obj_solvent.trj_file=None
	### ADD TASK
	task = lmp_large_job(fpl_obj_solvent, task1)
	fpl_obj_solvent.add_task(task)

	## Task 5 - Small Lammps Simulation
	### PARAMETERS
	task2 = run_name + "_small_lammps_solv"
	fpl_obj_solvent.queue=None
	fpl_obj_solvent.procs=1
	fpl_obj_solvent.lmp_run_len=10000
	### ADD TASK
		# Note, you can always overwrite the callback function
		# tsk = lmp_small_job(fpl_obj_solvent, task2)
		# tsk.callback = None
		# fpl_obj_solvent.add_task( task2, tsk )
	task = lmp_small_job(fpl_obj_solvent, task2)
	fpl_obj_solvent.add_task(task)

	# Store solvent system here
	solvent_system = fpl_obj_solvent.system

	########################################

	# Run the simulation here
	fpl_obj_solvent.start(save=False)

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
	tasks = fpl_calc.enthalpy_solvation(fpl_obj, task4, solvent_system)
	fpl_obj.add_task(tasks, parallel=True)

	fpl_obj.start(save=False)

	H_solv = fpl_calc.post_enthalpy_solvation(fpl_obj, unit=unit)

	return H_solv


def get_MBO(halide, cation, solvent,
			ion="Pb",
			num_solvents = 1,
			route_lvls = [0,0,0,0],
			avg=True,
			criteria=[["O","C"],["O","N"],["O","S"]]):
	"""
	Get the mayer bond order.  The Mayer Bond Order (see :func:`get_UMBO`) for more
	details.

	**Parameters**

		halide: *str*
			The halide within the perovskite.

		cation: *str*
			The cation within the perovskite.

		solvent: *str*
			The solvent of the system.

		ion: *str, optional*
			The ion of the perovskite.  By default this is Pb, but can also
			be Sn.

		num_solvents: *int, optional*
			The number of solvents to model explicitly (implicit is always on
			in background).

		route_lvls: *list, int, optional*
			The level of theory to use.

		avg: *bool, optional*
			Whether to average together all UMBO's that match the given
			criteria.

		criteria: *list, list, str, optional*
			A list of lists, each list holding a list describing what bonds
			you want the UMBO for.  By default, it is every bond with an
			oxygen atom involved.

	**Return**

		MBO: *list, float, or float*
			The Mayer Bond Order. If avg is False and there are more than one
			MBO matching criteria, a list is returned.
	"""	
	# Get string for solute
	solute = fpl_utils.reduce_to_name(ion, halide, cation)

	# Get molecules for solute and solvent
	#M_solute = fpl_utils.generate_lead_halide_cation(halide, cation, ion=ion)
	#if os.path.exists(fpl_constants.cml_dir+solvent.lower()+".cml"):
		#M_solvent = files.read_cml(fpl_constants.cml_dir+solvent.lower()+".cml")
	#elif os.path.exists(fpl_constants.cml_dir+solvent.upper()+".cml"):
		#M_solvent = files.read_cml(fpl_constants.cml_dir+solvent.upper()+".cml")
	#else:
		#raise Exception("Solvent file %s.cml does not exist in %s.  Ensure you gave the file exists and re-run." % (solvent, fpl_constants.cml_dir))

	# Setup workflow manager
	################################################################################
	run_name = "%s_%s" % (solute, solvent)
	fpl_obj = fpl.fpl_job(run_name, solvent, solute)
	## PARAMETERS SET HERE #########################################################
	fpl_obj.num_solvents = num_solvents
	fpl_obj.route_lvls = route_lvls
	fpl_obj.charge_and_multiplicity = "0 1"
	fpl_obj.charge_and_multiplicity_solute = "0 1"
	fpl_obj.charge_and_multiplicity_solvent = "0 1"
	fpl_obj.ion = ion
	## Generate the system
	fpl_obj.generate_system(halide, cation, ion=ion)
	########################################
	## Task 1 - Large Lammps Simulation for annealing solvents
	### PARAMETERS
	fpl_obj.queue=None
	fpl_obj.procs=1
	fpl_obj.lmp_run_len=10000
	fpl_obj.trj_file=None
	### ADD TASK
	task = lmp_large_job(fpl_obj, run_name + "_large_lammps")
	fpl_obj.add_task(task)

	## Task 2 - Small Lammps Simulation
	### PARAMETERS
	fpl_obj.queue=None
	fpl_obj.procs=1
	fpl_obj.lmp_run_len=10000
	### ADD TASK
	task = lmp_small_job(fpl_obj, run_name + "_small_lammps")
	fpl_obj.add_task(task)

	## Task 3 - Orca Simulation
	### PARAMETERS
	fpl_obj.queue="batch"
	fpl_obj.procs=4
	fpl_obj.route = route_lvls[0]
	### ADD TASK
	task = orca_job(fpl_obj, run_name + "_orca")
	fpl_obj.add_task(task)
	################################################################################

	# Run the simulation
	fpl_obj.start(save=False)

	# Read in the final results
	MBOs = fpl_obj.data.MBO
	# Find all MBOs for Oxygen and return them
	vals = []
	# Loop through all mbos
	for mbo in MBOs:
		# Loop through criteria
		for crit in criteria:
			# Get a list holding elements in MBO bond as c0,c1
			chk = [a.element for a in mbo[0]]
			# If wildcard for c0, or if c0 is in mbo bond
			if crit[0] == "*" or crit[0] in chk:
				# Remove c0 from mbo bond check (if it's there)
				if crit[0] in chk:
					del chk[chk.index(crit[0])]
				# And check if c1 is wildcard, or in mbo bond
				if crit[1] == "*" or crit[1] in chk:
					# If so, store this mbo
					vals.append(mbo[1])

	if avg:
		vals = sum(vals)/float(len(vals))

	return vals

def get_UMBO(halide, cation , solvent,
			ion="Pb",
			offset=2.0, 
			num_solvents = 1,
			route_lvls = [0,0,0,0],
			avg=True, criteria=[["O","C"],["O","N"],["O","S"]]):
	"""
	Get the unsaturation (average?) mayer bond order.  The Mayer Bond Order
	(MBO) is well described `here <http://pubs.rsc.org/en/Content/ArticleLanding/2001/DT/b102094n#!divAbstract>`_.  In short, it is a numerical representation
	of the probability of how many electrons partake in a bond.  For instance,
	a single bond would have a theoretical bond order of 1.0; however, in
	practice it may have more or less depending on how electrons distribute
	across the molecule.  The MBO helps describe this, and the Unsaturated MBO
	(UMBO) helps represent this in a more understandable fashion.  That is, 
	if the UMBO is larger than zero, the bond is weaker than theory.  If the
	UMBO is less than zero, then the bond is stronger than theory.

	**Parameters**

		halide: *str*
			The halide within the perovskite.

		cation: *str*
			The cation within the perovskite.

		solvent: *str*
			The solvent of the system.

		ion: *str, optional*
			The ion of the perovskite.  By default this is Pb, but can also
			be Sn.

		num_solvents: *int, optional*
			The number of solvents to model explicitly (implicit is always on
			in background).

		route_lvls: *list, int, optional*
			The level of theory to use.

		offset: *float, optional*
			The offset supplied to get the UMBO.  In most cases we consider,
			this is 2.0 as that is the theoretical bond order of a double
			bonded oxygen to sulfur.

		avg: *bool, optional*
			Whether to average together all UMBO's that match the given
			criteria.

		criteria: *list, list, str, optional*
			A list of lists, each list holding a list describing what bonds
			you want the UMBO for.  By default, it is every bond with an
			oxygen atom involved.

	**Return**

		UMBO: *list, float, or float*
			The Unsaturation Mayer Bond Order. If avg is False and there are
			more than one UMBO matching criteria, a list is returned.
	"""
	vals = get_MBO(halide, cation, solvent, ion=ion, avg=avg, num_solvents=num_solvents, route_lvls=route_lvls, criteria=criteria)
	umbos = offset - vals
	return umbos
