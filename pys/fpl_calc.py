# System Imports
import copy
import numpy as np

# Clancelot Imports
import structures
import geometry
import units
from joust_task import _jtask
from lammps_job import lmp_task
from orca import orca_task

# A function to get the index of a solute. In the case of no solute, we are calculating
# the self interaction enthalpy of solvation (ie, enthalpy of solvation of molecule x in x).
def _get_solute_index(fpl_obj):
	try:
		xyz = fpl_obj.data[-1]
	except TypeError:
		xyz = [fpl_obj.data.atoms]

	system = fpl_obj.system
	## Store end of last LAMMPs simulation to system.atoms variable
	for a,b in zip(system.atoms, xyz[-1]):
		a.x, a.y, a.z = b.x, b.y, b.z
		if any( [np.isnan(x) for x in (a.x,a.y,a.z)] ):
			return None

	## Grab only molecules we're interested in.  Here we find relative distances to the solute in question
	molecules_in_cluster = []
	m_solute = None
	if fpl_obj.solute:
		m_solute = structures.Molecule(fpl_obj.cml_dir+fpl_obj.solute, test_charges=False, allow_errors=True)
		diffs = []
		for molec in system.molecules:
			# NOTE, ORDER MATTERS! As procrustes WILL change the atomic positions of the
			# second list of atoms to best match the first.  We don't care if m_solute
			# changes, but if everything else overlaps with m_solute then we have an issue.
			chk = [molec.atoms, m_solute.atoms]
			if len(chk[0]) != len(chk[1]): continue
			#chk = [copy.deepcopy(molec.atoms), copy.deepcopy(m_solute.atoms)]
			geometry.procrustes(chk)
			diffs.append(geometry.motion_per_frame(chk)[-1])
		index_of_solute = diffs.index(min(diffs))
	else:
		index_of_solute = 0

	return index_of_solute

def enthalpy_solvation(fpl_obj, task_name, md_dft="dft"):
	# Enthalpy of solvation involves running 3 simulations in parallel
	#  1. Full System
	#  2. System with only Solute
	#  3. System with only Solvent
	# Should be noted that (1) may have been completed beforehand, as we are using
	# an optimized state to run this calculation.
	## NOTE! A second method SHOULD be written to do enthalpy of solvation using MD
	md_dft = md_dft.lower()
	if md_dft not in ["md","dft"]:
		raise Exception("md_dft must be either 'md' or 'dft'")

	full_system = copy.deepcopy(fpl_obj.system)
	solute_system = copy.deepcopy(fpl_obj.system)
	solvent_system = copy.deepcopy(fpl_obj.system)

	# Generate solute only system and solvent only system
	index_of_solute = _get_solute_index(fpl_obj)
	solute_system = solute_system.molecules[index_of_solute]
	solvent_system.Remove(solvent_system.molecules[index_of_solute])

	if md_dft == "md":
		full_task = lmp_task(task_name+"_full", full_system, queue=fpl_obj.queue, procs=fpl_obj.procs,
			mem=fpl_obj.mem, priority=fpl_obj.priority, xhosts=fpl_obj.xhosts, callback=None)
		solute_task = lmp_task(task_name+"_solute", solute_system, queue=fpl_obj.queue, procs=fpl_obj.procs,
			mem=fpl_obj.mem, priority=fpl_obj.priority, xhosts=fpl_obj.xhosts, callback=None)
		solvent_task = lmp_task(task_name+"_solvent", solvent_system, queue=fpl_obj.queue, procs=fpl_obj.procs,
			mem=fpl_obj.mem, priority=fpl_obj.priority, xhosts=fpl_obj.xhosts, callback=None)
	else:
		full_task = orca_task(task_name+"_full", full_system, queue=fpl_obj.queue, procs=fpl_obj.procs,
			mem=fpl_obj.mem, priority=fpl_obj.priority, xhosts=fpl_obj.xhosts, callback=None)
		solute_task = orca_task(task_name+"_solute", solute_system, queue=fpl_obj.queue, procs=fpl_obj.procs,
			mem=fpl_obj.mem, priority=fpl_obj.priority, xhosts=fpl_obj.xhosts, callback=None)
		solvent_task = orca_task(task_name+"_solvent", solvent_system, queue=fpl_obj.queue, procs=fpl_obj.procs,
			mem=fpl_obj.mem, priority=fpl_obj.priority, xhosts=fpl_obj.xhosts, callback=None)

		# Set parameters
		full_task.set_parameters(
			route = fpl_obj.route,
			extra_section = fpl_obj.extra_section,
			charge_and_multiplicity = fpl_obj.charge_and_multiplicity
		)
		solute_task.set_parameters(
			route = fpl_obj.route_solute,
			extra_section = fpl_obj.extra_section_solute,
			charge_and_multiplicity = fpl_obj.charge_and_multiplicity_solute
		)
		solvent_task.set_parameters(
			route = fpl_obj.route_solvent,
			extra_section = fpl_obj.extra_section_solvent,
			charge_and_multiplicity = fpl_obj.charge_and_multiplicity_solvent
		)

		full_task.persist_system=False
		solute_task.persist_system=False
		solvent_task.persist_system=False
	return [full_task, solute_task, solvent_task]

def post_enthalpy_solvation(fpl_obj, md_dft="dft", unit="kT_300"):
	energies = [data.energy for data in fpl_obj.data]
	return units.convert_energy("Ha",unit,energies[0]-energies[1]-energies[2])