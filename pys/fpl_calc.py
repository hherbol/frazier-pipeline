# System Imports
import copy
import numpy as np

# FPL imports
import fpl_constants
import fpl_utils

# Clancelot Imports
import structures
import geometry
import units
import lammps_job
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

def _minimize_solvent(fpl_obj, run_name):

	input_script = '''units real
atom_style full
pair_style lj/cut/coul/dsf 0.05 10.0 10.0
bond_style harmonic
angle_style harmonic
dihedral_style opls

boundary p p p
read_data $RUN_NAME$.data

dump 1 all xyz 100 $RUN_NAME$.xyz

fix av all ave/time 1 1000 1000 c_thermo_pe
thermo_style custom step f_av pe temp press
thermo 1000

minimize 1.0e-4 1.0e-6 100 1000

velocity all create 300.0 $SEED$ rot yes dist gaussian

fix motion all npt temp 300.0 300.0 100.0 aniso 1.0 1.0 1000.0

timestep 1.0
run $RUN_LEN$

write_restart $RUN_NAME$.restart'''
	
	# Setup input script
	input_script = fpl_utils.input_variable("$RUN_NAME$", run_name, input_script)
	input_script = fpl_utils.input_variable("$SEED$", fpl_obj.seed, input_script)
	input_script = fpl_utils.input_variable("$RUN_LEN$", fpl_obj.lmp_run_len, input_script)

	## Generate empty system
	system = structures.System(box_size=(25, 25, 25), name=fpl_obj.run_name)
	## Get structures for solvent and solute
	solvent = structures.Molecule(fpl_obj.cml_dir+fpl_obj.solvent_name, extra_parameters=fpl_obj.extra, allow_errors=True)

	## Pack the system
	system.packmol((solvent,), (1,), fpl_constants.solvent[fpl_obj.solvent_name]["density"], fpl_obj.seed, number=fpl_obj.num_solvents)
	system.name = run_name

	# Run simulation
	running_job = lammps_job.job(run_name, input_script, system, queue=None, procs=1, email=None, pair_coeffs_included=True, hybrid_pair=False, hybrid_angle=False, TIP4P=False)
	running_job.wait()

	# Read in data
	data = lammps_job.read(run_name, trj_file=fpl_obj.trj_file,
			   xyz_file=fpl_obj.xyz_file,
			   read_atoms=fpl_obj.read_atoms,
			   read_timesteps=fpl_obj.read_timesteps,
			   read_num_atoms=fpl_obj.read_num_atoms,
			   read_box_bounds=fpl_obj.read_box_bounds)
	xyz = data[-1]

	## Store end of last LAMMPs simulation to system.atoms variable
	for a,b in zip(system.atoms, xyz[-1]):
		a.x, a.y, a.z = b.x, b.y, b.z
		if any( [np.isnan(x) for x in (a.x,a.y,a.z)] ):
			return None

	return system

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

	print("Starting solvent minimization")
	solvent_system = _minimize_solvent(fpl_obj, task_name+"_pre_solv")
	print("Ending solvent minimization")

	if md_dft == "md":
		full_task = lmp_task(task_name+"_full", full_system, queue=fpl_obj.queue, procs=fpl_obj.procs,
			mem=fpl_obj.mem, priority=fpl_obj.priority, xhosts=fpl_obj.xhosts, callback=None)
		solute_task = lmp_task(task_name+"_solute", solute_system, queue=fpl_obj.queue, procs=fpl_obj.procs,
			mem=fpl_obj.mem, priority=fpl_obj.priority, xhosts=fpl_obj.xhosts, callback=None)
		solvent_task = lmp_task(task_name+"_solvent", solvent_system, queue=fpl_obj.queue, procs=fpl_obj.procs,
			mem=fpl_obj.mem, priority=fpl_obj.priority, xhosts=fpl_obj.xhosts, callback=None)
	else:
		print("Generating DFT tasks for enthalpy calc")
		full_task = orca_task(task_name+"_full", full_system, queue=fpl_obj.queue, procs=fpl_obj.procs,
			mem=fpl_obj.mem, priority=fpl_obj.priority, xhosts=fpl_obj.xhosts, callback=None)
		solute_task = orca_task(task_name+"_solute", solute_system, queue=fpl_obj.queue, procs=fpl_obj.procs,
			mem=fpl_obj.mem, priority=fpl_obj.priority, xhosts=fpl_obj.xhosts, callback=None)
		solvent_task = orca_task(task_name+"_solvent", solvent_system, queue=fpl_obj.queue, procs=fpl_obj.procs,
			mem=fpl_obj.mem, priority=fpl_obj.priority, xhosts=fpl_obj.xhosts, callback=None)

		# Set parameters
		if fpl_obj.implicit_solvent:
			fpl_obj.route += " COSMO"
			fpl_obj.extra_section = "%cosmo SMD true epsilon " + str(fpl_constants.solvent[fpl_obj.solvent_name]["dielectric"]) + " end " + fpl_obj.extra_section
			fpl_obj.route_solvent += " COSMO"
			fpl_obj.extra_section_solvent = "%cosmo SMD true epsilon " + str(fpl_constants.solvent[fpl_obj.solvent_name]["dielectric"]) + " end " + fpl_obj.extra_section

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
