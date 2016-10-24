# Code for large lammps simulation

# Python Imports
import os
import numpy as np

# fpl Imports
import fpl_constants
import fpl_utils

# Clancelot Imports
import files
import structures
import geometry
import lammps_job

# Requirements:
#     1. The simulation is run from a folder with a subfolder "cml" containing the structure
#        of the solvent in the file "$solvent_name.cml"
#     2. Similarly, if a solute is to be used its structure can be found in "cml/$solute.cml"
#

# Returns:
#     job() will return the system generated
#
#

def callback_strip_solvents(fpl_obj):
	xyz = fpl_obj.data[-1]
	system = fpl_obj.system
	## Store end of last LAMMPs simulation to system.atoms variable
	for a,b in zip(system.atoms, xyz[-1]):
		a.x, a.y, a.z = b.x, b.y, b.z
		if any( [np.isnan(x) for x in (a.x,a.y,a.z)] ):
			return None

	## Grab only molecules we're interested in.  Here we find relative distances to the solute in question
	if fpl_obj.solute:
		molecules_in_cluster = []

		# Generate a list of all ion atoms
		ion_list = [a for m in system.molecules for a in m.atoms if a.element == fpl_obj.ion]
		all_solutes = [m for m in system.molecules if fpl_obj.ion in [a.element for a in m.atoms]]
		if len(ion_list) == 0: raise Exception("Could not find %s!" % fpl_obj.ion)

		# Get distance between every ion and oxygen per molecule, saving molecule and min r
		for m in system.molecules:
			# Calculate the distance from the oxygen atom to all pb atoms.
			oxypos = [a for a in m.atoms if a.type.element in [8,"O"]]
			if oxypos == []: continue

			# Get a list of the closest oxygen of molecule m to the ions
			r_min = min([geometry.dist(i,o) for i in ion_list for o in oxypos])
			
			# If close enough, save
			if r_min < fpl_obj.R_cutoff:
				molecules_in_cluster.append((m, r_min))
	else:
		all_solutes = []
		origin = structures.Atom('X',0.0,0.0,0.0)
		molecules_in_cluster = [(m, geometry.dist(origin,m.atoms[0])) for m in system.molecules]

	# Now, if we want only N solvents, grab closest N
	molecules_in_cluster = sorted(molecules_in_cluster, key=lambda x: x[1])
	if fpl_obj.num_solvents is not None:
		molecules_in_cluster = all_solutes + [m[0] for m in molecules_in_cluster[:fpl_obj.num_solvents]]
	else:
		molecules_in_cluster = all_solutes + [m[0] for m in molecules_in_cluster]

	# Set indices
	for j,m in enumerate(molecules_in_cluster):
		for i,a in enumerate(m.atoms):
			a.index = i+1
	
	## Generate the new system
	system = None
	system = structures.System(box_size=(100, 100, 100), name=fpl_obj.run_name)
	for m in molecules_in_cluster:
		system.add(m)

	fpl_obj.system = system

def job(fpl_obj, task_name):

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

group mobile id > $MOBILE$
group immobile subtract all mobile

$IMOBILE$

velocity mobile create 300.0 $SEED$ rot yes dist gaussian
velocity immobile set 0.0 0.0 0.0

fix relax mobile nve/limit 0.1
run 10000
unfix relax

fix motion mobile npt temp 300.0 300.0 100.0 iso 1.0 1.0 1000.0

timestep 1.0
run $RUN_LEN$

write_restart $RUN_NAME$.restart'''
	
	# Setup input script
	solute = None
	if fpl_obj.solute is not None:
		solute = structures.Molecule(fpl_constants.cml_dir+fpl_obj.solute, test_charges=False, allow_errors=True)
	mobile = str(len(solute.atoms) if solute else 0)
	input_script = fpl_utils.input_variable("$MOBILE$", mobile, input_script)

	input_script = fpl_utils.input_variable("$RUN_NAME$", task_name, input_script)
	input_script = fpl_utils.input_variable("$SEED$", fpl_obj.seed, input_script)
	input_script = fpl_utils.input_variable("$RUN_LEN$", fpl_obj.lmp_run_len, input_script)

	imobile = ""
	if solute is not None:
		imobile = "velocity immobile zero linear\nfix freeze immobile setforce 0.0 0.0 0.0"
	input_script = fpl_utils.input_variable("$IMOBILE$", imobile, input_script)

	# Now we can generate the task
	# NOTE! Because the data file is written by the system name, we want to overwrite the system name here
	fpl_obj.system.name = task_name

	large_lammps_task = lammps_job.lmp_task(task_name,
		fpl_obj.system, queue=fpl_obj.queue, procs=fpl_obj.procs,
		priority=fpl_obj.priority, xhosts=fpl_obj.xhosts)

	large_lammps_task.set_parameters(input_script, email=fpl_obj.email,
		pair_coeffs_included=fpl_obj.pair_coeffs_included,
		hybrid_pair=fpl_obj.hybrid_pair, hybrid_angle=fpl_obj.hybrid_angle,
		trj_file=fpl_obj.trj_file, xyz_file=fpl_obj.xyz_file, 
		read_atoms=fpl_obj.read_atoms, 
		read_timesteps=fpl_obj.read_timesteps,
		read_num_atoms=fpl_obj.read_num_atoms,
		read_box_bounds=fpl_obj.read_box_bounds)

	large_lammps_task.callback = callback_strip_solvents

	return large_lammps_task
