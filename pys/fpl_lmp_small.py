# Code for small lammps simulation

# Python Imports
import numpy as np

# fpl Imports
import fpl_utils
import fpl_constants

# Clancelot Imports
import structures
import lammps_job

# Requirements:
#     1. The simulation is run from a folder with a subfolder "cml" containing the structure
#        of the solvent in the file "$solvent_name.cml"
#     2. Similarly, if a solute is to be used its structure can be found in "cml/$solute.cml"

# Returns:
#     job() will return None if failure at any point
#
#

def callback_grab_final(fpl_obj):
	elements_by_index = fpl_utils.get_xyz_elems(fpl_obj.task_order[0])
	xyz = fpl_obj.data[-1]
	## Convert Ba to Pb (we had used Ba in LAMMPs, but in reality we want Pb)
	for a,b in zip(fpl_obj.system.atoms,xyz[-1]):
		a.x, a.y, a.z, a.element = b.x, b.y, b.z, elements_by_index[b.element] if elements_by_index[b.element] != "Ba" else fpl_obj.ion
		if any( [np.isnan(x) for x in (a.x,a.y,a.z)] ):
			raise Exception("Error on small lammps callback.")

# TODO - Generalize where we find the solute. Currently we assume "Pb" is the solute and we find the molecule with the "Pb" element. What if it isn't in the solute though?
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

fix av all ave/time 1 100 100 c_thermo_pe
thermo_style custom step f_av pe temp press
thermo 100

group mobile id > $MOBILE$
group immobile subtract all mobile

$IMOBILE$

minimize 1.0e-4 1.0e-6 1000 10000

velocity mobile create 100.0 $SEED$ rot yes dist gaussian
velocity immobile set 0.0 0.0 0.0

fix motion mobile nvt temp 100.0 100.0 100.0

timestep 1.0
run $RUN_LEN$
#run 300

minimize 1.0e-4 1.0e-6 1000 10000

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

	small_lammps_task = lammps_job.lmp_task(task_name,
		fpl_obj.system, queue=fpl_obj.queue, procs=fpl_obj.procs,
		priority=fpl_obj.priority, xhosts=fpl_obj.xhosts)

	small_lammps_task.set_parameters(input_script, email=fpl_obj.email,
		pair_coeffs_included=fpl_obj.pair_coeffs_included,
		hybrid_pair=fpl_obj.hybrid_pair, hybrid_angle=fpl_obj.hybrid_angle,
		trj_file=fpl_obj.trj_file, xyz_file=fpl_obj.xyz_file, 
		read_atoms=fpl_obj.read_atoms, 
		read_timesteps=fpl_obj.read_timesteps,
		read_num_atoms=fpl_obj.read_num_atoms,
		read_box_bounds=fpl_obj.read_box_bounds)

	small_lammps_task.callback = callback_grab_final

	return small_lammps_task