# Code for large lammps simulation

# Python Imports
import os

# fpl Imports
import fpl_constants, fpl_utils

# Clancelot Imports
import files, structures

# Requirements:
#     1. The simulation is run from a folder with a subfolder "cml" containing the structure
#        of the solvent in the file "$solvent_name.cml"
#     2. Similarly, if a solute is to be used its structure can be found in "cml/$solute.cml"
#

# Returns:
#     job() will return the system generated
#
#

def job(run_name, solvent_name, solute=None, seed=1, run_len=10000, path=os.getcwd()+"/", extra={}, cml_dir=os.getcwd()+"/cml/", debug=False):

	# Step 0 - Ensure proper variables
	if not os.path.exists(path+"lammps"):
		os.mkdir(path+"lammps")
	if not path.endswith("/"): path += "/"
	LAMMPS_SIMULATION = '''units real
atom_style full
pair_style lj/cut/coul/dsf 0.05 10.0 10.0
bond_style harmonic
angle_style harmonic
dihedral_style opls

boundary p p p
read_data $RUN_NAME$.data

dump 1 all xyz 100 $RUN_NAME$.xyz

fix av all ave/time 1 1000 1000 c_thermo_pe
thermo_style custom f_av pe temp press
thermo 1000

group mobile id > $MOBILE$
group immobile subtract all mobile

$IMOBILE$

minimize 1.0e-4 1.0e-6 100 1000

velocity mobile create 300.0 $SEED$ rot yes dist gaussian
velocity immobile set 0.0 0.0 0.0

fix motion mobile npt temp 300.0 300.0 100.0 aniso 1.0 1.0 1000.0

timestep 1.0
run $RUN_LEN$
#run 5000

write_restart $RUN_NAME$.restart'''

	# Step 1 - Use Packmol to randomly solvate the solute
	os.chdir(path)
	## Generate empty system
	system = structures.System(box_size=(25, 25, 25), name=run_name)
	## Get structures for solvent and solute
	solvent = structures.Molecule(cml_dir+solvent_name, extra_parameters=extra, allow_errors=True)
	if solute:
		solute = structures.Molecule(cml_dir+solute, test_charges=False, allow_errors=True)
		system.add(solute)
	## Pack the system
	system.packmol((solvent,), (1,), fpl_constants.solvent[solvent_name]["density"], seed)
	
	# Step 2 - Run LAMMPs minimization for a large solvated box
	os.chdir(path+'lammps')
	files.write_lammps_data(system,True,default_angles=fpl_constants.default_angles)

	mobile = str(len(solute.atoms) if solute else 0)
	LAMMPS_SIMULATION = fpl_utils.input_variable("$RUN_NAME$", run_name, LAMMPS_SIMULATION)
	LAMMPS_SIMULATION = fpl_utils.input_variable("$MOBILE$", mobile, LAMMPS_SIMULATION)
	LAMMPS_SIMULATION = fpl_utils.input_variable("$SEED$", seed, LAMMPS_SIMULATION)
	LAMMPS_SIMULATION = fpl_utils.input_variable("$RUN_LEN$", run_len, LAMMPS_SIMULATION)

	imobile = ""
	if solute is not None:
		imobile = "velocity immobile zero linear\nfix freeze immobile setforce 0.0 0.0 0.0"
	LAMMPS_SIMULATION = fpl_utils.input_variable("$IMOBILE$", imobile, LAMMPS_SIMULATION)

	open(run_name+'.in', 'w').write(LAMMPS_SIMULATION)
	if not debug:
		os.system(fpl_constants.path_to_lammps+' -in '+run_name+'.in -log '+run_name+'.log')

	return system
