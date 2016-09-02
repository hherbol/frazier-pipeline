# Code for small lammps simulation

# Python Imports
import os
import numpy as np

# fpl Imports
import fpl_constants, fpl_utils

# Clancelot Imports
import files, utils

# Requirements:
#     1. The simulation is run from a folder with a subfolder "cml" containing the structure
#        of the solvent in the file "$solvent_name.cml"
#     2. Similarly, if a solute is to be used its structure can be found in "cml/$solute.cml"

# Returns:
#     job() will return None if failure at any point
#
#

# TODO - Generalize where we find the solute. Currently we assume "Pb" is the solute and we find the molecule with the "Pb" element. What if it isn't in the solute though?
def job(run_name, prev_run_name, system, solvent_name, solute=None, seed=1, num_solvents=25, path=os.getcwd()+"/", extra={}, cml_dir=os.getcwd()+"/cml/", debug=False):

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

fix av all ave/time 1 100 100 c_thermo_pe
thermo_style custom f_av pe temp press
thermo 100

group mobile id > $MOBILE$
group immobile subtract all mobile

$IMOBILE$

minimize 1.0e-4 1.0e-6 1000 10000

velocity mobile create 100.0 $SEED$ rot yes dist gaussian
velocity immobile set 0.0 0.0 0.0

fix motion mobile nvt temp 100.0 100.0 100.0

timestep 1.0
#run 15000
run 300

minimize 1.0e-4 1.0e-6 1000 10000

write_restart $RUN_NAME$.restart'''

	# Step 1 - Read in previous large lammps simulation
	xyz = files.read_xyz(prev_run_name+'.xyz')
	
	## Store end of last LAMMPs simulation to system.atoms variable
	for a,b in zip(system.atoms, xyz[-1]):
		a.x, a.y, a.z = b.x, b.y, b.z
		if any( [np.isnan(x) for x in (a.x,a.y,a.z)] ):
			return None

	## Grab only molecules we're interested in.  Here we find relative distances to the solute in question
	molecules_in_cluster = []
	m_solute = None
	if solute:
		m_solute = utils.Molecule(cml_dir+solute, test_charges=False, allow_errors=True)
		diffs = []
		for molec in system.molecules:
			# NOTE, ORDER MATTERS! As procrustes WILL change the atomic positions of the
			# second list of atoms to best match the first.  We don't care if m_solute
			# changes, but if everything else overlaps with m_solute then we have an issue.
			chk = [molec.atoms, m_solute.atoms]
			if len(chk[0]) != len(chk[1]): continue
			#chk = [copy.deepcopy(molec.atoms), copy.deepcopy(m_solute.atoms)]
			utils.procrustes(chk)
			diffs.append(utils.motion_per_frame(chk)[-1])
		index_of_solute = diffs.index(min(diffs))

		m_solute = system.molecules[index_of_solute]

		for m in system.molecules: # Get list of molecules and distances from solute molecule
			R = sum([(a-b)**2 for a,b in zip(m_solute.get_com(skip_H=True),m.get_com(skip_H=True))])/3.0
			molecules_in_cluster.append( (R,m) )
	else:
		origin = utils.Atom('X',0.0,0.0,0.0)
		for m in system.molecules:
			R = utils.dist(origin,m.atoms[0])
			molecules_in_cluster.append( (R,m) )

	## Grab closest x solvent molecules.  In the case of solutes existing, we have x+1 in total
	molecules_in_cluster.sort()
	molecules_in_cluster = [m[1] for m in molecules_in_cluster[ : (num_solvents+1 if solute else num_solvents) ]] 
	for j,m in enumerate(molecules_in_cluster):
		for i,a in enumerate(m.atoms):
			a.index = i+1
	
	## Generate the new system
	system = None
	system = utils.System(box_size=(100, 100, 100), name=run_name)
	for m in molecules_in_cluster:
		system.add(m)
		for a,b in zip(system.molecules[-1].atoms, m.atoms):
			a.x, a.y, a.z = b.x, b.y, b.z
	#files.write_lammps_data(system,True,default_angles=fpl_constants.default_angles)

	# Step 2 - Run LAMMPs minimization for a large solvated box
	os.chdir(path+'lammps')
	files.write_lammps_data(system,True,default_angles=fpl_constants.default_angles)

	mobile = str(len(m_solute.atoms) if m_solute else 0)
	LAMMPS_SIMULATION = fpl_utils.input_variable("$RUN_NAME$", run_name, LAMMPS_SIMULATION)
	LAMMPS_SIMULATION = fpl_utils.input_variable("$MOBILE$", mobile, LAMMPS_SIMULATION)
	LAMMPS_SIMULATION = fpl_utils.input_variable("$SEED$", seed, LAMMPS_SIMULATION)

	imobile = ""
	if solute is not None:
		imobile = "velocity immobile zero linear\nfix freeze immobile setforce 0.0 0.0 0.0"
	LAMMPS_SIMULATION = fpl_utils.input_variable("$IMOBILE$", imobile, LAMMPS_SIMULATION)

	open(run_name+'.in', 'w').write(LAMMPS_SIMULATION)
	if not debug:
		os.system(fpl_constants.path_to_lammps+' -in '+run_name+'.in -log '+run_name+'.log')

	return system
