# code for packmol simulations
import os, sys, random, math, re, pickle
from units import elem_sym_from_weight
from merlin import *
import numpy

# NOTE!
## There will be two lammps simulations.  One with a large system with the name "run_name"
## and another with a smaller system with the name "run_name_1"

def get_xyz_elems(run_name):
	data = open(run_name+'.data').read()
	start = data.index('Masses')
	try:
		end = data.index('Pair Coeffs')
	except:
		end = data.index('Bond Coeffs')

	elements_by_index = {}
	for line in data[start:end].splitlines():
		if line and line[0].isdigit():
			index, mass = line.split()
			elements_by_index[index] = elem_sym_from_weight(float(mass))

	return elements_by_index

def run(run_name, solvent_name, solute=None, seed=1, LAMMPS_LBL="Pb",  path=os.getcwd()+"/"):

	extra = {  #(47, 3, 46):(85.00, 120.00), (47, 47, 3, 46):(0.0, 14.0, 0.0, 0.0),
		#Pb: utils.Struct(index=Pb, index2=Pb_, element_name='Pb', element=82, mass=207.2, charge=2.0, vdw_e=0., vdw_r=4.0),
		#Cl: utils.Struct(index=Cl, index2=Cl_, element_name='Cl', element=17, mass=35.45, charge=-0.0, vdw_e=0.01, vdw_r=4.0),
	}

	densities = dict( zip(['THTO','gbl', 'DMF','dmso','nmp', 'acetone'], [1.2,1.1,.95,1.0,1.1,.78]) ) 

	# Packmol together
	system = None
	system = utils.System(box_size=(25, 25, 25), name=run_name)
	# NOTE!!!!
	# To not have to swap out Pb and Ba for MD and DFT calculations, we instead have in the CML file
	# the element Pb with label associated with the OPLS line for Ba
	s = "\nDEBUGGING:\n\t%s, %s, %s\n\t%s\n" % (run_name,solvent_name,solute,os.getcwd())
	s = utils.colour_set(s,"RED")
	s = utils.colour_set(s,"BOLD")
	print("%s" % s)

	os.chdir(cwd)
	solvent=utils.Molecule(cwd+'cml/'+solvent_name, extra_parameters=extra, allow_errors=True)
	#for a in solvent.bonds: print("\n\t%s, %s, %s\n\n" % (a.type,a.atoms[0],a.atoms[1]))
	if solute:
		solute = utils.Molecule(cwd+'cml/'+solute, test_charges=False, allow_errors=True)
		system.add(solute)

	files.packmol(system, (solvent,), (1,), densities[solvent_name], seed)
	
	# Run LAMMPs minimization and 1000 timestep simulation to minimize positions -> 
	os.chdir(cwd+'lammps')
	files.write_lammps_data(system,True)
	commands = ('''units real
atom_style full
pair_style lj/cut/coul/dsf 0.05 10.0 10.0
bond_style harmonic
angle_style harmonic
dihedral_style opls

boundary p p p
read_data	'''+run_name+'''.data

dump 1 all xyz 100 '''+run_name+'''.xyz

fix av all ave/time 1 1000 1000 c_thermo_pe
thermo_style custom f_av pe temp press
thermo 1000

group mobile id > '''+str(len(solute.atoms) if solute else 0)+'''

group immobile subtract all mobile
velocity immobile zero linear
fix freeze immobile setforce 0.0 0.0 0.0

minimize 1.0e-4 1.0e-6 100 1000
velocity mobile create 300.0 '''+str(seed)+''' rot yes dist gaussian
velocity immobile set 0.0 0.0 0.0
fix motion mobile npt temp 300.0 300.0 100.0 aniso 1.0 1.0 1000.0
#fix motion mobile nvt temp 300.0 300.0 100.0
timestep 1.0
run 100000

write_restart '''+run_name+'''.restart''')
	open(run_name+'.in', 'w').write(commands)
	os.system('/fs/home/bas348/lammps/15May2015/src/lmp_serial -in '+run_name+'.in -log '+run_name+'.log')

	# Read in output from previous LAMMPs simulation	
	xyz = files.read_xyz(run_name+'.xyz')
	print len(xyz), 'frames'
	# Store end of last LAMMPs simulation to system.atoms variable
	for a,b in zip(system.atoms,xyz[-1]):
		a.x, a.y, a.z = b.x, b.y, b.z
		if any( [numpy.isnan(x) for x in (a.x,a.y,a.z)] ):
			return None
	molecules_in_cluster = []

	if solute:
		pb = None
		for m in system.molecules:
			if LAMMPS_LBL not in [a.element for a in m.atoms]: continue
			pb = m
		if pb is None: raise Exception("Could not find Pb!")
		#pb = [a for a in system.atoms if a.element=='Pb'][0] # Grab solute molecule
		for m in system.molecules: # Get list of molecules and distances from solute molecule
			R = sum([(a-b)**2 for a,b in zip(pb.get_com(),m.get_com())])/3.0
			molecules_in_cluster.append( (R,m) )
	else:
		origin = utils.Atom('X',0.0,0.0,0.0)
		for m in system.molecules:
			R = utils.dist(origin,m.atoms[0])
			molecules_in_cluster.append( (R,m) )

	# Grab closest x solvent molecules.  In the case of solutes existing, we have x+1 in total
	molecules_in_cluster.sort()
	molecules_in_cluster = [m[1] for m in molecules_in_cluster[ : (26 if solute else 25) ]] 
	for j,m in enumerate(molecules_in_cluster):
		for i,a in enumerate(m.atoms):
			a.index = i+1

	# Setup a new system with this reduced size
	system = None
	system = utils.System(box_size=(100, 100, 100), name=run_name+'_1')
	for m in molecules_in_cluster:
		system.add(m)
		for a,b in zip(system.molecules[-1].atoms, m.atoms):
			a.x, a.y, a.z = b.x, b.y, b.z
	files.write_lammps_data(system,True)
	
	# Reoptimize in LAMMPs for the smaller system
	commands = ('''units real
atom_style full
pair_style lj/cut/coul/dsf 0.05 10.0 10.0
bond_style harmonic
angle_style harmonic
dihedral_style opls

boundary p p p
read_data	'''+system.name+'''.data

dump 1 all xyz 100 '''+system.name+'''.xyz

fix av all ave/time 1 100 100 c_thermo_pe
thermo_style custom f_av pe temp press
thermo 100

group mobile id > '''+str(len(solute.atoms) if solute else 0)+'''

group immobile subtract all mobile
velocity immobile zero linear
fix freeze immobile setforce 0.0 0.0 0.0

minimize 1.0e-4 1.0e-6 1000 10000
velocity mobile create 100.0 '''+str(seed)+''' rot yes dist gaussian
velocity immobile set 0.0 0.0 0.0
#fix motion mobile npt temp 300.0 300.0 100.0 aniso 1.0 1.0 1000.0
fix motion mobile nvt temp 100.0 100.0 100.0
timestep 1.0
run 15000
minimize 1.0e-4 1.0e-6 1000 10000

write_restart '''+system.name+'''.restart''')
	open(system.name+'.in', 'w').write(commands)
	os.system('/fs/home/bas348/lammps/15May2015/src/lmp_serial -in '+system.name+'.in -log '+system.name+'.log')

	# Grab last frame of this smaller, newly optimized system and store
	elements_by_index = get_xyz_elems(run_name+"_1")

	xyz = files.read_xyz(system.name+'.xyz')
	print len(xyz), 'frames'
	for a,b in zip(system.atoms,xyz[-1]):
		a.x, a.y, a.z, a.element = b.x, b.y, b.z, elements_by_index[b.element] if elements_by_index[b.element] != "Ba" else "Pb"
		if any( [numpy.isnan(x) for x in (a.x,a.y,a.z)] ):
			return None
	files.write_xyz(system.atoms, run_name+'.cluster')	
	pickle.dump(system, open(run_name+'.pickle', 'w'))

	os.chdir(cwd)

	dielectrics = dict( zip(['THTO','gbl', 'DMF','dmso','nmp', 'acetone'], [42.84,40.24,36.7,46.7,32.2,20.7]) )   
	extra_section='%%cosmo  SMD true  epsilon %f  end  %%basis aux auto NewECP Pb "def2-SD" "def2-TZVP" end NewECP Cs "def2-SD" "def2-TZVP" end NewGTO S "def2-TZVP" end end' % dielectrics[solvent_name]
	#orca.job(run_name, '! OPT B97-D3 SV GCP(DFT/TZ) ECP{def2-TZVP} COSMO Grid7 SlowConv LooseOpt', system.atoms, queue='long', mem=5000, procs=1, extra_section=extra_section, charge=0 if solute else 0)


#run(sys.argv[1], sys.argv[2], None if sys.argv[3]=='None' else sys.argv[3], int(sys.argv[4]) )

cwd = os.getcwd()+"/" # makes sure you are in the correct directory

for i in range(1):
	for solvent in ['acetone']: #,'gbl', 'DMF','dmso','nmp', 'acetone'
		for solute in ['PbI3MA']: #'PbI3FA','PbBr3FA','PbCl3FA','PbI3Cs','PbBr3Cs','PbCl3Cs','PbI3MA','PbBr3MA','PbCl3MA','None'
			name = '%s_%s_%d' % (solute, solvent, i)
			run(name, solvent, None if solute=='None' else solute, seed=i+1, path=cwd)
