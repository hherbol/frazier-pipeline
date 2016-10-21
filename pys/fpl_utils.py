# Shared code
import os

# Clancelot Imports
import units
import structures
import geometry
import files
import orca
from constants import PERIODIC_TABLE

# FPL code
import fpl_constants

def input_variable(s_id, var, s):
	while s_id in s:
		s = s.replace(s_id, str(var))
	return s

def get_xyz_elems(run_name):
	data = open("lammps/%s/%s.data" % (run_name, run_name)).read()
	start = data.index('Masses')
	try:
		end = data.index('Pair Coeffs')
	except:
		end = data.index('Bond Coeffs')

	elements_by_index = {}
	for line in data[start:end].splitlines():
		if line and line[0].isdigit():
			index, mass = line.split()
			elements_by_index[index] = units.elem_sym_from_weight(float(mass))

	return elements_by_index

def generate_lead_halide(halide, ion="Pb"):
	PbX = structures.Molecule([structures.Atom(ion,0,0,0)])
	if type(halide) is str:
		halide = [halide, halide, halide]
	vdw = lambda y: PERIODIC_TABLE[units.elem_s2i(y)]['vdw_r']
	for x in halide:
		v = vdw(x)
		PbX.atoms.append(structures.Atom(x,v,0,0.5*v))
		R = geometry.rotation_matrix([0,0,1],120,units="deg")
		PbX.rotate(R)
	return PbX

def generate_lead_halide_cation(halide, cation, ion="Pb", run_opt=True):
	cml_path = fpl_constants.cml_dir
	# Check if system exists
	fname = reduce_to_name(ion, halide, cation)
	if not cml_path.endswith("/"):
		cml_path += "/"

	if os.path.exists(cml_path+fname+".cml"):
		print("Found system in cml folder, returning system")
		system = structures.Molecule(files.read_cml(cml_path+fname+".cml", test_charges=False, allow_errors=True)[0])
		return system

	vdw = lambda y: PERIODIC_TABLE[units.elem_s2i(y)]['vdw_r']
	# Get the PbX3 system
	PbX3 = generate_lead_halide(halide, ion=ion)
	# Get the cation from the cml file
	atoms, bonds, _, _ = files.read_cml(cml_path + cation + ".cml", test_charges=False, allow_errors=True)
	system = structures.Molecule(atoms)
	# Align along X axis
	system.atoms = geometry.align_centroid(system.atoms)[0]
	# Rotate to Z axis
	# NOTE! In case of FA, we want flat so only translate to origin instead
	# NOTE! We have exactly 3 cations we observe: Cs, MA, FA. If 2 N, then FA
	elems = [a.element for a in system.atoms]
	if elems.count("N") == 2:
		system.translate(system.get_center_of_mass())
	else:
		R = geometry.rotation_matrix([0,1,0],90,units="deg")
		system.rotate(R)
	# If N and C in system, ensure N is below C (closer to Pb)
	if "N" in elems and "C" in elems:
		N_index = [i for i,a in enumerate(system.atoms) if a.element == "N"][0]
		C_index = [i for i,a in enumerate(system.atoms) if a.element == "C"][0]
		if system.atoms[N_index].z > system.atoms[C_index].z:
			# Flip if needed
			R = geometry.rotation_matrix([0,1,0],180,units="deg")
			system.rotate(R)
	# Offset system so lowest point is at 0 in the z dir
	z_offset = min([a.z for a in system.atoms])*-1
	system.translate([0,0,z_offset])

	# Add to the PbX3 system with an offset of vdw(Pb)
	system.translate([0,0,vdw(ion)])
	system.atoms += PbX3.atoms
	
	# Run a geometry optimization of this system
	if run_opt:
		PbXY = orca.job(fname,fpl_constants.default_routes[0],
			atoms=system.atoms,
			extra_section=fpl_constants.extra_section,
			queue="batch",
			procs=2)
		PbXY.wait()
		new_pos = orca.read(fname).atoms
		for a,b in zip(system.atoms, new_pos):
			a.x, a.y, a.z = [b.x, b.y, b.z]

	# Set OPLS types
	for a in system.atoms:
		if a.element in [ion,"Cl","Br","I"]:
			a.type = fpl_constants.atom_types[a.element]
			a.type_index = a.type["index"]

	# Write cml file so we don't re-generate, and return system
	files.write_cml(system, bonds=bonds, name=cml_path+fname+".cml")
	return system

def reduce_to_name(i,h,c):
	if type(h) is str:
		fname = i+h+"3"+c
	else:
		hh = [x + str(h.count(x)) if h.count(x) > 1 else x for x in geometry.reduce_list(h)]
		hh.sort()
		hh = "".join(hh)
		fname = i+hh+c
	return fname