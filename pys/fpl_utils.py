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

def generate_lead_halide(halide):
	PbX = structures.Molecule([structures.Atom("Pb",0,0,0)])
	if type(halide) is str:
		halide = [halide, halide, halide]
	vdw = lambda y: PERIODIC_TABLE[units.elem_s2i(y)]['vdw_r']
	for x in halide:
		v = vdw(x)
		PbX.atoms.append(structures.Atom(x,v,0,0.5*v))
		R = geometry.rotation_matrix([0,0,1],120,units="deg")
		PbX.rotate(R)
	return PbX

def generate_lead_halide_cation(halide, cation, run_opt=True):
	cml_path = fpl_constants.cml_dir
	# Check if system exists
	if type(halide) is str:
		fname = "Pb"+halide+"3"+cation
	else:
		fname = "Pb"+"".join([x + str(halide.count(x)) if halide.count(x) > 1 else x for x in geometry.reduce_list(halide)])+cation
	if not cml_path.endswith("/"):
		cml_path += "/"

	if os.path.exists(cml_path+fname+".cml"):
		print("Found system in cml folder, returning system")
		system = structures.Molecule(files.read_cml(cml_path+fname+".cml", test_charges=False, allow_errors=True)[0])
		return system

	vdw = lambda y: PERIODIC_TABLE[units.elem_s2i(y)]['vdw_r']
	# Get the PbX3 system
	PbX3 = generate_lead_halide(halide)
	# Get the cation from the cml file
	system = structures.Molecule(files.read_cml(cml_path + cation + ".cml", test_charges=False, allow_errors=True)[0])
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
	system.translate([0,0,vdw("Pb")])
	system.atoms += PbX3.atoms
	
	# Run a geometry optimization of this system
	if run_opt:
		PbXY = orca.job(fname,"! OPT B97-D3 SV GCP(DFT/TZ) ECP{def2-TZVP} Grid7 SlowConv LooseOpt",
			atoms=system.atoms,
			extra_section="%basis aux auto NewECP Pb \"def2-SD\" \"def2-TZVP\" end NewECP Cs \"def2-SD\" \"def2-TZVP\" end NewGTO S \"def2-TZVP\" end end",
			queue="batch",
			procs=2)
		PbXY.wait()
		system.atoms = orca.read(fname).atoms

	files.write_cml(system, name=cml_path+fname+".cml")

	return system