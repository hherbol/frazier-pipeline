# Shared code

# Clancelot Imports
import units

def input_variable(s_id, var, s):
	while s_id in s:
		s = s.replace(s_id, str(var))
	return s

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
			elements_by_index[index] = units.elem_sym_from_weight(float(mass))

	return elements_by_index