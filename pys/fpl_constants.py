# Clancelot Imports
import structures
import sysconst

path_to_lammps = sysconst.lmp_path

solvent = {"THTO":{"density":1.2, "dielectric":42.84},
           "GBL":{"density":1.1, "dielectric":40.24},
           "DMF":{"density":0.95, "dielectric":36.7},
           "DMSO":{"density":1.0, "dielectric":46.7},
           "NMP":{"density":1.1, "dielectric":32.3},
           "ACETONE":{"density":0.78, "dielectric":20.7}}
for key in solvent.keys():
  solvent[key.lower()] = solvent[key]

extra = {  #(47, 3, 46):(85.00, 120.00), (47, 47, 3, 46):(0.0, 14.0, 0.0, 0.0),
		#Pb: utils.Struct(index=Pb, index2=Pb_, element_name='Pb', element=82, mass=207.2, charge=2.0, vdw_e=0., vdw_r=4.0),
		#Cl: utils.Struct(index=Cl, index2=Cl_, element_name='Cl', element=17, mass=35.45, charge=-0.0, vdw_e=0.01, vdw_r=4.0),
	}

default_angles = {"type":structures.Struct(),
                  "angle":110.7,
                  "style":"harmonic",
                  "e":37.5,
                  "index2s":(13, 13, 46)}

default_routes = ["! HF SV ECP{def2-TZVP}",
         "! OPT HF SV ECP{def2-TZVP} LooseOpt",
         "! OPT B97-D3 SV GCP(DFT/TZ) ECP{def2-TZVP} Grid7 SlowConv LooseOpt",
         "! OPT B97-D3 SV GCP(DFT/TZ) ECP{def2-TZVP} Grid7 SlowConv LooseOpt"]

cml_dir="/fs/europa/g_pc/fpl/cml/"

atom_types = {
  "Pb" : {'bond_count': 0, 'style': 'lj/cut', 'vdw_r': 3.81661, 'index': 356, 'notes': 'Barium Ion Ba+2', 'element': 56, 'vdw_e': 0.047096, 'charge': 2.0, 'mass': 137.33, 'index2': 76, 'element_name': 'Ba'},
  "Cl" : {'bond_count': 0, 'style': 'lj/cut', 'vdw_r': 4.02, 'index': 344, 'notes': 'Chloride Ion Cl-', 'element': 17, 'vdw_e': 0.71, 'charge': -1.0, 'mass': 35.453, 'index2': 21, 'element_name': 'Cl'},
  "Br" : {'bond_count': 0, 'style': 'lj/cut', 'vdw_r': 4.28, 'index': 345, 'notes': 'Bromide Ion Br-', 'element': 35, 'vdw_e': 0.71, 'charge': -1.0, 'mass': 79.904, 'index2': 65, 'element_name': 'Br'},
  "I" : {'bond_count': 0, 'style': 'lj/cut', 'vdw_r': 4.81, 'index': 346, 'notes': 'Iodide Ion I-', 'element': 53, 'vdw_e': 0.71, 'charge': -1.0, 'mass': 126.905, 'index2': 66, 'element_name': 'I'}
}