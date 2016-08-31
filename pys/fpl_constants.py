
path_to_lammps = "/fs/home/bas348/lammps/15May2015/src/lmp_serial"

solvent = {"THTO":{"density":1.2, "dielectric":42.84},
           "GBL":{"density":1.1, "dielectric":40.24},
           "DMF":{"density":0.95, "dielectric":36.7},
           "DMSO":{"density":1.0, "dielectric":46.7},
           "NMP":{"density":1.1, "dielectric":32.3},
           "ACETONE":{"density":0.78, "dielectric":20.7}}

extra = {  #(47, 3, 46):(85.00, 120.00), (47, 47, 3, 46):(0.0, 14.0, 0.0, 0.0),
		#Pb: utils.Struct(index=Pb, index2=Pb_, element_name='Pb', element=82, mass=207.2, charge=2.0, vdw_e=0., vdw_r=4.0),
		#Cl: utils.Struct(index=Cl, index2=Cl_, element_name='Cl', element=17, mass=35.45, charge=-0.0, vdw_e=0.01, vdw_r=4.0),
	}