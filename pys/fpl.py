# Code for submission of jobs

# Python Imports
import os, time

# fpl Imports
import fpl_lmp_large, fpl_lmp_small, fpl_orca #, fpl_calc

class fpl_job:
	def __init__(self, run_name, solvent_name, solute=None, seed=1, num_solvents=25, path=os.getcwd()+"/", cml_dir=os.getcwd()+"/cml/", extra={},
                 route = "! OPT B97-D3 SV GCP(DFT/TZ) ECP{def2-TZVP} Grid7 SlowConv LooseOpt",
                 extra_section = "%%basis aux auto NewECP Pb \"def2-SD\" \"def2-TZVP\" end NewECP Cs \"def2-SD\" \"def2-TZVP\" end NewGTO S \"def2-TZVP\" end end",
                 implicit_solvent = True,
                 dft_params = {"queue":"long", "mem":5000, "procs":1, "charge":0},
		         debug=False):
		self.run_name = run_name
		self.solvent_name = solvent_name
		self.solute = solute
		self.seed = seed
		self.num_solvents = num_solvents
		self.path = path
		self.cml_dir = cml_dir
		self.extra = extra
		self.route = route
		self.extra_section = extra_section
		self.dft_params = dft_params
		self.debug = debug

		self.finished = False
		self.enthalpy_of_solvation = None

	def wait_till_finished(self, t=60):
		while self.finished is False:
			time.sleep(60)

	def is_finished(self):
		return self.finished

	def start(self):
		# Run jobs here
		print("\n\t\t\tRUNNING BIG\n")
		system = fpl_lmp_large.job(self.run_name, self.solvent_name, solute=self.solute, seed=self.seed, path=self.path, cml_dir=self.cml_dir, extra=self.extra, debug=self.debug)
		print("\n\t\t\tRUNNING SMALL\n")
		system = fpl_lmp_small.job(self.run_name+"small", system, self.solvent_name, solute=self.solute, seed=self.seed, num_solvents=self.num_solvents, path=self.path, cml_dir=self.cml_dir, extra=self.extra, debug=self.debug)
		print("\n\t\t\tRUNNING DFT\n")
		fpl_orca.job(self.run_name, prev_run_name+"small", system, self.solvent_name, path=self.path,
					route = self.route,
					extra_section = self.extra_section,
					implicit_solvent = self.implicit_solvent,
					dft_params = self.dft_params,
					debug=self.debug)

		# Calculate all output here
		#self.enthalpy_of_solvation = fpl_calc.get_solvation()

		# Done, so specify we're done
		self.finished = True