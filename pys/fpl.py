# Code for submission of jobs

# Python Imports
import os, time
import cPickle as pickle

# fpl Imports
import fpl_lmp_large, fpl_lmp_small, fpl_orca, fpl_utils #, fpl_calc

# Clancelot Imports
import log, utils

class fpl_job:
	py_on_queue = '''import os, time
import cPickle as pickle

import fpl_lmp_large, fpl_lmp_small, fpl_orca

import log

os.system("rm $THIS_PY_FILE$")
fptr = open("$FPL_JOBS_PICKLE$")
fpl_job = pickle.load(fptr)
fpl_job.start(on_queue=None)
'''

	def __init__(self, run_name, solvent_name, solute=None, seed=1, num_solvents=25, path=os.getcwd()+"/", cml_dir=os.getcwd()+"/cml/", extra={},
                 route = "! OPT B97-D3 SV GCP(DFT/TZ) ECP{def2-TZVP} Grid7 SlowConv LooseOpt",
                 extra_section = "%basis aux auto NewECP Pb \"def2-SD\" \"def2-TZVP\" end NewECP Cs \"def2-SD\" \"def2-TZVP\" end NewGTO S \"def2-TZVP\" end end",
                 implicit_solvent = True,
                 dft_params = {"queue":"long", "mem":5000, "procs":1, "charge":0, "xhost":None},
                 pysub_params = {"queue":"batch", "procs":1, "xhost":None},
		         debug=False):
		self.run_name = run_name
		self.solvent_name = solvent_name
		self.solute = solute
		self.seed = seed
		self.num_solvents = num_solvents

		self.path = path
		if not self.path.endswith("/"): self.path += "/"
		self.cml_dir = cml_dir
		if not self.cml_dir.endswith("/"): self.cml_dir += "/"

		self.extra = extra
		self.route = route
		self.extra_section = extra_section
		self.implicit_solvent = implicit_solvent
		self.dft_params = dft_params
		self.pysub_params = pysub_params
		self.debug = debug

		self.finished = False
		self.enthalpy_of_solvation = None

	def wait_till_finished(self, t=60):
		while self.run_name in log.get_jlist():
			time.sleep(60)
		while self.run_name+"_dft" in log.get_jlist():
			time.sleep(60)
		self.finished = True

	def is_finished(self):
		jlist = log.get_jlist()
		self.finished = (self.run_name+"_dft" not in jlist) and (self.run_name not in jlist)
		return self.finished

	def start(self, on_queue=False):
		if not self.path.endswith("/"): self.path += "/"

		if not on_queue:
			# Run jobs here
			print("\n\t\t\tRUNNING BIG\n")
			system = fpl_lmp_large.job(self.run_name+"_large", self.solvent_name, solute=self.solute, seed=self.seed, path=self.path, cml_dir=self.cml_dir, extra=self.extra, debug=self.debug)
			print("\n\t\t\tRUNNING SMALL\n")
			system = fpl_lmp_small.job(self.run_name+"_small", self.run_name+"_large", system, self.solvent_name, solute=self.solute, seed=self.seed, num_solvents=self.num_solvents, path=self.path, cml_dir=self.cml_dir, extra=self.extra, debug=self.debug)
			print("\n\t\t\tRUNNING DFT\n")
			fpl_orca.job(self.run_name+"_dft", self.run_name+"_small", system, self.solvent_name, path=self.path,
						route = self.route,
						extra_section = self.extra_section,
						implicit_solvent = self.implicit_solvent,
						dft_params = self.dft_params,
						debug=self.debug)
		else:
			if not os.path.exists(self.path+"fpl_jobs"):
				os.mkdir(self.path+"fpl_jobs")
			f_pickle = self.path+"fpl_jobs/"+self.run_name+".pickle"
			py_file = self.path+self.run_name+".py"
			pickle.dump(self, open(f_pickle, 'w'))
			self.py_on_queue = fpl_utils.input_variable("$FPL_JOBS_PICKLE$", f_pickle, self.py_on_queue)
			self.py_on_queue = fpl_utils.input_variable("$THIS_PY_FILE$", py_file, self.py_on_queue)
			f_py_file = open(py_file, 'w')
			f_py_file.write(self.py_on_queue)
			f_py_file.close()
			utils.pysub(self.run_name, nprocs=self.pysub_params["procs"], queue=self.pysub_params["queue"], xhost=self.pysub_params["xhost"], path=os.getcwd(), remove_nbs=True)

		# Calculate all output here
		#self.enthalpy_of_solvation = fpl_calc.get_solvation()

		# Done, so specify we're done
		