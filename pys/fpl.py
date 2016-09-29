# Code for submission of jobs

# Python Imports
import os, time
import cPickle as pickle

# Clancelot Imports
from joust import Joust
import structures

# Frazier Pipeline Imports
import fpl_constants

class fpl_job(Joust):
	py_on_queue = '''import os, time
import cPickle as pickle

os.system("rm $THIS_PY_FILE$")
fptr = open("$FPL_JOBS_PICKLE$")
fpl_job = pickle.load(fptr)
fpl_job.start(on_queue=None)
'''

	def __init__(self, run_name, solvent_name, solute, system=None, queue=None, procs=1, mem=1000,
				 priority=100, xhosts=None):
		"""
		The object housing all parameters for simulations, as well as the
		workflow_manager.

		**Parameters**

		**Returns**

			This :class:`fpl.fpl_job` object.
		"""
		Joust.__init__(self, run_name, system=None, queue=None, procs=1, mem=1000,
				 priority=100, xhosts=None)

		# Main system information is here
		self.run_name = run_name
		self.solute = solute
		self.solvent_name = solvent_name
		self.num_solvents = 10

		# Paths are here
		self.HOME_DIR = os.getcwd()+"/"		
		self.cml_dir = os.getcwd()+"/cml/"

		# Simulation parameters are here
		self.queue = None
		self.procs = 1
		self.mem = 1000
		self.priority = 100
		self.xhosts = None
		self.email = None
		self.pair_coeffs_included = True
		self.hybrid_pair = False
		self.hybrid_angle = False
		self.xyz_file = ''
		self.trj_file = ''
		self.read_atoms = True
		self.read_timesteps = True
		self.read_num_atoms = True
		self.read_box_bounds = True

		# LAMMPs parameters are here
		self.lmp_run_len = 2000
		self.seed = 12345
		self.extra = {}

		# DFT Parameters are here
		self.route = "! OPT B97-D3 SV GCP(DFT/TZ) ECP{def2-TZVP} Grid7 SlowConv LooseOpt"
		self.extra_section = "%basis aux auto NewECP Pb \"def2-SD\" \"def2-TZVP\" end NewECP Cs \"def2-SD\" \"def2-TZVP\" end NewGTO S \"def2-TZVP\" end end"
		self.implicit_solvent = True
		self.dft_params = {"queue":"long", "mem":5000, "procs":1, "charge":0, "xhost":None}
		self.pysub_params = {"queue":"batch", "procs":1, "xhost":None}

		# Other parameters are here
		self.debug = False

	def generate_system(self):
		"""
		Generate a system of solute + solvents using Packmol.

		**Returns**

			None
		"""
		## Generate empty system
		system = structures.System(box_size=(25, 25, 25), name=self.run_name)
		## Get structures for solvent and solute
		solvent = structures.Molecule(self.cml_dir+self.solvent_name, extra_parameters=self.extra, allow_errors=True)
		if self.solute is not None:
			solute = structures.Molecule(self.cml_dir+self.solute, test_charges=False, allow_errors=True)
			system.add(solute)
		## Pack the system
		system.packmol((solvent,), (1,), fpl_constants.solvent[self.solvent_name]["density"], self.seed)

		self.system = system