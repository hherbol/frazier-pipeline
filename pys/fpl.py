# Code for submission of jobs

# Python Imports
import os, time
import cPickle as pickle

# Clancelot Imports
from fpl_joust import Joust
import structures

# Frazier Pipeline Imports
import fpl_utils
import fpl_constants

class fpl_job(Joust):
	py_on_queue = '''import os, time
import cPickle as pickle

os.system("rm $THIS_PY_FILE$")
fptr = open("$FPL_JOBS_PICKLE$")
fpl_job = pickle.load(fptr)
fpl_job.start(on_queue=None)
'''

	def __init__(self, run_name, solvent_name, solute, ion="Pb",
				 system=None, queue=None, procs=1, mem=1000,
				 priority=100, xhosts=None):
		"""
		The object housing all parameters for simulations, as well as the
		workflow_manager.

		**Parameters**

		**Returns**

			This :class:`fpl.fpl_job` object.
		"""
		Joust.__init__(self, run_name, system=None, queue=None, procs=1,
				 mem=1000, priority=100, xhosts=None)

		# Main system information is here
		self.run_name = run_name
		self.solute = solute
		self.solvent_name = solvent_name
		self.num_solvents = 10

		# Paths are here
		self.HOME_DIR = os.getcwd()+"/"

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
		self.ion = "Pb"

		# DFT Parameters are here
		self.route = "! OPT B97-D3 SV GCP(DFT/TZ) ECP{def2-TZVP} Grid7 SlowConv LooseOpt"
		self.extra_section = fpl_constants.extra_section
		self.grad = False
		self.charge = None
		self.multiplicity = None
		self.charge_and_multiplicity = '0 1'
		self.previous = None   
		self.implicit_solvent = True

		# Other parameters are here
		self.debug = False

	def generate_system(self, halide, cation, ion="Pb"):
		"""
		Generate a system of solute + solvents using Packmol.

		**Returns**

			None
		"""
		## Generate empty system
		system = structures.System(box_size=(25, 25, 25), name=self.run_name)
		## Get structures for solvent and solute
		# Check if upper case or lower case file exists and use the one that does
		if os.path.exists(fpl_constants.cml_dir+self.solvent_name.lower()+".cml"):
			solvent = structures.Molecule(fpl_constants.cml_dir+self.solvent_name.lower()+".cml", extra_parameters=self.extra, allow_errors=True)
		elif os.path.exists(fpl_constants.cml_dir+self.solvent_name.upper()+".cml"):
			solvent = structures.Molecule(fpl_constants.cml_dir+self.solvent_name.upper()+".cml", extra_parameters=self.extra, allow_errors=True)
		else:
			raise Exception("Solvent file %s.cml does not exist in %s.  Ensure you gave the file exists and re-run." % (self.solvent_name, fpl_constants.cml_dir))
		if self.solute is not None:
			fpl_utils.generate_lead_halide_cation(halide, cation, ion=ion)
			solute = structures.Molecule(fpl_constants.cml_dir+self.solute, test_charges=False, allow_errors=True, default_angles=fpl_constants.default_angles)
			system.add(solute)
		## Pack the system
		system.packmol((solvent,), (1,), fpl_constants.solvent[self.solvent_name]["density"], self.seed)

		self.system = system