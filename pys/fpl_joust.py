"""
Job Organizer for User Simulation Tasks (JOUST).

- :class:`Joust`

------------

"""
# Job Organizer for User Simulation Tasks (JOUST)

# Desired Capabilities:
#    Generate list of simulations to run
#    Potential fall back simulations if conditionals are met
#    Read in and write out data in a consistent manner
#    Save and load states (using pickle)

# System Imports
import cPickle as pickle
import dill
# Clancelot Imports
import structures

class Joust(object):
	"""
	The Job Organizer for User Simulation Tasks (JOUST) is the Clancelot
	workflow manager, used for automating the process of simulating
	consecutive jobs.
	"""

	def __init__(self, name, system=None, queue=None, procs=1, mem=1000,
				 priority=100, xhosts=None):
		self.task_list = {}
		self.task_order = []

		self.name = name

		self.system = system
		self.queue = queue
		self.procs = procs
		self.mem = mem
		self.priority = priority
		self.xhosts = xhosts

		self.parameters = structures.Struct()
		self.data = None

	def add_task(self, task, append_to_run_list=True, parallel=False):
		"""
		Append a task to be run. This is an order dependent process and thus,
		tasks added first will run first.  Note, not all tasks added need be run.
		Any task that may be called from another, but only when a specific
		conditional is met, should also be added here.  In those cases, use:

		.. code-block:: python

			append_to_run_list = False

		**Parameters**

			task_name: *str*
				Name of the task to be run.
			task:
				Task object.
			append_to_run_list: *bool, optional*
				Whether this will add the task to the run list, or not.

		**Returns**

			None
		"""
		if type(task) is not list:
			task = [task]

		parallel_list = []
		for t in task:
			task_name = t.task_name
			if task_name in self.task_list:
				raise Exception("Task name already exists. Remove task or rename.")
			self.task_list[task_name] = t
			if append_to_run_list:
				if not parallel:
					self.task_order.append(task_name)
				else:
					parallel_list.append(task_name)
		if parallel:
			self.task_order.append(parallel_list)


	def del_task(self, task_name):
		"""
		Remove a task from the task list.  This will delete the task, as well as
		any instances in which it would have been called by JOUST.  NOTE! This
		does NOT remove the task from other tasks however.

		**Parameters**

			task_name: *str*
				Name of the task to be removed.
		"""
		if task_name not in self.task_list:
			raise Exception("Task not in list.")
		del self.task_list[task_name]
		ii = [i for i,t in enumerate(self.task_order) if t == task_name][::-1]
		for i in ii:
			del self.task_order[i]

	def save_state(self, fname=None):
		"""
		Save the current workflow state.

		**Parameters**

			fname: *str, optional*
				The name of the pickle file to save.

		**Returns**

			None
		"""
		if fname is None:
			fname = self.name
		if not fname.endswith(".pickle"):
			fname += ".pickle"
		dill.dump(self, open(fname, "wb"))

	def load_state(self, fname):
		"""
		Loads a workflow state.

		**Parameters**

			fname: *str*
				The name of the pickle file to load.

		**Returns**

			None
		"""
		self = dill.load(open(fname,"rb"))

	def start(self, save=False):
		"""
		Start the workflow manager.

		**Parameters**

			save: *bool, optional*
				Whether to save the state of the workflow manager at each
				iteration.

		**Returns**

			None
		"""
		save_counter = 0
		while len(self.task_order) > 0:
			# If saving at each state, do so
			if save:
				save_counter += 1
				self.save_state("%s_%d" % (self.name, save_counter))

			# Get the task to run, set it up, and run it
			task = self.task_order[0]
			## In the case of a sublist, we'll run all in parallel
			if type(task) is list:
				running_jobs = []
				job_handles = []
				for sub_task in task:
					running_jobs.append(self.task_list[sub_task])
					## Setup
					if running_jobs[-1].persist_system:
						running_jobs[-1].system = self.system
					running_jobs[-1].system.name = running_jobs[-1].task_name

					## Run
					job_handles.append(running_jobs[-1].run())
				for j in job_handles:
					j.wait()

				# Read in the data
				self.data = []
				for j in running_jobs:
					j.read_results()
					self.data.append(j.data)

				# Check for conditionals
				conditional_jobs = []
				for j in running_jobs:
					if j.conditional(j.data):
						conditional_jobs.append(j.conditional_sim_name)
				if len(conditional_jobs) > 0:
					if len(conditional_jobs) == 1:
						conditional_jobs = conditional_jobs[0]
					self.task_order[0] = conditional_jobs
					continue

				# Check callbacks
				for j in running_jobs:
					if j.callback is not None:
						j.callback(self)

				# Remove the last simulations
				del self.task_order[0]
			else:
				running_job = self.task_list[task]
				## Setup
				if running_job.persist_system:
					running_job.system = self.system
				running_job.system.name = running_job.task_name
				## Run
				job_handle = running_job.run()

				job_handle.wait()

				# Read in the results of the simulation
				running_job.read_results()

				# If we have a conditional simulation to run, check and do so
				# Note, in the case of a conditional, callback is not run!
				if running_job.conditional(running_job.data):
					self.task_order[0] = running_job.conditional_sim_name
					self.data = running_job.data
					continue

				# Store the data from the last simulation here
				self.data = running_job.data

				if running_job.callback is not None:
					running_job.callback(self)
				
				# Else, remove the finished simulation and continue
				del self.task_order[0]

			if save:
				save_counter += 1
				self.save_state("%s_%d" % (self.name, save_counter))
