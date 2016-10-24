# Code for orca's dft simulation

# Python Imports
import os
import numpy as np
import cPickle as pickle

# fpl Imports
import fpl_constants, fpl_utils

# Clancelot Imports
import files, orca

def job(fpl_obj, task_name):
	# If route is specified by "level of theory index" take care of it here
	if type(fpl_obj.route) is int:
		fpl_obj.route = fpl_constants.default_routes[fpl_obj.route]
		#raise Exception("Currently all 0, 1, and 2 are the same (lvl 0, lowest)")

	if fpl_obj.implicit_solvent:
		if "COSMO" not in fpl_obj.route:
			fpl_obj.route += " COSMO"
		if "cosmo SMD true epsilon" not in fpl_obj.extra_section:
			fpl_obj.extra_section = "%cosmo SMD true epsilon " + str(fpl_constants.solvent[fpl_obj.solvent_name]["dielectric"]) + " end " + fpl_obj.extra_section

	if fpl_obj.debug:
		return None

	fpl_obj.system.name = task_name

	orca_task = orca.orca_task(task_name,
		fpl_obj.system, queue=fpl_obj.queue, procs=fpl_obj.procs, mem=fpl_obj.mem,
		priority=fpl_obj.priority, xhosts=fpl_obj.xhosts)

	orca_task.set_parameters(route=fpl_obj.route,extra_section=fpl_obj.extra_section, grad=fpl_obj.grad,
					   charge=fpl_obj.charge, multiplicity=fpl_obj.multiplicity,
					   charge_and_multiplicity=fpl_obj.charge_and_multiplicity, previous=fpl_obj.previous
		)

	#orca_task.callback = callback_grab_final

	return orca_task