import os
import fpl
from fpl_lmp_large import job as lmp_large_job

import time, datetime

####################
solute = "pb2+"
solvent = "THTO"
run_name = "%s_%s" % (solvent, solute)

# Generate initial object
fpl_obj = fpl.fpl_job(run_name, solvent, solute)

# Set simulation parameters for system here
fpl_obj.num_solvents=8

# Generate system
fpl_obj.generate_system()

# Set simulation parameters for lammps here
fpl_obj.queue=None
fpl_obj.procs=1

# Add the task here
fpl_obj.add_task(
	lmp_large_job(fpl_obj)
	)

fpl_obj.start()

####################

