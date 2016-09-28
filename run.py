import os
import fpl

import time, datetime

main_dir = os.getcwd()

####################
solute = "pb2+"
for solvent in ["THTO"]:
	os.chdir(main_dir)
	run_name = "debug_2_"+solvent

	if not os.path.exists(run_name):
		os.mkdir(run_name)
	os.chdir(run_name)

	run_simulation = fpl.fpl_job(run_name, solvent, solute)

	run_simulation.cml_dir = "/home/hherbol/frazier-pipeline/cml/"
	run_simulation.path = os.getcwd()
	run_simulation.dft_params["queue"]=None
	run_simulation.dft_params["procs"]=2
	run_simulation.num_solvents=1
	run_simulation.bl_rl = 2000
	run_simulation.ll_rl = 2000

	run_simulation.start()
	#run_simulation.start(on_queue=True)

#print("Started at %s" % datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
#run_simulation.wait_till_finished()
#print("Ended at %s" % datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
####################

