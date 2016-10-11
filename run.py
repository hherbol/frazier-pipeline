import fpl_auto
import time
import units
import sys

solute = "pb2+"
solvent = "acetone"

print("RUNNING TEST BENCH FOR %s in %s" % (solute, solvent))
for i in range(1,10):
	t0 = time.time()
	e_solv = fpl_auto.get_enthalpy_solvation(solute,solvent,num_solvents=i,on_queue=True,unit="Ha",charge_and_multiplicity="2 1",charge_and_multiplicity_solvent="2 1",charge_and_multiplicity_solute="2 1")
	e_solv.wait()
	H = e_solv.enthalpy()
	t = time.time()
	print("Enthalpy of Solvation for pb2+ in 10 acetone is %lg Ha = %lg kcal/mol" % (H, units.convert_energy("Ha","kcal/mol",H)))
	print("\tNum Solvents = %d, time to complete = %lg s" % (i, t-t0))
	sys.stdout.flush()
