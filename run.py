import fpl_auto

#e_solv = fpl_auto.get_enthalpy_solvation("pb2+","THTO")
#print e_solv

e_solv = fpl_auto.get_enthalpy_solvation("pb2+","THTO",on_queue=True)
e_solv.wait()
H = e_solv.enthalpy(e_solv)
print H