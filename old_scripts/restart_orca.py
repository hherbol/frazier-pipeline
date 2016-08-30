# for each converged simulation in orca, restars another orca simulation
from merlin import *


def restart_on_TZ_():
	for i in range(18,19):
		for solvent,dielectric in zip(['THTO'],[42.84]):
			for solute in ['None', 'pb2+']:
				extra_section='%%cosmo  SMD true  epsilon %f  end' % dielectric
				old_name = '%s_%s_%d' % (solute,solvent,i)
				new_name = '%s_%s_tz_%d' % (solute,solvent,i)
				if os.path.exists('SMRFF/orca/%s/%s.out'%(old_name, old_name)):
					orca.job(new_name, '! Opt  B97-D3 def2-TZVP GCP(DFT/TZ) ECP{def2-TZVP} COSMO Grid3 FinalGrid5 SlowConv LooseOpt printbasis', previous='/fs/home/bas348/perovskites/SMRFF/orca/%s/%s.orca.gbw'%(old_name, old_name),queue='long',mem=3000,procs=12, priority=200, extra_section=extra_section)
				elif os.path.exists('orca/%s/gbw_files/%s.orca.out'%(old_name, old_name)):
					orca.job(new_name, '! Opt B97-D3 def2-TZVP GCP(DFT/TZ) ECP{def2-TZVP} COSMO Grid3 FinalGrid5 SlowConv LooseOpt printbasis', previous='/fs/home/bas348/perovskites/SMRFF/orca/%s/gbw_files/%s.orca.out'%(old_name, old_name),queue='long',mem=3000,procs=12,priority=200, extra_section=extra_section)


def TZ_new_func():
	for i in range(20):
		for solvent,dielectric in zip(['THTO'],[42.84]):
			for solute in ['None','pb2+']:
				extra_section='%%cosmo  SMD true  epsilon %f  end' % dielectric
				old_name = '%s_%s_tz2_%d' % (solute,solvent,i)
				new_name = '%s_%s_tz3_%d' % (solute,solvent,i)
				if os.path.exists('SMRFF/orca/%s/%s.out'%(old_name, old_name)):
					orca.job(new_name, '! Opt RIJCOSX PW6B95 D3BJ def2-TZVP GCP(DFT/TZ) ECP{def2-TZVP} COSMO Grid3 FinalGrid5 SlowConv LooseOpt printbasis', previous='/fs/home/bas348/perovskites/SMRFF/orca/%s/%s.orca.gbw'%(old_name, old_name),queue='long',mem=5000,procs=4, priority=250, extra_section=extra_section)
				elif os.path.exists('orca/%s/gbw_files/%s.orca.out'%(old_name, old_name)):
					orca.job(new_name, '! Opt RIJCOSX PW6B95 D3BJ def2-TZVP GCP(DFT/TZ) ECP{def2-TZVP} COSMO Grid3 FinalGrid5 SlowConv LooseOpt printbasis', previous='/fs/home/bas348/perovskites/SMRFF/orca/%s/gbw_files/%s.orca.out'%(old_name, old_name),queue='long',mem=5000,procs=4,priority=250, extra_section=extra_section)



def restart_on_dh():
	for i in range(0,1):
		for solvent,dielectric in zip(['THTO'],[42.84]):
			for solute in ['pb2+']:
				extra_section='%%cosmo  SMD true  epsilon %f  end' % dielectric
				old_name = '%s_%s_tz_%d' % (solute,solvent,i)
				new_name = '%s_%s_dh_%d' % (solute,solvent,i)
				if os.path.exists('orca/%s/%s.out'%(old_name, old_name)):
					orca.job(new_name, '! RIJCOSX RI-PWPB95 D3BJ Def2-QZVPP(-g,-f) ECP{def2-TZVP} COSMO TIGHTSCF Grid5 FinalGrid6 SlowConv printbasis', previous='/fs/home/bas348/perovskites/orca/%s/%s.orca.gbw'%(old_name, old_name),queue='long',mem=5000, procs=12, extra_section=extra_section)
				elif os.path.exists('orca/%s/gbw_files/%s.orca.out'%(old_name, old_name)):
					orca.job(new_name, '! RIJCOSX RI-PWPB95 D3BJ Def2-QZVPP(-g,-f) ECP{def2-TZVP} COSMO TIGHTSCF Grid5 FinalGrid6 SlowConv printbasis', previous='/fs/home/bas348/perovskites/orca/%s/gbw_files/%s.orca.out'%(old_name, old_name),queue='long',mem=5000, procs=12, extra_section=extra_section)



def TZ():
	old_names=['pb2+_THTO_0']#,'pb2+_THTO_3','pb2+_THTO_4','pb2+_THTO_9'
	dielectrics=[42.84]
	for old_name,dielectric in zip(old_names, dielectrics):
		extra_section='%%cosmo  SMD true  epsilon %f  end' % dielectric
		for name in old_name:
			new_name=old_name[:-1]+'tz_0'		
			orca.job(new_name, '! Opt B97-D3 def2-TZVP GCP(DFT/TZ) ECP{def2-TZVP} COSMO Grid3 FinalGrid5 SlowConv LooseOpt printbasis', previous='/fs/home/bas348/perovskites/SMRFF/orca/%s/%s.orca.gbw'%(old_name, old_name),queue='long',mem=3000,procs=8, extra_section=extra_section)
  

def restart_on_TZ_wrong_charge():
	for i in range(20):
		for solvent,dielectric in zip(['benzene','odcb','nitromethane'],[2.27,9.93,35.87]):
			for solute in ['None', 'pb2+']:
				extra_section='%%cosmo  SMD true  epsilon %f  end' % dielectric
				old_name = '%s_%s_%d' % (solute,solvent,i)
				new_name = '%s_%s_tz_%d' % (solute,solvent,i)
				if os.path.exists('SMRFF/orca/%s/%s.out'%(old_name, old_name)):
					orca.job(new_name, '! Opt B97-D3 def2-TZVP GCP(DFT/TZ) ECP{def2-TZVP} COSMO printbasis', previous='/fs/home/bas348/perovskites/SMRFF/orca/%s/%s.orca.gbw'%(old_name, old_name),queue='medium', extra_section=extra_section)

def restart_on_TZ_right_charge():
	for i in range(1):
		for solvent,dielectric in zip(['benzene','odcb','nitromethane'],[2.27,9.93,35.87]):
			for solute in ['None', 'pb2+']:
				extra_section='%%cosmo  SMD true  epsilon %f  end' % dielectric
				old_name = '%s_%s_%d' % (solute,solvent,i)
				new_name = '%s_%s_tz.1_%d' % (solute,solvent,i)
				if os.path.exists('SMRFF/orca/%s/%s.out'%(old_name, old_name)):
					orca.job(new_name, '! Opt LooseOpt B97-D3 GCP(DFT/TZ) def2-TZVP ECP{def2-TZVP} Grid7  COSMO printbasis', previous='/fs/home/bas348/perovskites/SMRFF/orca/%s/%s.orca.gbw'%(old_name, old_name), queue='medium', extra_section=extra_section, charge=2 if solute=='pb2+' else 0)




def restart_double_hybrid():

	for i in range(20):
		#for solvent,dielectric in zip(['acetone','ACN','DMF','dmso', 'gbl', 'methacrolein', 'nmp'], [20.7, 37.5, 36.7, 46.7, 40.24, 10.9 ,32.2],):
		for solvent,dielectric in zip(['THTO'],[42.84]):
			for solute in ['None', 'pb2+']:
				extra_section='%%cosmo  SMD true  epsilon %f  end' % dielectric # originals converged too soon, showing no concern for gradient convergence criteria! Replace with LooseOpt.
				old_name = '%s_%s_tz_%d' % (solute,solvent,i)
				name = '%s_%s_dh_%d' % (solute,solvent,i)

				if os.path.exists('orca/%s/%s.out'%(old_name, old_name)):
					if 'ABORTING THE RUN' not in open('orca/%s/%s.out'%(old_name, old_name)).read():
						print  old_name, 'finished successfully'
						continue

					if solute=='None':
						orca.job(name, '! RIJCOSX RI-PWPB95 D3BJ Def2-QZVPP(-g,-f) ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6 SlowConv', previous=old_name,procs=8, queue='long', mem=3000, extra_section=extra_section)# charge_and_multiplicity='0 1'
					else:
						orca.job(name, '! RIJCOSX RI-PWPB95 D3BJ Def2-QZVPP(-g,-f) ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6 SlowConv', previous=old_name,procs=8, queue='long', mem=3000, extra_section=extra_section)
				




def mem_test():

	extra_section='%%cosmo  SMD true  epsilon %f  end' % 46.7 # originals converged too soon, showing no concern for gradient convergence criteria! Replace with LooseOpt.
	old_name = 'pb2+_dmso.4_0'
	atoms = orca.read('pb2+_dmso.3_0').atoms
	for mem in [100,250,500,1000,1500,2000,3000,4000]:
		name = 'pb2+_dmso.4_0_mem_%d' % mem
		orca.job(name, '! RIJCOSX RI-PWPB95 D3BJ Def2-QZVPP(-g,-f) ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6 SlowConv', atoms=atoms, previous=old_name, queue='batch', mem=mem, extra_section=extra_section, charge_and_multiplicity='2 1')
				


def restart_double_hybrid_again():

	for solvent,dielectric in zip([ 'gbl', 'nmp'], [ 40.24,32.2],):
		for solute,charge in zip(['pb2+', 'None'], ['2 1', '0 1']):
			count = 0
			for i in range(50):
				extra_section='%%cosmo  SMD true  epsilon %f  end' % dielectric # originals converged too soon, showing no concern for gradient convergence criteria! Replace with LooseOpt.
				extra_section += '\n%mp2 Q1Opt -1 end\n%scf maxiter 0 end'
				old_name = '%s_%s.5_%d' % (solute,solvent,i)
				name = '%s_%s.6_%d' % (solute,solvent,i)

				if os.path.exists('orca/%s/%s.out'%(old_name, old_name)):
					contents = open('orca/%s/%s.out'%(old_name, old_name)).read()
					if 'ABORTING THE RUN' not in contents:
						print  old_name, 'finished successfully'
						count += 1
						continue
					elif 'An error has occured in the MP2 module' in contents:
						orca.job(name, '! RIJCOSX RI-PWPB95 D3BJ Def2-QZVPP(-g,-f) ECP{def2-TZVP}', previous=old_name, queue='long', mem=8000, extra_section=extra_section, charge_and_multiplicity=charge)
						continue

				old_name = '%s_%s.4_%d' % (solute,solvent,i)
				if os.path.exists('orca/%s/%s.out'%(old_name, old_name)):
					contents = open('orca/%s/%s.out'%(old_name, old_name)).read()
					if 'ABORTING THE RUN' not in contents:
						print  old_name, 'finished successfully'
						count += 1
						continue
					elif 'An error has occured in the MP2 module' in contents:
						orca.job(name, '! RIJCOSX RI-PWPB95 D3BJ Def2-QZVPP(-g,-f) ECP{def2-TZVP}', previous=old_name, queue='long', mem=8000, extra_section=extra_section, charge_and_multiplicity=charge)


			print solvent, solute, count


def restart_double_hybrid_nmp():

	for solvent,dielectric in zip(['nmp'], [32.2],):
		for solute,charge in zip(['pb2+'], ['2 1']):
			count = 0
			for i in range(50):
				if i in [19, 20, 9, 15, 13, 10, 12, 11, 3, 8, 5, 7, 2, 1, 21, 17, 16, 25, 18]:
					continue #these already worked
				extra_section='%%cosmo  SMD true  epsilon %f  end' % dielectric # originals converged too soon, showing no concern for gradient convergence criteria! Replace with LooseOpt.
				extra_section += '\n%mp2 Q1Opt -1 end\n%scf maxiter 0 end'
				old_name = '%s_%s.5_%d' % (solute,solvent,i)
				name = '%s_%s.6_%d' % (solute,solvent,i)
				queue = 'long'

				if os.path.exists('orca/%s/%s.out'%(old_name, old_name)):
					contents = open('orca/%s/%s.out'%(old_name, old_name)).read()
					if 'ABORTING THE RUN' not in contents:
						print  old_name, 'finished successfully'
						count += 1
						continue
					elif 'An error has occured in the MP2 module' in contents:
						orca.job(name, '! RIJCOSX RI-PWPB95 D3BJ Def2-QZVPP(-g,-f) ECP{def2-TZVP}', previous=old_name, queue=queue, mem=8000, extra_section=extra_section, charge_and_multiplicity=charge)
						continue

				old_name = '%s_%s.4_%d' % (solute,solvent,i)
				if os.path.exists('orca/%s/%s.out'%(old_name, old_name)):
					contents = open('orca/%s/%s.out'%(old_name, old_name)).read()
					if 'ABORTING THE RUN' not in contents:
						print  old_name, 'finished successfully'
						count += 1
						continue
					elif 'An error has occured in the MP2 module' in contents:
						orca.job(name, '! RIJCOSX RI-PWPB95 D3BJ Def2-QZVPP(-g,-f) ECP{def2-TZVP}', previous=old_name, queue=queue, mem=8000, extra_section=extra_section, charge_and_multiplicity=charge)


			print solvent, solute, count

def read_nmp():
	for i in range(50):
		old_name = 'pb2+_nmp.6_%d' % i
		if os.path.exists('orca/%s/%s.out'%(old_name, old_name)):
			contents = open('orca/%s/%s.out'%(old_name, old_name)).read()
			if 'Cannot open gbw file' in contents:
				print i, 'gbw error'
			else:
				print i



def restart_double_hybrid_None_nmp():

	for solvent,dielectric in zip(['nmp'], [32.2],):
		for solute,charge in zip(['None'], ['0 1']):
			count = 0
			for i in range(50):
				
				extra_section='%%cosmo  SMD true  epsilon %f  end' % dielectric # originals converged too soon, showing no concern for gradient convergence criteria! Replace with LooseOpt.
				extra_section += '\n%mp2 Q1Opt -1 end\n%scf maxiter 0 end'
				old_name = '%s_%s.5_%d' % (solute,solvent,i)
				name = '%s_%s.6_%d' % (solute,solvent,i)
				queue = 'long'

				if os.path.exists('orca/%s/%s.out'%(old_name, old_name)):
					contents = open('orca/%s/%s.out'%(old_name, old_name)).read()
					if 'ABORTING THE RUN' not in contents:
						print  old_name, 'finished successfully'
						count += 1
						continue
					elif 'An error has occured in the MP2 module' in contents:
						orca.job(name, '! RIJCOSX RI-PWPB95 D3BJ Def2-QZVPP(-g,-f) ECP{def2-TZVP}', previous=old_name, queue=queue, mem=8000, extra_section=extra_section, charge_and_multiplicity=charge)
						continue

				old_name = '%s_%s.4_%d' % (solute,solvent,i)
				if os.path.exists('orca/%s/%s.out'%(old_name, old_name)):
					contents = open('orca/%s/%s.out'%(old_name, old_name)).read()
					if 'ABORTING THE RUN' not in contents:
						print  old_name, 'finished successfully'
						count += 1
						continue
					elif 'An error has occured in the MP2 module' in contents:
						orca.job(name, '! RIJCOSX RI-PWPB95 D3BJ Def2-QZVPP(-g,-f) ECP{def2-TZVP}', previous=old_name, queue=queue, mem=8000, extra_section=extra_section, charge_and_multiplicity=charge)


			print solvent, solute, count


def restart_double_hybrid_nmp_TZVP():

	for run,basis in [(7,'SVP'), (8,'TZVP')]:
		#for digit,extra in enumerate(['%mp2 density relaxed end', '%scf maxiter 0 end', '']):
		for digit,extra in enumerate(['%mp2 density relaxed end']):
			for solvent,dielectric in zip(['nmp'], [32.2],):
				for solute,charge,i in zip(['pb2+', 'None'], ['2 1', '0 1'], [20,20]):
					count = 0

					extra_section='%%cosmo  SMD true  epsilon %f  end\n' % dielectric
					extra_section += extra
					old_name = '%s_%s.5_%d' % (solute,solvent,i)
					name = '%s_%s.%d.%d_%d' % (solute,solvent,run,digit,i)
					queue = 'batch'

					if os.path.exists('orca/%s/%s.out'%(old_name, old_name)):
						orca.job(name, '! RIJCOSX RI-PWPB95 D3BJ Def2-'+basis+' ECP{def2-TZVP}', previous=old_name, queue=queue, mem=4000, extra_section=extra_section, charge_and_multiplicity=charge)
						continue

					old_name = '%s_%s.4_%d' % (solute,solvent,i)
					if os.path.exists('orca/%s/%s.out'%(old_name, old_name)):
						orca.job(name, '! RIJCOSX RI-PWPB95 D3BJ Def2-'+basis+' ECP{def2-TZVP}', previous=old_name, queue=queue, mem=4000, extra_section=extra_section, charge_and_multiplicity=charge)

TZ_new_func()

