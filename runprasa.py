#!/dls/science/groups/i23/pyenvs/ctrl_conda python3
import os
import sys
import numpy as np

atomin = sys.argv[1]
sitesin = sys.argv[2]
trialsin = sys.argv[3]
lowres = sys.argv[4]
highres = sys.argv[5]
minpeaks = str(np.round((int(sitesin) * 0.5), 0))
natoms = str(int(sitesin) * 3)

def afro():
    os.system(f'''
afro HKLIN truncate.mtz HKLOUT afro.mtz << eof
TARG SAD
Xtal S-SAD
ATOM {atomin}
NUMB {sitesin}
DNAME PEAK
COLUmn F+=F(+) SF+=SIGF(+) F-=F(-) SF-=SIGF(-)
EXCLUDE SIGF 3
eof
            ''')
    
def prasa():
    os.system(f"prasa -mtzin afro.mtz -colin-fo */*/[EA,SIGEA] -atom {atomin} -natoms {sitesin} -chargeflip 2 -ncycles 200 -pdbout prasa.pdb -natoms {natoms} -minpeaks {minpeaks} -ntrials {trialsin} -maxstatrescut {highres} -minstatrescut {lowres} -numthreads {str(os.cpu_count)} > prasaout.txt")
    
afro()
prasa()
