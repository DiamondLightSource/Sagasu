afro HKLIN truncate.mtz HKLOUT afro.mtz << eof
TARG SAD
Xtal S-SAD
ATOM "$1"
NUMB "$2"
DNAMe PEAK
GIAC
COLUmn F+=F(+) SF+=SIGF(+) F-=F(-) SF-=SIGF(-)
EXCLuse SIGF 3 SANO 1
eof

sleep 5

prasa -mtzin afro.mtz -colin-fo */*/[EA,SIGEA] -chargeflip 2 -ncycles 200 -shannon 1.15 -numthreads 40 -pdbout prasa.pdb -histmatch 1 -atom "$1" -rescut "$3" -ntrials "$4" -maxstatrescut "$5" -minstatrescut "$6" -natoms "$7" -minpeaks "$8" > prasa.txt