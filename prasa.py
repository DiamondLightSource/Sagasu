#!
import os


class prasa_run(object):
    def __init__(self):
        os.system("module load ccp4")

    def inputs(self):
        self.datain = str(input("Data file in (hkl/sca/mtz): "))
        self.atomin = str(input("Heavy atom: "))
        self.sitesin = str(input("Sites: "))
        self.trials = str(input("How many trials? "))

    def pointless(self):
        if self.datain.endswith(".hkl" or ".HKL"):
            self.pointless = (
                "pointless HKLOUT pointless.mtz HKLIN "
                + str(self.datain)
                + " > /dev/null 2>&1"
            )
        elif self.datain.endswith(".sca" or ".SCA"):
            self.pointless = "pointless HKLOUT pointless.mtz SCAIN " + str(
                self.datain + " > /dev/null 2>&1"
            )
        else:
            print("\nNo data given? Guessing...\n")
            self.pointless = (
                f"pointless HKLOUT pointless.mtz HKLIN "
                + str(self.datain)
                + " > /dev/null 2>&1"
            )
        os.system(str(self.pointless))

    def aimless(self):
        os.system(
            """
aimless HKLIN pointless.mtz HKLOUT aimless.mtz > /dev/null 2>&1 << eof
ANOMALOUS ON
eof
                  """
        )

    def ctruncate(self):
        os.system(
            "ctruncate -hklin aimless.mtz -hklout truncate.mtz -colin '/*/*/[I(+),SIGI(+),I(-),SIGI(-)]' > /dev/null 2>&1"
        )

    def afro(self):
        os.system(
            f"""
afro HKLIN truncate.mtz HKLOUT afro.mtz << eof
TARG SAD
Xtal S-SAD
ATOM {self.atomin}
NUMB {self.sitesin}
DNAME PEAK
COLUmn F+=F(+) SF+=SIGF(+) F-=F(-) SF-=SIGF(-)
EXCLUDE SIGF 3
eof
                """
        )

    def prasa(self):
        os.system(
            f"prasa -mtzin afro.mtz -colin-fo */*/[EA,SIGEA] -atom {self.atomin} -natoms {self.sitesin} -chargeflip 2 -ncycles 200 -pdbout prasa.pdb -ntrials {self.trials} -maxstatrescut 3.8 -numthreads 31"
        )


prasa = prasa_run()
prasa.inputs()
prasa.pointless()
prasa.aimless()
prasa.ctruncate()
prasa.afro()
prasa.prasa()


if __name__ == "__main__":
    prasa = prasa_run()
    prasa.inputs()
    prasa.pointless()
    prasa.aimless()
    prasa.ctruncate()
    prasa.afro()
    prasa.prasa()
