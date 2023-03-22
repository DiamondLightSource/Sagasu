import re
import os

def replace(file, pattern, subst):
    file_handle = open(file, "r")
    file_string = file_handle.read()
    file_handle.close()
    file_string = re.sub(pattern, subst, file_string)
    file_handle = open(file, "w")
    file_handle.write(file_string)
    file_handle.close()
    
res = 2.5
file = "prasa.txt"
resfile = os.path.join(os.getcwd(), "prasaout.txt")
with open(file, "r") as infile, open(resfile, "w") as outfile:
    datalines = []
    for line in infile:
        if line.startswith("End of trial"):
            if not line.endswith(" 0\n"):
                print(line)
                outfile.write(line)
        else:
            pass
        
replace(resfile, "End of trial ", "")
replace(resfile, "final CC is ", "")
replace(resfile, "CCrange is ", "")
replace(resfile, "CCall is ", "")
replace(resfile, "(candidate for a solution)", "")

with open(resfile, "r") as file, open("strip.txt", "w") as out:
    for line in file:
        lineout = line.replace(" ()", "")
        lineout = str(res) + ", " + lineout
        out.write(lineout)