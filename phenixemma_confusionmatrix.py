import subprocess
import re
import os
from multiprocessing import Pool
from itertools import combinations

# emma works with paths
# can use the scaled.mtz or similar input file to extract symmetry information eg.
# phenix.emma --symmetry=scaled.mtz --diffraction_index_equivalent 19_25.pdb 19_26.pdb

pairs_pattern = r"Pairs:\s*(\d+)"
singles_model1_pattern = r"Singles model 1:\s*(\d+)"
singles_model2_pattern = r"Singles model 2:\s*(\d+)"
projname = "testing"
highres = 19
lowres = 27
highsites = 30
lowsites = 20

os.system("module load phenix")

input_files = [os.path.join(os.getcwd(), f"{projname}/{str(i)}/{str(j)}/{projname}_fa.pdb") for i in range(highres, lowres) for j in range(lowsites, highsites)]
results = []

def run_command(file1, file2):
    command = ["phenix.emma", "--symmetry=scaled.mtz", file1, file2]
    result = subprocess.run(command, stdout=subprocess.PIPE)
    return (file1, file2, result.stdout.decode('utf-8'))

def collect_result(result):
    global results
    results.append(result)

parallel_filelist = list(combinations(input_files, 2))

with Pool() as pool:
    for file1, file2 in parallel_filelist:
        pool.apply_async(run_command, args=(file1, file2), callback=collect_result)
    pool.close()
    pool.join()

for file1, file2, output in results:
    pairs_match = re.search(pairs_pattern, output)
    singles_model1_match = re.search(singles_model1_pattern, output)
    singles_model2_match = re.search(singles_model2_pattern, output)

    if pairs_match:
        pairs_number = pairs_match.group(1)
        print("Pairs:", pairs_number)

    if singles_model1_match:
        singles_model1_number = singles_model1_match.group(1)
        print("Singles model 1:", singles_model1_number)

    if singles_model2_match:
        singles_model2_number = singles_model2_match.group(1)
        print("Singles model 2:", singles_model2_number)