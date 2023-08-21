import subprocess
import re
import os
from multiprocessing import Pool
from itertools import combinations
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import plotly.express as px



def run_command(file1, file2):
    command = f"module load phenix && phenix.emma --symmetry=scaled.mtz {str(file1)} {str(file2)}"  # Separate command and arguments
    result = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
    return (file1, file2, result.stdout.decode('utf-8'))

if __name__ == "__main__":

    pool = Pool(os.cpu_count() - 1)

    pairs_pattern = r"Pairs:\s*(\d+)"
    singles_model1_pattern = r"Singles model 1:\s*(\d+)"
    singles_model2_pattern = r"Singles model 2:\s*(\d+)"
    projname = "testing"
    highres = 22
    lowres = 26
    highsites = 29
    lowsites = 23

    input_files = [os.path.join(os.getcwd(), f"{projname}/{str(i)}/{str(j)}/{projname}_fa.pdb") for i in range(highres, lowres) for j in range(lowsites, highsites)]
    results = []
    parallel_filelist = list(combinations(input_files, 2))

    results = pool.starmap(run_command, parallel_filelist)

    percentages = []
    
    for file in input_files:
        filename = str(os.path.basename(os.path.dirname(os.path.dirname(file)))) + "_" + str(os.path.basename(os.path.dirname(file)))
        diagonalval = ([str(filename), str(filename), str(1)])
        percentages.append(diagonalval)

    for file1, file2, output in results:
        pairs_match = re.search(pairs_pattern, output)
        singles_model1_match = re.search(singles_model1_pattern, output)
        singles_model2_match = re.search(singles_model2_pattern, output)

        if pairs_match and singles_model1_match and singles_model2_match:
            pairs_number = pairs_match.group(1)
            singles_model1_number = singles_model1_match.group(1)
            singles_model2_number = singles_model2_match.group(1)
        else:
            print("Cannot match text")
        
        if pairs_number and singles_model1_number and singles_model2_number:
            percentage = np.around((float(pairs_number) / (float(pairs_number) + float(singles_model1_number) + float(singles_model2_number))), 2)
            filename1 = str(os.path.basename(os.path.dirname(os.path.dirname(file1)))) + "_" + str(os.path.basename(os.path.dirname(file1)))
            filename2 = str(os.path.basename(os.path.dirname(os.path.dirname(file2)))) + "_" + str(os.path.basename(os.path.dirname(file2)))
            percentages.append([str(filename1), str(filename2), str(percentage)])
        else:
            print("Cannot get numbers out")

    df = pd.DataFrame(percentages, columns=["filename_1", "filename_2", "percentage"])
    print(df)
    df.to_csv("df.csv")
    matrix_df = df.pivot(index="filename_1", columns="filename_2", values="percentage")
    matrix_df.fillna(0, inplace=True)
    print(matrix_df)
    matrix_df.to_csv("matrixdf.csv")
    matrix_df.drop('filename_1', axis=1)
    matrix_df.iloc[:, 1:] = matrix_df.iloc[:, 1:] * 100
    matrix_df = matrix_df.corr()

    plt.figure(figsize=(10, 8))
    sns.heatmap(matrix_df, annot=True, cmap="YlGnBu")
    plt.title("Percentage Correlation Matrix")
    plt.xlabel("Filename 2")
    plt.ylabel("Filename 1")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig("heat.png")
    
    numeric_data = matrix_df.iloc[:, 1:].values