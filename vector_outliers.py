import pandas as pd
import numpy as np
import matplotlib.pylab as plt

def ccalloutliers(filename, resolution, sitessearched):
    df = pd.read_csv(
        filename,
        sep=",",
        names=[
            "linebeg",
            "TRY",
            "CPUNO",
            "CCALL",
            "CCWEAK",
            "CFOM",
            "BEST",
            "PATFOM",
        ],
    )
    pd.DataFrame.drop(df, labels=["linebeg", "CPUNO", "BEST", "TRY", "CFOM", "PATFOM"] , axis=1, inplace=True)
    ccall_mean = df["CCALL"].mean()
    ccweak_mean = df["CCWEAK"].mean()
    df["CCALL_VEC"] = df["CCALL"] - ccall_mean
    df["CCWEAK_VEC"] = df["CCWEAK"] - ccweak_mean
    df = df[df["CCALL_VEC"] > 0]
    df = df[df["CCWEAK_VEC"] > 0]
    df["VEC_DIFF"] = np.absolute(df["CCALL_VEC"] - df["CCWEAK_VEC"])
    df["COMB_VEC"] = np.sqrt(np.square(df["CCALL_VEC"]) + np.square(df["CCWEAK_VEC"]))
    df["NORM_VEC_DIFF"] = (df["VEC_DIFF"] / df["VEC_DIFF"].abs().max()) + 0.000001
    df["NORM_COMB_VEC"] = (df["COMB_VEC"] / df["COMB_VEC"].abs().max()) + 0.000001 
    #df["WEIGHTED"] = df["NORM_COMB_VEC"] / df["NORM_VEC_DIFF"]
    df["WEIGHTED"] = np.power(df["NORM_COMB_VEC"], 18) / np.cbrt(df["NORM_VEC_DIFF"])
    print(df)
    plt.scatter(df["CCWEAK"], df["CCALL"], c=df["WEIGHTED"], cmap="Blues", marker="o")
    plt.show()
    # allmad = open("vector_all.csv", "a")
    # allmad.write(
    #     str(int(resolution) / 10)
    #     + ","
    #     + str(sitessearched)
    #     + ","
    #     + str(mad5)
    #     + ","
    #     + str(mad6)
    #     + ","
    #     + str(mad7)
    #     + ","
    #     + str(mad8)
    #     + ","
    #     + str(mad9)
    #     + ","
    #     + str(mad10)
    #     + "\n"
    # )
    # allmad.close()

if __name__ == "__main__":
    for x in range(30,45):
        for y in range(20,65):
            name = str(x) + "_" + str(y) + ".csv"
            ccalloutliers(name, x, y)