import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import plotly.express as px


def ccalloutliers(filename, resolution, sitessearched):
    df = pd.read_csv(
        filename,
        sep=",",
        names=["linebeg", "TRY", "CPUNO", "CCALL", "CCWEAK", "CFOM", "BEST", "PATFOM",],
    )
    pd.DataFrame.drop(
        df,
        labels=["linebeg", "CPUNO", "BEST", "TRY", "CFOM", "PATFOM"],
        axis=1,
        inplace=True,
    )
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
    # df["WEIGHTED"] = df["NORM_COMB_VEC"] / df["NORM_VEC_DIFF"]
    df["WEIGHTED"] = np.power(df["NORM_COMB_VEC"], 18) / np.cbrt(df["NORM_VEC_DIFF"])
    # plt.scatter(df["CCWEAK"], df["CCALL"], c=df["WEIGHTED"], cmap="Blues", marker="o")
    # plt.show()
    df = df[df["WEIGHTED"] > 0.1]
    df["RES"] = resolution / 10
    df["SITES"] = sitessearched
    # print(df)
    return df


if __name__ == "__main__":
    all_data = pd.DataFrame()
    for x in range(30, 46):
        for y in range(3, 10):
            name = str(x) + "_" + str(y) + ".csv"
            data = ccalloutliers(name, x, y)
            all_data = pd.concat([all_data, data], axis=0, ignore_index=True)
            print(f"Res - {x / 10}   Sites - {y}")
    # plt.scatter(all_data["CCWEAK"], all_data["CCALL"], c=all_data["WEIGHTED"], cmap="Blues", marker="o")
    all_data.sort_values(by=["COMB_VEC"], axis=0, inplace=True, ascending=False)
    print(all_data)
    customdata = np.stack(
        (all_data["RES"], all_data["SITES"], all_data["COMB_VEC"]), axis=1
    )
    fig = px.scatter(
        all_data,
        x="CCWEAK",
        y="CCALL",
        color="COMB_VEC",
        color_continuous_scale="Bluered_r",
    )
    hovertemplate = (
        "Res: %{customdata[0]} Å<br>"
        + "Sites: %{customdata[1]}<br>"
        + "Distance: %{customdata[2]:,.3f}<br>"
        + "CCWeak: %{x} <br>"
        + "CCAll: %{y}"
        + "<extra></extra>"
    )
    fig.update_traces(customdata=customdata, hovertemplate=hovertemplate)
    fig.write_html("output.html")
    fig.show()
    # plt.colorbar()
    # plt.savefig("out2.png")
