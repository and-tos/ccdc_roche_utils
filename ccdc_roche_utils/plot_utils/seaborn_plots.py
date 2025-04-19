import seaborn as sns
from sklearn.metrics import mean_absolute_error
from scipy.stats import pearsonr
import numpy as np
from matplotlib.font_manager import fontManager, FontProperties
from matplotlib import pyplot as plt
from ccdc_roche_utils.plot_utils import roche_branding

path = "/Roche_fonts/Roche Sans/RocheSansLight-Medium.ttf"
fontManager.addfont(path)
prop = FontProperties(fname=path)
sns.set(font=prop.get_name())


def plot_unity(xdata, ydata, **kwargs):
    mn = min(xdata.min(), ydata.min())
    mx = max(xdata.max(), ydata.max())
    points = np.linspace(mn, mx, 100)
    plt.gca().plot(
        points, points, color="k", marker=None, linestyle="--", linewidth=1.0
    )


def regression_plot(
    data,
    x,
    y,
    outname=None,
    mae=True,
    hue=None,
    palette=[
        roche_branding.RocheColours().roche_colours_rgb_dict["roche_blue"],
        roche_branding.RocheColours().roche_colours_rgb_dict["dark_orange"],
    ],
    **kwargs,
):
    pearson_r2 = pearsonr(data[x], data[y])[0] ** 2
    line_kws = {"label": f"R2={pearson_r2:.3}   N={data.shape[0]}"}
    scatter_kws = {}
    print(pearson_r2)
    if mae:
        mae = mean_absolute_error(data[x], data[y])
        label = f"R2={pearson_r2:.3}   MAE={mae:.3}   N={data.shape[0]}"
        line_kws = {"label": label}

    if "line_kws" in kwargs:
        line_kws.update(kwargs.pop("line_kws"))
    if "scatter_kws" in kwargs:
        scatter_kws = kwargs.pop("scatter_kws")

    ax = sns.scatterplot(
        data,
        x=x,
        y=y,
        hue=hue,
        palette=palette,
        edgecolor="black",
        **scatter_kws,
        **kwargs,
    )
    ax = sns.regplot(
        data,
        x=x,
        y=y,
        scatter=False,
        ax=ax,
        label=label,
        line_kws=line_kws,
        scatter_kws=scatter_kws,
        **kwargs,
    )
    if not outname:
        return ax
    else:
        ax.legend()
        ax.figure.savefig(outname, dpi=600)
        ax.figure.clf()
        return


def heatmap_plot():
    return


def main():
    import pandas as pd

    df = pd.read_csv("ml_predictions.csv")
    for c in ["py_pred", "py_pre_pred", "sklearn_pred"]:
        regression_plot(df, x="delta_min_min [kcal/mol]", y=c, outname=f"{c}.png")


if __name__ == "__main__":
    main()
