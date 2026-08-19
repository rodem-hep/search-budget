import matplotlib as mpl
import matplotlib.ticker as mticker

TEXTW = 6.2

SURF, INK, INK2, GRID = "#ffffff", "#000000", "#1a1a1a", "#d9d9d9"

BLUE = {3: "#c6dbef", 4: "#6baed6", 5: "#2171b5", 6: "#08306b"}

C_BKG, C_ARG, C_THR, C_ALT = "#0072b2", "#d55e00", "#009e73", "#e69f00"

MARK, MARK_PALE = "#0072b2", "#e8eef4"

SHOW_TITLES = False

mpl.rcParams.update({
    "figure.facecolor": SURF,
    "axes.facecolor": SURF,
    "savefig.facecolor": SURF,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "savefig.dpi": 400,
    "font.family": "serif",
    "font.serif": ["cmr10", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
    "axes.formatter.use_mathtext": True,
    "font.size": 9,
    "axes.labelsize": 9.5,
    "axes.titlesize": 9.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8.5,
    "legend.frameon": False,
    "legend.handlelength": 1.9,
    "legend.labelspacing": 0.3,
    "legend.borderaxespad": 0.6,
    "lines.linewidth": 1.3,
    "lines.markersize": 4.0,
    "patch.linewidth": 0.7,
    "axes.linewidth": 0.7,
    "axes.edgecolor": INK,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": INK, "ytick.color": INK,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.major.size": 5.0, "ytick.major.size": 5.0,
    "xtick.minor.size": 2.6, "ytick.minor.size": 2.6,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.minor.width": 0.55, "ytick.minor.width": 0.55,
    "xtick.major.pad": 3.0, "ytick.major.pad": 3.0,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
})


def size(frac=1.0, aspect=0.62):
    w = frac * TEXTW
    return (w, w * aspect)


def style(ax, grid=False, minor="both"):
    ax.set_facecolor(SURF)
    for s in ("top", "right", "bottom", "left"):
        ax.spines[s].set_visible(True)
        ax.spines[s].set_color(INK)
        ax.spines[s].set_linewidth(0.7)
    ax.tick_params(which="both", direction="in", top=True, right=True, colors=INK)
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator()
                               if minor in ("both", "x") else mticker.NullLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator()
                               if minor in ("both", "y") else mticker.NullLocator())
    ax.set_axisbelow(True)
    if grid:
        ax.grid(color=GRID, lw=0.5, ls=":", alpha=0.9)
    return ax


def panels(ncol=2, figsize=None, wspace=0.24, **kw):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, ncol, figsize=figsize or size(1.0, 0.42),
                             facecolor=SURF, gridspec_kw=dict(wspace=wspace, **kw))
    for ax in (axes if ncol > 1 else [axes]):
        style(ax)
    return fig, axes


def title(ax, text, size=9.5):
    if SHOW_TITLES:
        ax.set_title(text, color=INK, fontsize=size, pad=6, loc="left")


def labels(ax, x=None, y=None):
    if x:
        ax.set_xlabel(x, color=INK)
    if y:
        ax.set_ylabel(y, color=INK)

