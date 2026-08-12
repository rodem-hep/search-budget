#!/usr/bin/env python3
"""Shared figure style for the statistics plots (see docs/METHOD_NOTES.md).

Conventions follow those of a published collider-physics figure rather than a report
graphic, because that is where these end up: Computer Modern to match the LaTeX body
text, a closed frame with major and minor ticks turned inward on all four sides, no
grid, and no in-axes title -- the LaTeX caption carries the description.

Figures are drawn at their *final printed size*. TEXTW is the \\textwidth of both
documents in inches, so a figure included at 0.9\\textwidth should be 0.9*TEXTW wide
and its labels come out at the point size set here. Drawing oversized and letting
\\includegraphics shrink it is what makes axis text illegible.

Colours: the categorical slots are Okabe--Ito, a published palette designed to stay
separable under the common colour-vision deficiencies; the BLUE ramp is single-hue
light-to-dark for magnitude. Keep them in those roles -- a sequential ramp used for
identity, or a categorical hue used for magnitude, both misread.
"""
import re
import matplotlib as mpl
import matplotlib.ticker as mticker

TEXTW = 6.2          # \textwidth of a4paper with the documents' margins, in inches

SURF, INK, INK2, GRID = "#ffffff", "#000000", "#1a1a1a", "#d9d9d9"

# sequential: one hue, light -> dark (ColorBrewer Blues)
BLUE = {3: "#c6dbef", 4: "#6baed6", 5: "#2171b5", 6: "#08306b"}

# categorical: Okabe-Ito. blue / vermillion / bluish-green / orange
C_BKG, C_ARG, C_THR, C_ALT = "#0072b2", "#d55e00", "#009e73", "#e69f00"

# single-series mark, and its pale partner for an "absent" cell
MARK, MARK_PALE = "#0072b2", "#e8eef4"

# The axes description belongs in the LaTeX caption, so title() draws nothing.
# Set True when running a script standalone and the figure has to speak for itself.
SHOW_TITLES = False

mpl.rcParams.update({
    "figure.facecolor": SURF,
    "axes.facecolor": SURF,
    "savefig.facecolor": SURF,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "savefig.dpi": 400,
    # Computer Modern, to match the body text of both documents
    "font.family": "serif",
    "font.serif": ["cmr10", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,          # cmr10 has no U+2212
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
    """figsize for a figure included at `frac`*\\textwidth, height/width = aspect."""
    w = frac * TEXTW
    return (w, w * aspect)


def style(ax, grid=False, minor="both"):
    """Close the frame, turn the ticks inward, drop the grid.

    minor: which axes carry minor ticks -- "both", "x", "y" or None. Categorical
    axes want them off, since there is nothing between two adjacent categories.
    """
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


# ---------------------------------------------------------------- observable labels
# The observable keys are ASCII ("m(gammagamma)") so that tables stay terminal-readable.
# On a figure they should look like the physics: m(gamma gamma). Longest token first,
# so that mumu is not eaten by mu.
_PARTS = [("multi", r"\mathrm{multi}"),          # before "mu", which it contains
          ("gammagamma", r"\gamma\gamma"), ("gamma", r"\gamma"),
          ("mutau", r"\mu\tau"), ("etau", r"e\tau"), ("emu", r"e\mu"),
          ("tautau", r"\tau\tau"), ("taunu", r"\tau\nu"), ("mumu", r"\mu\mu"),
          ("mub", r"\mu b"), ("muj", r"\mu j"), ("muZ", r"\mu Z"), ("muv", r"\mu\nu"),
          ("mu", r"\mu"), ("tau", r"\tau"), ("nu", r"\nu"), ("ev", r"e\nu")]
_MAP = dict(_PARTS)
# One pass, alternatives longest-first: a second pass would re-read the backslash
# it just wrote and turn \gamma into \\gamma.
_TOK = re.compile("|".join(re.escape(k) for k, _ in _PARTS))
_SAFE = re.compile(r"^[A-Za-z0-9\\^_{}+\-/() ]*$")


def textsafe(s):
    """cmr10 is OT1-encoded, where '<' and '>' are the inverted exclamation and question
    marks -- 'h->aa' comes out as 'h-!aa'. Push them through mathtext, which uses the CM
    math fonts and has the real glyphs."""
    s = s.replace("->", r"$\to$")
    return s.replace("<", r"$<$").replace(">", r"$>$")


def _greek(s):
    # trailing space closes the command name, so \mu + jj is not read as \mujj;
    # math mode ignores the space itself
    out = _TOK.sub(lambda m: _MAP[m.group(0)] + " ", s)
    return re.sub(r" +", " ", out).strip()


def mathify(lab):
    """'m(gammagamma)' -> '$m(\\gamma\\gamma)$', 'm(ee) SS' -> '$m(ee)$ SS'.

    Anything that does not parse is returned untouched: a wrong label is worse than
    an ASCII one.
    """
    m = re.match(r"^(m|mT|m_T)\((.*?)\)(.*)$", lab)
    if not m:
        return lab
    head, inner, tail = m.group(1), m.group(2), m.group(3)
    head = r"m_{\mathrm{T}}" if head in ("mT", "m_T") else "m"
    inner = _greek(inner)
    # a second m(...) in the tail, as in m(tt)/m(jj)
    tail_m = re.match(r"^/(m|mT)\((.*?)\)(.*)$", tail)
    if tail_m:
        inner2 = _greek(tail_m.group(2))
        h2 = r"m_{\mathrm{T}}" if tail_m.group(1) == "mT" else "m"
        body, tail = f"{head}({inner})/{h2}({inner2})", tail_m.group(3)
    else:
        body = f"{head}({inner})"
    if not _SAFE.match(body):
        return lab
    return f"${body}$" + (f" {tail.strip()}" if tail.strip() else "")
