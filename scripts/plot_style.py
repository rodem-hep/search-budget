#!/usr/bin/env python3
"""Shared figure palette for the statistics plots (see docs/METHOD_NOTES.md)."""

SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e0e0dc"
BLUE = {3: "#a5c9ec", 4: "#5f9ae0", 5: "#2a78d6", 6: "#1a4f8f"}
C_BKG, C_ARG, C_THR, C_ALT = "#4a3aa7", "#d03b3b", "#5f9e6e", "#e06a2a"

# single-series mark, and its pale partner for an "absent" cell
MARK, MARK_PALE = "#4c78a8", "#eef1f4"


def style(ax, grid=True):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("bottom", "left"):
        ax.spines[s].set_color(GRID)
        ax.spines[s].set_linewidth(1)
    ax.tick_params(colors=INK2, labelsize=10, length=3)
    ax.set_axisbelow(True)
    if grid:
        ax.grid(color=GRID, lw=0.8, alpha=0.6)
    return ax


def panels(ncol=2, figsize=(13.4, 5.1), wspace=0.22, **kw):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, ncol, figsize=figsize, facecolor=SURF,
                             gridspec_kw=dict(wspace=wspace, **kw))
    for ax in (axes if ncol > 1 else [axes]):
        style(ax)
    return fig, axes


def title(ax, text, size=13):
    ax.set_title(text, color=INK, fontsize=size, pad=12, loc="left")


def labels(ax, x=None, y=None):
    if x:
        ax.set_xlabel(x, color=INK2, fontsize=10.5)
    if y:
        ax.set_ylabel(y, color=INK2, fontsize=10.5)
