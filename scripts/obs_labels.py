#!/usr/bin/env python3
"""Rendering the ASCII observable keys as physics.

The keys are ASCII ("m(gammagamma)") so that tables stay terminal-readable; a figure or a LaTeX
table should show m(gamma gamma). Both consumers need the same map, so it lives here rather than
in plot_style: matplotlib's mathtext and LaTeX accept the same $...$ output, and this module is
pure standard library, so the document generators keep working without numpy installed.

Imported by plot_style.py (figures) and export_model_map_tex.py (the appendix table).
"""
import re

# Longest token first, so that mumu is not eaten by mu.
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
    marks -- 'h->aa' comes out as 'h-!aa'. Push them through math mode, which uses the CM
    math fonts and has the real glyphs. OT1 LaTeX inverts them the same way."""
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
