import re

_PARTS = [("multi", r"\mathrm{multi}"),
          ("gammagamma", r"\gamma\gamma"), ("gamma", r"\gamma"),
          ("mutau", r"\mu\tau"), ("etau", r"e\tau"), ("emu", r"e\mu"),
          ("tautau", r"\tau\tau"), ("taunu", r"\tau\nu"), ("mumu", r"\mu\mu"),
          ("mub", r"\mu b"), ("muj", r"\mu j"), ("muZ", r"\mu Z"), ("muv", r"\mu\nu"),
          ("mu", r"\mu"), ("tau", r"\tau"), ("nu", r"\nu"), ("ev", r"e\nu")]
_MAP = dict(_PARTS)
_TOK = re.compile("|".join(re.escape(k) for k, _ in _PARTS))
_SAFE = re.compile(r"^[A-Za-z0-9\\^_{}+\-/() ]*$")


def textsafe(s):
    s = s.replace("->", r"$\to$")
    return s.replace("<", r"$<$").replace(">", r"$>$")


def _greek(s):
    out = _TOK.sub(lambda m: _MAP[m.group(0)] + " ", s)
    return re.sub(r" +", " ", out).strip()


def mathify(lab):
    m = re.match(r"^(m|mT|m_T)\((.*?)\)(.*)$", lab)
    if not m:
        return lab
    head, inner, tail = m.group(1), m.group(2), m.group(3)
    head = r"m_{\mathrm{T}}" if head in ("mT", "m_T") else "m"
    inner = _greek(inner)
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
