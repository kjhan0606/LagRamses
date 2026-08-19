#!/usr/bin/env python3
"""Phase-1 census v2: joins & continuations and accepts array-element child
operands, the two blind spots Fable's review exposed. Reports every remaining
hand-rolled cell-index site in the sources bin/Makefile actually compiles."""
import re, sys, pathlib

def logical_lines(text):
    """Join Fortran free-form '&' continuations, keeping the first line number."""
    raw = text.split('\n')
    out, buf, start = [], None, None
    for i, l in enumerate(raw, 1):
        s = l.rstrip()
        body = s.split('!', 1)[0] if not s.lstrip().lower().startswith('!$omp') else s
        cont = body.rstrip().endswith('&')
        piece = body.rstrip().rstrip('&')
        if buf is None:
            start, buf = i, piece
        else:
            buf += ' ' + piece.lstrip().lstrip('&')
        if not cont:
            out.append((start, buf, l)); buf = None
    if buf is not None: out.append((start, buf, ''))
    return out

ATOM = r'[A-Za-z_]\w*(?:\((?:[^()]|\([^()]*\))*\))?'
PATS = [
 ('stride',      re.compile(r'ncoarse\s*\+\s*\(\s*'+ATOM+r'\s*-\s*1\s*\)\s*\*\s*ngridmax\b', re.I)),
 ('stride_swap', re.compile(r'ncoarse\s*\+\s*ngridmax\b\s*\*\s*\(\s*'+ATOM+r'\s*-\s*1\s*\)', re.I)),
 ('rev_child',   re.compile(r'\(\s*'+ATOM+r'\s*-\s*ncoarse\s*(?:-\s*1\s*)?\)\s*/\s*ngridmax\b', re.I)),
 ('rev_grid',    re.compile(r'-\s*ncoarse\s*-\s*\(\s*'+ATOM+r'\s*-\s*1\s*\)\s*\*\s*ngridmax\b', re.I)),
 ('rev_grid_sw', re.compile(r'-\s*ncoarse\s*-\s*ngridmax\b\s*\*\s*\(\s*'+ATOM+r'\s*-\s*1\s*\)', re.I)),
 ('stride_0based', re.compile(r'ncoarse\s*\+\s*'+ATOM+r'\s*\*\s*ngridmax\b\s*\+', re.I)),
 ('rev_child_0b', re.compile(r'\(\s*'+ATOM+r'\s*-\s*ncoarse\s*\)\s*/\s*ngridmax\b', re.I)),
 ('stride_gridfirst', re.compile(ATOM+r'\s*\+\s*ncoarse\s*\+\s*'+ATOM+r'\s*\*\s*ngridmax\b', re.I)),
 ('mod_form',    re.compile(r'mod\s*\(\s*'+ATOM+r'\s*-\s*ncoarse[^)]*ngridmax\b', re.I)),
]
CAP = re.compile(r'twotondim\s*\*\s*ngridmax\b|ngridmax\b\s*\*\s*twotondim', re.I)

def sources():
    root = pathlib.Path('.')
    mk = (root/'bin/Makefile').read_text()
    # Makefile comments list alternative object sets (e.g. the .kisti variants).
    # Parsing them pulls in sources the build never compiles, which cannot be
    # gate-verified, so drop commented lines before extracting object names.
    mk = '\n'.join(l for l in mk.split('\n') if not l.lstrip().startswith('#'))
    vp = [p.replace('$(PATCH)','../patch/lagRamses').replace('../$(SOLVER)','../hydro')
          for p in re.search(r'^VPATH\s*=\s*(.+)$', mk, re.M).group(1).split(':')]
    vp = [(root/'bin'/p).resolve() for p in vp]
    objs = sorted(set(re.findall(r'([A-Za-z0-9_.]+)\.o\b', mk)))
    out = {}
    for o in objs:
        for d in vp:
            f = d/(o+'.f90')
            if f.exists(): out[o] = f; break
    return out

def main():
    only = set(sys.argv[1:])
    total = 0; per = {}
    for o, f in sorted(sources().items()):
        if only and str(f) not in only and o not in only: continue
        # amr_index.f90 IS the one legitimate home of the arithmetic
        if f.name == 'amr_index.f90': continue
        hits = []
        for ln, joined, raw in logical_lines(f.read_text(errors='replace')):
            s = raw.lstrip()
            if s.startswith('!') and not s.lower().startswith('!$omp'): continue
            if CAP.search(joined): continue
            for name, pat in PATS:
                if pat.search(joined):
                    hits.append((ln, name, joined.strip()[:100])); break
        if hits:
            per[str(f)] = hits; total += len(hits)
    for f, hits in sorted(per.items(), key=lambda kv: -len(kv[1])):
        print(f"\n{f}  ({len(hits)})")
        for ln, name, txt in hits:
            print(f"  {ln:>5} {name:<12} {txt}")
    print(f"\nTOTAL remaining sites: {total} in {len(per)} files")
    return 1 if total else 0

sys.exit(main())
