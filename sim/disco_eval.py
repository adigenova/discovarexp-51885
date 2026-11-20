#!/usr/bin/env python3
"""Disco-style evaluation of an assembly FASTA against a reference.

Reproduces the metrics of the DiscovarDenovo `a.final/stats` block
for any FASTA (the pipeline computes them internally from the graph;
this works from the sequence), plus consensus-quality details:

  - bases in 1 kb+ / 10 kb+ scaffolds (and % of genome)
  - scaffold N50
  - N50 PERFECT STRETCH (error-free run N50)  <- consensus quality
  - mismatch + indel rates per 100 kb over aligned bases
  - optional: per-ONT-fill identity (--fills FASTA), the quality of
    the read-consensus gap sequences

Usage:
  disco_eval.py REF_FASTA ASSEMBLY_FASTA [MINIMAP2_BIN]
                [--fills FILLS_FASTA] [--preset asm5]
"""

import subprocess
import sys
import argparse
import os
import tempfile


def fasta_sizes(path):
    out, name, ln = [], None, 0
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                out.append((name, ln))
            name, ln = line[1:].strip(), 0
        else:
            ln += len(line.strip())
    if name is not None:
        out.append((name, ln))
    return out


def n50(lens):
    if not lens:
        return 0
    total = sum(lens)
    acc = 0
    for l in sorted(lens, reverse=True):
        acc += l
        if acc >= total / 2:
            return l
    return 0


def align(mm2, preset, ref, query):
    """Return list of (qname, qsize, strand, tstart, tend, cs_ops)."""
    with tempfile.NamedTemporaryFile("w", delete=False,
                                     suffix=".paf") as tf:
        paf = tf.name
    try:
        with open(paf, "w") as out:
            subprocess.run([mm2, "-x", preset, "-c", "--cs",
                            "--secondary=no", ref, query],
                           stdout=out, stderr=subprocess.DEVNULL,
                           check=True)
        rows = []
        for line in open(paf):
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or f[5] == "*":
                continue
            cs = None
            for tag in f[12:]:
                if tag.startswith("cs:Z:"):
                    cs = tag[5:]
                    break
            if cs is None:
                continue
            rows.append((f[0], int(f[1]), f[4], int(f[7]),
                          int(f[8]), cs))
    finally:
        os.unlink(paf)
    return rows


def parse_cs(cs):
    """Yield (kind, length) ops: kind in {'=', '*', '+', '-'}."""
    i = 0
    while i < len(cs):
        c = cs[i]
        if c == ":":
            j = i + 1
            while j < len(cs) and cs[j].isdigit():
                j += 1
            yield ("=", int(cs[i + 1:j]))
            i = j
        elif c == "*":
            yield ("*", 1)
            i += 3
        elif c in "+-":
            j = i + 1
            while j < len(cs) and cs[j] in "ACGTNacgtn":
                j += 1
            yield (c, j - i - 1)
            i = j
        else:
            i += 1


def perfect_stretches(rows):
    """Error-free runs over aligned bases; returns (runs, matched,
    mismatches, ins, dels).  Runs are per-alignment (they cannot
    span scaffold boundaries)."""
    runs, matched, mm, ins, dele = [], 0, 0, 0, 0
    for _, _, _, _, _, cs in rows:
        cur = 0
        for kind, n in parse_cs(cs):
            if kind == "=":
                cur += n
                matched += n
            else:
                if cur > 0:
                    runs.append(cur)
                cur = 0
                if kind == "*":
                    mm += 1
                elif kind == "+":
                    ins += n
                else:
                    dele += n
        if cur > 0:
            runs.append(cur)
    return runs, matched, mm, ins, dele


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("assembly")
    ap.add_argument("mm2", nargs="?",
                    default="minimap2/minimap2/minimap2")
    ap.add_argument("--fills", default=None,
                    help="FASTA of ONT-consensus gap fills")
    ap.add_argument("--preset", default="asm5")
    args = ap.parse_args()

    _, ref_size = fasta_sizes(args.ref)[0]
    scafs = fasta_sizes(args.assembly)
    scaf_lens = [l for _, l in scafs]
    total1 = sum(l for l in scaf_lens if l >= 1000)
    total10 = sum(l for l in scaf_lens if l >= 10000)

    print("# assembly statistics (Disco style)")
    print(f"scaffolds:            {len(scafs)}")
    print(f"scaffold N50:         {n50(scaf_lens):,}")
    print(f"total bases in 1 kb+ scaffolds:  {total1:,} "
          f"({100.0*total1/ref_size:.1f}% of genome)")
    print(f"total bases in 10 kb+ scaffolds: {total10:,} "
          f"({100.0*total10/ref_size:.1f}% of genome)")

    rows = align(args.mm2, args.preset, args.ref, args.assembly)
    runs, matched, mm, ins, dele = perfect_stretches(rows)
    print(f"aligned bases:        {matched:,} "
          f"({100.0*matched/ref_size:.1f}% of genome)")
    print(f"mismatch rate:        {100000.0*mm/max(1,matched):.2f} "
          f"per 100 kb")
    print(f"indel rate:           "
          f"{100000.0*(ins+dele)/max(1,matched):.2f} per 100 kb")
    print(f"N50 perfect stretch:  {n50(runs):,}")

    if args.fills and os.path.isfile(args.fills):
        print()
        print("# ONT-consensus gap fills")
        fills = []
        name, buf = None, []
        for line in open(args.fills):
            if line.startswith(">"):
                if name:
                    fills.append((name, "".join(buf)))
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
        if name:
            fills.append((name, "".join(buf)))
        if not fills:
            print("  (none)")
        with tempfile.NamedTemporaryFile("w", delete=False,
                                         suffix=".fa") as tf:
            tmpfa = tf.name
        try:
            with open(tmpfa, "w") as out:
                for name, s in fills:
                    out.write(">" + name.split()[0] + "\n")
                    for i in range(0, len(s), 70):
                        out.write(s[i:i + 70] + "\n")
            frows = align(args.mm2, "map-ont", args.ref, tmpfa)
        finally:
            os.unlink(tmpfa)
        # per-fill identity
        byfill = {}
        for q, _, _, _, _, cs in frows:
            m = e = 0
            for kind, n in parse_cs(cs):
                if kind == "=":
                    m += n
                elif kind == "*":
                    e += 1
                elif kind in "+-":
                    e += n
            if q not in byfill or m > byfill[q][0]:
                byfill[q] = (m, e)
        stats = []
        for name, s in fills:
            q = name.split()[0]
            m, e = byfill.get(q, (0, len(s)))
            ident = 100.0 * m / max(1, m + e)
            stats.append((ident, len(s), q))
        stats.sort()
        idents = [s[0] for s in stats]
        tot = sum(s[1] for s in stats)
        print(f"  fills:              {len(fills)} "
              f"({tot:,} bp)")
        print(f"  identity vs ref:    mean "
              f"{sum(idents)/len(idents):.2f}%  worst "
              f"{min(idents):.2f}%  best {max(idents):.2f}%")
        print("  worst five:")
        for ident, ln, q in stats[:5]:
            print(f"    {q}: {ident:.2f}% ({ln} bp)")


if __name__ == "__main__":
    sys.exit(main())
