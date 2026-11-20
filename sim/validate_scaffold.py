#!/usr/bin/env python3
"""Validate a scaffold FASTA against a reference genome.

Aligns with minimap2 (asm5 preset) and checks:
  - aligned coverage of the reference
  - block collinearity along the scaffold (circular genomes get one
    wrap for free)
  - misjoins = blocks whose reference span OVERLAPS the traversed
    frontier by more than --overlap-slack (unaligned slivers and
    small gaps between blocks are normal, overlaps are not)

Usage:
  validate_scaffold.py REF_FASTA SCAFFOLD_FASTA [MINIMAP2_BIN]
                       [--overlap-slack 500] [--origin-slack 50000]
"""

import subprocess
import sys
import argparse
import os
import tempfile


def read_fasta_sizes(path):
    tot, n = 0, 0
    ln = 0
    for line in open(path):
        if line.startswith(">"):
            if ln:
                tot += ln
                n += 1
            ln = 0
        else:
            ln += len(line.strip())
    if ln:
        tot += ln
        n += 1
    return n, tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("scaffold")
    ap.add_argument("mm2", nargs="?", default="minimap2/minimap2/minimap2")
    ap.add_argument("--overlap-slack", type=int, default=500,
                    help="block-frontier overlap tolerated before "
                         "calling a misjoin [500]")
    ap.add_argument("--origin-slack", type=int, default=50000,
                    help="reference-position window that counts as "
                         "the circular origin [50000]")
    args = ap.parse_args()

    n_ref, ref_size = read_fasta_sizes(args.ref)
    n_scaf, scaf_size = read_fasta_sizes(args.scaffold)

    with tempfile.NamedTemporaryFile("w", delete=False,
                                     suffix=".paf") as tf:
        paf_path = tf.name
    try:
        with open(paf_path, "w") as out:
            subprocess.run([args.mm2, "-x", "asm5", "-c",
                            "--secondary=no", args.ref, args.scaffold],
                           stdout=out, stderr=subprocess.DEVNULL,
                           check=True)
        rows = []
        for line in open(paf_path):
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or f[5] == "*":
                continue
            rows.append((int(f[2]), int(f[3]), f[4],
                         int(f[7]), int(f[8]), int(f[9])))
    finally:
        os.unlink(paf_path)

    rows.sort()
    aligned = sum(r[5] for r in rows)

    # frontier walk: misjoins only when a block OVERLAPS the
    # traversed reference beyond slack
    prev_t = None
    wraps = 0
    misjoins = []
    for i, (qs, qe, strand, ts, te, mlen) in enumerate(rows):
        if prev_t is None:
            prev_t = te if strand == "+" else ts
            continue
        bad = False
        if strand == "+":
            if ts + args.overlap_slack < prev_t:
                if wraps == 0 and ts < args.origin_slack:
                    wraps = 1          # the circular origin
                else:
                    bad = True
            prev_t = max(prev_t, te)
        else:
            if te > prev_t + args.overlap_slack:
                if wraps == 0 and te > ref_size - args.origin_slack:
                    wraps = 1
                else:
                    bad = True
            prev_t = min(prev_t, ts)
        if bad:
            misjoins.append((qs, prev_t, ts if strand == "+" else te,
                             mlen))

    print(f"reference:            {ref_size:,} bp ({n_ref} seqs)")
    print(f"scaffold:             {scaf_size:,} bp ({n_scaf} seqs)")
    print(f"aligned blocks:       {len(rows)}")
    print(f"aligned bp:           {aligned:,} "
          f"({100.0 * aligned / ref_size:.1f}% of reference)")
    print(f"circular wraps:        {wraps}")
    print(f"misjoins (overlaps):  {len(misjoins)}")
    for qs, front, pos, mlen in misjoins:
        print(f"    at scaffold {qs:,}: frontier {front:,} "
              f"vs block at {pos:,} ({mlen:,} matches)")

    return 0 if not misjoins else 1


if __name__ == "__main__":
    sys.exit(main())
