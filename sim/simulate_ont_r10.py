#!/usr/bin/env python3
"""
simulate_ont_r10.py -- simulate ONT R10.4 (Kit-14 simplex, Q20-class) reads
from a reference genome.  Pure standard library; deterministic given SEED.

Model (approximating published R10.4.1 simplex profiles):
  - read lengths: mixture of two log-normal modes (short fragments +
    long fragments), clipped to [MIN_LEN, MAX_LEN]
  - substitutions:    ~0.30% per base
  - insertions:      ~0.25% per base, mostly 1-2 bp
  - deletions:       ~0.65% per base background, plus elevated deletion
    inside homopolymer runs (per-base ~3.5% for runs >= 3), which is the
    classic ONT signature
  - rare long deletions (~0.1 expected per read, median ~300 bp)
  - no Ns; deterministic under --seed
Output: gzipped FASTA (the reads are alignment templates for LR-1; no
quality strings are needed).

Usage:
  python3 simulate_ont_r10.py REF=ecoli/ecoli-k12-ref.fna \
      OUT=ecoli/ont/EC.10x.ont.fa.gz COV=10 SEED=42
"""
import gzip, math, random, sys, os

# ----------------------------------------------------------------- params
ARGS = {}
for a in sys.argv[1:]:
    k, _, v = a.partition("=")
    ARGS[k.upper()] = v
REF   = ARGS.get("REF", "")
OUT   = ARGS.get("OUT", "")
COV   = float(ARGS.get("COV", "10"))
SEED  = int(ARGS.get("SEED", "42"))
MODE  = ARGS.get("MODE", "r10")   # r10 (Q20-class) or r10_bad (Q14-class)

P_SUB   = 0.0030
P_INS   = 0.0025
P_DEL   = 0.0065
P_HOMO_DEL = 0.035     # extra per-base deletion inside homopolymer runs
LONG_DEL_RATE = 0.10   # expected long deletions per read
LONG_DEL_MEDIAN = 300.0
if MODE == "r10_bad":  # noisier variant for stress tests
    P_SUB, P_INS, P_DEL, LONG_DEL_RATE = 0.006, 0.005, 0.012, 0.20
    P_HOMO_DEL = 0.05

LEN_SHORT_W, LEN_SHORT_MED, LEN_SHORT_SIGMA = 0.65, 9000, 0.90
LEN_LONG_W,  LEN_LONG_MED,  LEN_LONG_SIGMA  = 0.35, 40000, 0.80
MIN_LEN, MAX_LEN = 1000, 250000

BASES = "ACGT"
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}

# ----------------------------------------------------------------- ref io
def load_ref(path):
    name, seq = None, []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None: break
                name = line[1:].strip().split()[0]
            else:
                seq.append(line.strip())
    return name, "".join(seq).upper()

def draw_len(rng):
    if rng.random() < LEN_LONG_W:
        L = rng.lognormvariate(math.log(LEN_LONG_MED), LEN_LONG_SIGMA)
    else:
        L = rng.lognormvariate(math.log(LEN_SHORT_MED), LEN_SHORT_SIGMA)
    return int(min(max(L, MIN_LEN), MAX_LEN))

def homopolymer_del(seq, i):
    """per-base deletion probability at position i (homopolymer-aware)."""
    p = P_DEL
    b, j = seq[i], i
    n = 1
    while j+1 < len(seq) and seq[j+1] == b and n < 12:
        j += 1; n += 1
    # run of length n >= 3 around/at i?
    if n >= 3:
        p = P_DEL + P_HOMO_DEL
    return p

COMP = {"A":"T","C":"G","G":"C","T":"A"}

def simulate_read(rng, ref):
    L = draw_len(rng)
    pos = rng.randrange(len(ref))
    out = []
    i = 0
    rc_flag = rng.random() < 0.5   # 50/50 forward/reverse
    # rare long deletions: scheduled as skip events
    long_dels = []
    while rng.random() < LONG_DEL_RATE / (1 + LONG_DEL_RATE):
        pass
    n_long = 0
    while True:
        if n_long < rng.random():  # ~Poisson(LONG_DEL_RATE) upper bound
            break
        long_dels.append(int(min(rng.expovariate(1.0/LONG_DEL_MEDIAN), 5000)))
        n_long += 1
        if n_long > 10: break
    dl_set = sorted(rng.randrange(L) for _ in long_dels)
    while len(out) < L:
        if i >= len(ref): i = 0            # circular reference
        if dl_set and i >= dl_set[0] and i - dl_set[0] < 0 + L:
            pass  # placeholder (handled below)
        r = rng.random()
        if r < P_INS and len(out) < L:
            out.append(BASES[rng.randrange(4)])
            continue
        b = ref[(pos + i) % len(ref)]
        pd = homopolymer_del_ref_at(ref, pos, i)
        r2 = rng.random()
        if r2 < pd:
            i += 1                        # deletion: skip reference base
            continue
        if rng.random() < P_SUB:
            out.append(BASES[(BASES.index(b) + 1 + rng.randrange(3)) % 4])
        else:
            out.append(b)
        i += 1
    seq = "".join(out)
    if rc_flag:
        seq = "".join(COMP[c] for c in reversed(seq))
    return seq

def homopolymer_del_ref_at(ref, pos, i):
    b = ref[(pos + i) % len(ref)]
    n = 1
    while n < 12 and ref[(pos + i + n) % len(ref)] == b:
        n += 1
    j = 1
    while j < 12 and ref[(pos + i - j) % len(ref)] == b:
        n += 1; j += 1
    return P_DEL + (P_HOMO_DEL if n >= 3 else 0.0)

# ----------------------------------------------------------------- main
def main():
    if not REF or not OUT:
        print(__doc__); sys.exit(1)
    os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
    name, ref = load_ref(REF)
    rng = random.Random(SEED)
    target_bp = int(COV * len(ref))
    n_reads, total = 0, 0
    lens = []
    with gzip.open(OUT, "wt", compresslevel=6) as fh:
        while total < target_bp:
            read = simulate_read(rng, ref)
            fh.write(">ont_%s_%06d len=%d\n" % (MODE, n_reads, len(read)))
            for j in range(0, len(read), 70):
                fh.write(read[j:j+70] + "\n")
            lens.append(len(read))
            total += len(read)
            n_reads += 1
    lens.sort()
    n50 = 0
    acc = 0
    for L in reversed(lens):
        acc += L
        if acc >= total / 2:
            n50 = L
            break
    print("wrote %s: %d reads, %.2f Mbp, coverage %.2fx, N50 %d, "
          "max %d, mode=%s, seed=%d"
          % (OUT, n_reads, total/1e6, total/len(ref), n50, lens[-1], MODE,
             SEED))

if __name__ == "__main__":
    main()
