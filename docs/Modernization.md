# Discovar modernization roadmap: shallow long-read hybrid assembly and phasing

This document records the modernization plan for the DiscovarDenovo
(`discovarexp-51885`) assembler: the strategic direction, the concrete
integration points in this codebase, a phased workplan, and the feasibility
math behind the target data model.

Companion documents:

* `docs/GraphConstruction.md` — how the K=60 → K=200 unitig graph is built,
  and the verified memory profile of the a.200 stage (improvements M1–M7).
* `docs/ErrorCorrection.md` — the read error-correction machinery
  (trusted-kmer filtering, PreCorrect, the local CorrectionSuite).

---

## 1. Strategic vision

**Goal: a chromosome-scale, haplotype-aware assembler built on the Discovar
short-read graph core, using wengan-style shallow long-read integration.**

Target data model (the working hypothesis, H1):

> **60× Illumina PE150** (PCR-free, ~350–450 bp fragments) **+ ~10× ONT**
> (Q20-class chemistry) — with PacBio HiFi (~10–15×) as an alternative or
> complementary long-read layer.

The combination is the point: short reads are cheap, accurate, and dense
enough to build a near-perfect *sequence* graph; even shallow long reads
provide *structure* — repeat spanning, scaffold links, gap spans, and phase
links — that no amount of short-read coverage can supply. Wengan (the
author's hybrid assembler; "Wengan: a powerhouse for chromosome-scale genome
assembly", Nat. Commun. 2022) demonstrated that this division of labor can be
implemented extremely efficiently. The vision here is to run the same hybrid
playbook on top of a much stronger short-read core.

### 1.1 Why Discovar is the right core

| Discovar property | Why it matters for the hybrid design | Where |
|---|---|---|
| Errors never enter the global graph | long reads can inform *structure* without polluting *sequence* — the same philosophy extends naturally | `docs/ErrorCorrection.md` |
| Graph-preserving output (lines/cells, `{A,C}` alternatives) | the substrate for phasing already exists: heterozygous sites are cells with two paths | `Lines.cc`, `FinalFiles.cc`, `docs/DiscovarExp` |
| Local gap assembly with full correction stack | gap closing upgrades from "short reads only" to "long-read-guided, short-read-finished" | `AssembleGaps.cc`, `GapToyTools.cc` (`MakeLocalAssembly1`), `LoadCorrectCore.cc` (`CorrectionSuite`) |
| Deterministic, symmetric graph (involution) | scaffold/phase modules can trust edge-level RC invariants | `HBVFromEdges.cc`, `TestInvolution` |
| Complete read→graph placement (`a.200.paths`) | fragment-pair phase links are already computable from existing artifacts | `Repath.cc`, `GapToyCore.cc:885-898` |

### 1.2 Design principles

1. **Long reads inform structure; short reads supply sequence.**
   ONT bases (2–5% raw error) never enter the trusted graph — they are used
   only through *alignments*: scaffold links, gap spans, repeat bridges,
   phase threads. All output sequence comes from short reads (optionally
   HiFi, at Q30, may be treated as near-trusted consensus for gap fill).
   This is the wengan pattern, and it is a direct extension of Discovar's
   existing error philosophy.
2. **Preserve alternatives, don't collapse them.** The lines/cells machinery
   is kept as the output substrate; phasing orders the alternatives rather
   than eliminating them.
3. **Stay cheap.** The value proposition is chromosome-scale assembly on
   affordable data and modest hardware. The memory fixes (M1–M2) and the
   efficient wengan-style integration are part of the same goal.

---

## 2. Hybrid assembly design (H1)

### 2.1 Current pipeline, overlay of new modules

```
 reads PE150 (60x)
    │ buildReadQGraph (K=60) → Repath (K=200) → Clean200     [unchanged core]
    ▼
 a.200 graph (unitigs + read paths)
    │ Lines.cc: lines and cells                                [unchanged core]
    ▼
 lines / cells (alternatives preserved)
    │
    ├── LR-1  long-read mapping to the graph          (new: minimap2-style)
    ├── LR-2  long-read scaffold graph over lines     (new: augments MakeGaps)
    ├── LR-3  long-read-guided gap filling            (new: feeds AssembleGaps)
    │
    ├── PH-1  pair-based phase links                  (new: from a.200.paths)
    ├── PH-2  ONT-read phase threads across cells     (new)
    ├── PH-3  phase block construction + QA           (new)
    ▼
 phased lines → a.final (haplotagged GFA, phased efasta/VCF)
```

### 2.2 Long-read ingestion (LR-0)

`ExtractReads` (`src/paths/long/large/ExtractReads.cc`) already accepts
`.bam`, `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz`, `.fastb` (line ~201) and
converts `N` → `A` with low quality — acceptable for 2-bit storage of ONT
bases used *only* as alignment templates. Long reads are loaded into their
own `vecbvec` and **never routed through** `buildReadQGraph`,
`GoodLenTailFinder`, or any trusted-kmer logic. At 10× on a 3.1 Gb genome
(~31 Gbp) this costs ~8 GB in 2-bit — modest.

### 2.3 Long-read mapping (LR-1)

Two options, recommend A first:

* **(A) External mapper as a subprocess** — minimap2 with `-x map-ont` /
  `-x map-hifi` against `a.200/a.fastb`, parsed from PAF (or BAM via the
  existing `src/bam` readers). Discovar already shells out to external tools
  via `SystemSucceed` (e.g. `MakeLookupTable`, `GapToyCore.cc:791-794`;
  SAMtools in `ExtractReads`), so the dependency pattern — and vendoring it
  under `deps/` like jemalloc — is established practice.
* **(B) Internal coarse aligner** seeded on trusted kmer anchors — avoids
  the dependency but is a significant engineering effort with no accuracy
  advantage over minimap2 for this purpose.

Mapping targets the *edge sequences* of `a.200`, then projects alignments
onto the graph (edges → vertices), retaining: read id, edge spans with
orientation, soft-clip/gap coordinates at edge transitions. Everything
downstream consumes this projection, never raw ONT bases.

### 2.4 Long-read scaffolding (LR-2)

Replace the pair-link model of `MakeGaps` (`src/paths/long/large/MakeGaps.cc`,
currently `MIN_LINK_COUNT=3` pairs, `MIN_LINE=5000`, `GapToyCore.cc:339-342`)
with a two-evidence scaffold graph over **lines** (not raw edges):

* nodes: lines (and standalone edges), with involution pairs;
* links: long-read threads — an ONT read spanning lines `A … B` in known
  orientation contributes a link with span estimate; multiple consistent
  reads accumulate link weight (coverage-robust threshold, Poisson-calibrated
  instead of a fixed `MIN_LINK_COUNT`);
* conflict resolution: links that disagree (repeat-induced) are resolved by
  uniqueness of flanks and short-read pair support; wengan-style
  conservative breaking when evidence conflicts.

Keep the existing gap representation (captured gaps, 100-N cells) so the
downstream patching/finishing machinery is reused unchanged.

### 2.5 Repeat bridging

K=200 unitigs already resolve repeats < 200 bp; the assembly's contiguity
ceiling is repeats ≥ 200 bp. Every repeat shorter than the ONT read length
is bridgeable in principle; expected *spanning* coverage for a repeat of
length `r` is `C·E[max(0, L−r)]/E[L]` (see §7 math). At 10× with N50 ~30 kb,
repeats up to ~15–20 kb bridge with several reads. This covers Alu, LINE1,
most microsatellite/segmental-duplication-class obstacles; centromeres and
>30 kb segmental duplications remain open (see risks, §6).

### 2.6 Long-read-guided gap filling (LR-3)

Discovar's gap closing is already local-assembly-based (`AssembleGaps`,
`MakeLocalAssembly1` → `CorrectionSuite` → local `LongHyper`). The upgrade:

1. find ONT reads spanning a captured gap (from the LR-1 projection);
2. use the read's gap-spanning segment to *recruit* short-read pairs
   straddling the gap (the read supplies coordinates, not sequence);
3. run the existing local assembly + correction stack on the recruited set;
4. where the local assembly fails (low coverage, cycle), fall back to the
   ONT consensus segment polished by short reads (Racon-style, one pass),
   marked as low-confidence in output.

This reuses the crown-jewel machinery and adds long-range recruitment the
current pair-window (`MAX_PROX_LEFT/RIGHT=400`) cannot provide.

---

## 3. Phasing design (PH-1..3)

### 3.1 What already exists

* Cells with two (occasionally more) paths across them *are* heterozygous
  sites (`Lines.cc`, `MAX_CELL_PATHS=50`, `MAX_DEPTH=10`,
  `GapToyCore.cc:334-335`). The manual notes partial phasing already causes
  >2 paths in some cells.
* Read pairing is fully tracked (`PairsManager`; path-level pair ids via
  `DefinePairs`, `GapToyTools.cc:1010-1046`, used by `Preclose.cc:179-223`),
  and every read's placement on the K=200 graph is in `a.200.paths`.

### 3.2 Phase links

* **Fragment pairs (PH-1)**: a pair whose reads traverse cells c and c′
  gives a short-range phase link between path choices. Range ≈ fragment size
  (~450 bp → adjacent cells). Frees existing artifacts; no new mapping.
* **ONT threads (PH-2)**: one 30 kb ONT read at human heterozygosity threads
  ~20–30 cells; its projection votes for one path choice per cell, and —
  critically — links the *choices across cells* into a consistent haplotype.
  This is sparse long-range phasing with the same data bought for
  scaffolding — no extra sequencing.

### 3.3 Phase block construction (PH-3)

Formally: a graph whose nodes are (cell, path) choices and whose edges are
pair links (weight: pair count) and ONT threads (weight: read count ×
aligned fraction). Resolve by a two-pass scheme:

1. greedy seed-and-extend along high-weight ONT threads (haplotype
   "spines");
2. local refinement (EM or max-cut heuristic) using pair links;
3. QA: switch-error proxy from pair consistency (pairs disagreeing with the
   assigned phase), reported per line; phase block boundaries where evidence
  runs out.

Expected phase block N50 is governed by ONT read N50 (hundreds of kb to
low Mb at 30–50 kb N50) — far beyond what fragment pairs alone achieve.

### 3.4 Output

* two haplotype-resolved line sets (`HP:1`/`HP:2` per block, hifiasm-style
  GFA convention);
* phased efasta: alternatives tagged with phase block id instead of
  unordered `{A,C}`;
* VCF of cell alternatives vs a chosen reference path, with phased
  genotypes.

---

## 4. PE150 adaptation

The pipeline is tuned for 2×250 MiSeq (`docs/DiscovarExp`). For PE150:

| Item | Location | Change |
|---|---|---|
| branch-walk depth | `Clean200.cc:52` (`max_rl = 250`) | parameterize to read length (150) |
| extension walk depth | `GetExtensions` call, `Clean200.cc:62-64` | same |
| read-path extension heuristics | `ExtendReadPath`, `HBVPather` algorithm-two trimming (≤5-kmer back-off) | re-tune for 91 60-mers/read |
| pair-selection windows | `MAX_PROX_LEFT/RIGHT=400` (`GapToyCore.cc:343-344`) | ok for 350–450 bp fragments; verify |
| K=60 suitability | `BuildReadQGraph.cc:51` | fine for 150 bp reads (91 consecutive 60-mers); K=200 places unaffected (places are edge lists, not read-length limited) |
| unplaced reads | `SHORT_KMER_READ_PATHER` (`GapToyCore.cc:327`) | expect more unplaced reads → enable and re-tune |

Acceptance: `X=tiny`/`X=small` runs (`GapToyCore.cc:509-590`) on PE150 test
data reproduce K=250-class edge accuracy statistics.

---

## 5. Supporting modernization workstreams

These unblock or amplify the strategic plan (details in
`docs/GraphConstruction.md` §7 where marked):

1. **Engineering unblock (P0, hard prerequisite)**
   * port the toolchain beyond GCC 4.8 (README constraint) — C++17/20,
     CMake build, containers;
   * self-contained tests + CI: the `X=` evaluation modes currently hard-code
     Broad fileserver paths (`GapToyCore.cc:671`) — replace with bundled
     mini-datasets and golden `a.200` outputs; property tests for involution
     symmetry and kmer→edge totality;
   * **memory fixes M1/M2** (dictionary right-sizing + unloading reads
     around `Repath`): peak ~250–450 GB → ~130–160 GB, enabling the
     60×-human data model on ordinary 256 GB machines.
2. **GFA writer (P0.5)**: emit `a.200` edges as segments (canonical edges;
   RC via GFA1.1 links; segment overlaps are exactly K−1 = 199, consistent
   with the circular-sequence convention in `docs/DiscovarExp`), lines as
   paths. Enables Bandage/odgi/ragtag interop and is the substrate for the
   haplotagged phased output.
3. **Translation rewrite (M5)**: direct edge-list emission from `Pather`,
   removing the "very ugly" `Repath.cc:139-209` machinery and its known
   indexing bug — should land before new consumers (PH-1) depend on
   `a.200.paths`.
4. **Adaptive multi-K** (later): data-driven K selection from the kmer
   spectrum (`CN1PeakFinder` exists) instead of the fixed 60→200 jump;
   helps non-human coverage profiles.
5. **Probabilistic `Clean200`** (later): coverage-model-based weak-branch
   removal instead of fixed vote thresholds, for robustness at 20×–90×.

Explicitly *not* pursued: rewriting the local-assembly machinery (it is
the core asset), and competing head-on with HiFi-only assemblers (hifiasm
owns that niche; the value here is the short-read + shallow-LR regime).

---

## 6. Risks and mitigations

| Risk | Mitigation |
|---|---|
| ONT Q20 errors degrade links inside repeats | link only via unique flanking anchors; require multi-read consistency; fall back to pair evidence |
| 10× coverage non-uniform (pore bias, ultra-low regions) | coverage ladder (5×/10×/20×) benchmarks; scaffold conservatively where links are sparse |
| "Chromosome-scale" not reached at 10× N50 30 kb for repeat-heavy genomes | levers: ultra-long ONT boost, optional Hi-C layer, HiFi complement; measure scaffold N50 vs read N50 in P6 to characterize the ceiling |
| PE150 reduces read-path density → more captured gaps | `SHORT_KMER_READ_PATHER`, `FILL_JOIN`, tuned local recruitment; long-read recruitment (LR-3) compensates |
| phase switch errors in low-heterozygosity regions | blocks simply end; report block N50 and pair-consistency QA, don't force |
| external minimap2 dependency | vendored under `deps/` (jemalloc precedent); internal coarse aligner (LR-1 option B) as fallback |
| memory at human scale | M1/M2 in P0 (see `docs/GraphConstruction.md`) |
| licensing of modified redistribution | review `LICENSE` (Broad terms) before public release of the modernized build |

---

## 7. Feasibility math (H1 assumptions)

Assume a 3.1 Gb diploid genome, PE150 at 60×: ~620 M read pairs, ~186 Gbp;
short-read storage ~47 GB (2-bit) — fits the post-M1/M2 memory envelope.

**ONT layer.** 10× ≈ 31 Gbp. With N50 30 kb (Q20-class chemistry, standard
prep; ultra-long preps push N50 ≥ 100 kb): ~1 M reads.

*Spanning coverage* of a repeat of length `r`:

```
E[spanning depth] ≈ C · E[max(0, L − r)] / E[L]
```

* `r = 5 kb`:  ≈ 8× spanning → solid bridges.
* `r = 15 kb`: ≈ 5× → several supporting reads.
* `r = 25 kb`: ≈ 1.7× → marginal; needs luck or higher N50.
* `r ≥ 30 kb` (centromeric arrays, long segmental duplications): open —
  same as for any read ≤ N50 technology.

*Scaffold N50* statistically scales with read N50; with consistent
multi-read threading it exceeds it (empirically hybrid scaffolders reach
scaffold N50 of several × read N50 at ≥10×). The honest position for H1:
**scaffold N50 in the multi-Mb range is the base case; true chromosome-scale
(few breaks per chromosome) requires the long tail of the read length
distribution — ultra-long ONT, or a complementary long-range layer — and
this is a central quantity to measure in P6, not to assert.** Wengan's
published results indicate the hybrid regime can approach chromosome scale
on real data; the exact threshold is genome- and prep-dependent.

*Phasing.* Human heterozygosity ~1/1–1.5 kb ⇒ a 30 kb ONT read threads
~20–30 het cells ⇒ phase block N50 ~ ONT read N50, extendable by
multi-read stitching; fragment pairs give the fine-grained links within
blocks. Phase data costs nothing extra beyond the scaffolding layer.

*HiFi alternative.* 10–15× HiFi (N50 15–20 kb, Q30) bridges the same repeat
classes with shorter spans but near-perfect bases; it can additionally serve
as trusted consensus for gap fill (opt-in, deviation from Principle 1).

---

## 8. Phased workplan

| Phase | Deliverable | Acceptance criteria |
|---|---|---|
| **P0** engineering unblock | modern toolchain, CI, mini-datasets, M1+M2 memory diffs, GFA writer | builds on modern GCC; tests green; `X=tiny` reproduces current `a.200` checksum; human-scale peak ≤ ~160 GB |

**P0 status (partially landed, 2026-11):** the tree now builds with
modern compilers (GCC 16, `-std=gnu++14`; legacy dead flags removed,
jemalloc auto-detected, `LT_INIT`) and the M1/M2 memory diffs are in
(see `docs/GraphConstruction.md` §7). Also fixed en route: friend
declarations with default arguments (defaults relocated to namespace-scope
declarations), the `::align`/`std::align` collision (~190 type positions
qualified), a private-access violation in `PQVec.h`, an explicit-constructor
conversion in `KmerShape.h`, the `ifstream`→bool conversion in `System.cc`,
a non-conforming `std::map` comparator/allocator in `SAM2CRD.cc`, dead
template code in `KmerParcelsClasses.h`, a genuine `std::cout`-for-
`std::endl` typo in `GapToyTools4.cc`, and missing explicit template
instantiations for `ProcessFrequentKmers`. The **GFA1 writer** has landed
as well (`WriteGFA` in `FinalFiles.{h,cc}`, emitting
`a.final/a.lines.gfa`: doubled segments with `INV` reverse-complement
tags and `KC` coverage, all vertex adjacencies as `199M` links, scaffold
lines as paths, split at captured gaps). It is validated on the E. coli
K-12 50× MiSeq test (52 s end-to-end, 0 misassemblies, N50-perfect
1.27 Mb): all 478 segments byte-identical to `a.fasta`, all RC pairs
sequence-consistent, and the path-length N50 reproduces the reported
133,298 exactly. **M8** (need-based MapReduce buffer pool) also landed:
the same E. coli run now peaks at 1.1 GB instead of 6.8 GB (bit-identical
result), so small-genome runs no longer inherit the memory profile of
human-scale ones. Still open in P0: CI with bundled
mini-datasets, M3 (16-byte kmer records).
| **P1** PE150 adaptation | parameterized read length, re-tuned pathing | PE150 `X=small` ≈ PE250-class edge accuracy |
| **P2** LR-0/LR-1 | long-read loading + minimap2 projection onto `a.200` | ≥95% of ≥20 kb reads project onto ≥2 edges with consistent coordinates |
| **P3** LR-2 scaffolder | line-level long-read scaffold graph, Poisson link thresholds | scaffold N50 ≥ 5× PE-only baseline on test genome; ≤1 misjoin/Mb |
| **P4** LR-3 gap filling | ONT-guided recruitment into `AssembleGaps`/`CorrectionSuite` | captured-gap count reduced ≥50% on test genome; no quality regression in filled gaps |
| **P5** PH-1..3 | pair + ONT phase links, phase blocks, haplotagged outputs | phase block N50 ≥ 500 kb at 10× ONT; pair-consistency QA reported |
| **P6** benchmarking | NA12878-scale run + report vs baselines (wengan, SPAdes+long-read scaffolder, hifiasm on HiFi-only) | characterize H1: scaffold N50 vs ONT coverage/N50 curves, switch-error rates, BUSCO/QUAST-LG, runtime/memory |

Sequencing note: P2–P4 and PH-2 share one ONT dataset; no incremental data
cost across the plan.

---

## 9. Concept mapping: wengan → Discovar

| wengan concept | Discovar counterpart / integration point |
|---|---|
| efficient short-read contig set | `buildReadQGraph` → `Repath` → `Clean200` (`a.200`), with graph alternatives preserved — stronger substrate |
| long-read mapping to contigs | LR-1 projection onto `a.200` edges (external minimap2, `SystemSucceed` pattern) |
| long-read link graph for scaffolding | LR-2 over lines, replacing pair-only `MakeGaps` links; Poisson thresholds |
| long-read-assisted gap closing | LR-3 feeding the existing `AssembleGaps`/`CorrectionSuite` local stack |
| conservative conflict breaking | line-level uniqueness checks + short-read pair verification |
| (new here) graph alternatives as phasing substrate | `Lines.cc` cells + `a.200.paths` pairs + LR-1 threads → PH-1..3 |

---

*Roadmap based on discovarexp-51885 plus the analysis in
`docs/GraphConstruction.md` and `docs/ErrorCorrection.md`. Line numbers
refer to this source tree.*
