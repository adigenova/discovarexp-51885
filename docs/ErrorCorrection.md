# Short-read error correction in DiscovarDenovo (discovarexp-51885)

This document describes, at the code level, how Illumina short-read errors are
handled by the `DiscovarExp` assembler in this repository. It is written for
developers/maintainers of the codebase; the user-facing manual is
`docs/DiscovarExp`.

## 1. Design overview

There is no single "error correction" executable. Instead, the assembler uses a
**layered, quality-score-driven strategy** in which errors are (a) kept out of
the assembly graph in the first place, and (b) actively fixed in *local* read
sets when local assemblies are built.  The general rules used throughout:

* Correction is **evidence-weighted by quality scores**: every algorithm sums
  base quality scores (Q) per candidate base rather than counting votes.
* Corrected bases are marked as unreliable: whenever a base is changed, its
  quality score is set to **0** (see `apply_corrections`,
  `pre_correct_parallel`, `readstack::CorrectAll`, and `ZeroCorrectedQuals`).
* Confirmed bases are protected from later changes ("locking"), so
  contradictory corrections become visible as *conflicts* rather than silently
  overwriting each other.
* Reads are never globally rewritten by the `DiscovarExp` main pipeline;
  global error suppression happens through **trusted-kmer filtering**, and
  base-level correction happens inside **local (gap) assemblies**.

Entry point of the executable:

```
src/paths/long/large/DiscovarExp.cc   -> GapToyCore(argc,argv)
src/paths/long/large/GapToyCore.cc    -> the actual pipeline
```

## 2. Where correction happens in the pipeline

```
BAM/fastq files
  |
  v
ExtractReads                      (input sanitization; N -> A)
  |
  v
buildReadQGraph                   (implicit error suppression: quality-trimmed
  |                                reads, K=60 trusted-kmer graph, read placement)
  |
  v  (optional PRECLOSE=True)
Preclose                          (local K=60 gap closing; each local read set is
  |                                pre-corrected with pre_correct_parallel,
  |                                FillPairs, Correct1Pre)
  v
Repath to K=200 -> Clean200        (graph hygiene; removes error-induced weak branches)
  |
  v
AssembleGaps2 -> MakeLocalAssembly1   (per-gap local assemblies)
  |                                      |
  |                                      v
  |                                 CorrectionSuite        (the full local
  |                                      |                  read-correction stack:
  |                                      |                  Precorrect -> FillPairs
  |                                      |                  -> Correct1Pre ->
  |                                      |                  CorrectPairs1)
  v                                      v
  LongHyper (local assembly from corrected reads) -> patch global graph
```

Key call sites:

| Stage | File | Function |
|---|---|---|
| Read extraction | `src/paths/long/large/ExtractReads.cc` | `ExtractReads` (called from `GapToyCore.cc:738`) |
| Trusted-kmer graph | `src/paths/long/BuildReadQGraph.cc` | `buildReadQGraph` (`GapToyCore.cc:799`) |
| K=60 pre-closing | `src/paths/long/large/Preclose.cc` | `Preclose` (`GapToyCore.cc:817`) |
| Local gap assembly | `src/paths/long/large/GapToyTools.cc` | `MakeLocalAssembly1` -> `CorrectionSuite` (line 1222) |
| Correction suite | `src/paths/long/LoadCorrectCore.cc` | `CorrectionSuite` (line 452) |

## 3. Input sanitization

`ExtractReads` (`src/paths/long/large/ExtractReads.cc`) parses BAMs (via
`bam/ReadBAM.h`, optionally `USE_PF_ONLY`) or paired fastq files into
2-bit `BaseVecVec` + `VecPQVec` quality vectors, written as
`frag_reads_orig.fastb` / `.qualp`.  No alignment information from the BAM is
used.  In the fastq path, `N` bases are replaced by `A`
(`ExtractReads.cc:392-395`), which together with low quality scores makes them
harmless downstream.  This is *sanitization*, not correction, but it is the
first place read errors are neutralized.

## 4. Implicit correction: the trusted-kmer graph (BuildReadQGraph)

`src/paths/long/BuildReadQGraph.cc` builds the K=60 unipath graph *only from
error-free-looking data*, which is the main reason raw sequencing errors never
enter the global assembly graph:

1. **Quality trimming of read tails** — `GoodLenTailFinder`
   (`BuildReadQGraph.cc:59-83`) walks each read from the 3' end backwards; a
   base with `Q < MIN_QUAL` (default **7**, `GapToyCore.cc:278`) resets the run
   length; the "good length" of a read is the prefix ending at the last window
   of K=60 consecutive good-quality bases.  Only that prefix is kmerized.
2. **Kmer frequency filter** — `Kmerizer::reduce`
   (`BuildReadQGraph.cc:129-133`) keeps only canonical K=60 kmers seen at
   least `MIN_FREQ` times (default **3**, `GapToyCore.cc:277`).  A random
   sequencing error creates kmers seen once, so error kmers are discarded;
   `pDict->recomputeAdjacencies()` repairs adjacencies broken by the removed
   kmers.
3. **Read placement instead of read sequence** — the trusted kmers are
   extended into unipath edges (`EdgeBuilder`), and reads are re-expressed as
   paths through those edges (`Pather`, `ShortKmerReadPather`, then
   `FixPaths`).  Read bases that disagree with the graph simply do not
   participate in the graph; disagreements show up later as "captured gaps"
   in read paths.

This is why the pipeline tolerates the raw error rate of 250-base reads
without a global read-corrector pass.  The construction of this graph itself
(unitig walks, HyperBasevector assembly, and the K=200 re-unipathing built
on top of it) is documented separately in `docs/GraphConstruction.md`.

## 5. PreCorrect: K=25 center-base voting (`pre_correct_parallel`)

The actively-correcting workhorse.  Code:

* driver: `src/paths/FindErrorsCore.{h,cc}` — `pre_correct_parallel()`
* kernel: `src/kmers/naif_kmer/KernelPreCorrector.h` — class `PreCorrector`
* framework: `src/kmers/naif_kmer/NaifKmerizer.h` — `naif_kmerize()`
* "hollow kmer": `src/kmers/naif_kmer/Kmers.h` — class `Kmer29H`

Algorithm (comment block at `FindErrorsCore.cc:31-53`, verified in code):

1. All K-mers (K=25, odd; `ForceAssert(_K & 1)`) of all reads are generated in
   canonical form, tagged with read id, start offset, rc flag and palindrome
   flag (`PreCorrector::parse_base_vec`).
2. `naif_kmerize` is a memory-bounded map-reduce: kmers are buffered into
   "parcels", sorted, and identical kmers are brought together.  The trick
   here: the record type is `KmerBVLocBase<Kmer29H>`, and **`Kmer29H`
   compares/hashes kmers with the central base masked out** (`flanks()`,
   `Kmers.h:275-297`).  Hence each sorted pile contains all reads' 25-mers
   that share the same flanking 12+12 bases but may differ at the central
   (13th) base.
3. `PreCorrector::summarize()` processes each pile (`KernelPreCorrector.h:350`):
   * Piles with fewer than `n_kmers_min` (**6**) rows, or piles whose flank is
     a palindrome, are ignored.
   * `BaseColumn` tallies, per candidate base A/C/G/T of the central column,
     the quality-sum `q_sum[4]`, the count of high-quality calls
     (`n_q_high[4]`, qual >= `q_high` = **20**) and the full
     per-base/per-quality histogram `nq[4][101]`.
   * Winner = base with max `q_sum`, valid only if its `q_sum >=
     q_sum_min_to_win` (**60**) (`BaseColumn::process`).
   * Loser (to be fixed): base != winner with `n_q_high <=
     n_q_high_max_to_lose` (**1**) and `q_sum < ratio_loser` (**0.25**) x
     winner `q_sum`.  Every read row carrying a loser base gets a
     *recommendation* (read id, position, replacement = winner base, rc-aware
     complementation) buffered per-thread, and periodically flushed under
     block locks into a copy of the reads (`flush_recommends`).
   * As a by-product, the same piles update a `KmerSpectrum` (kmer frequency
     histogram) which is merged at the end and used elsewhere to judge
     coverage/error rates.
4. Back in `pre_correct_parallel` (`FindErrorsCore.cc:82-145`) recommendations
   are applied with an **isolation rule**: proposed corrections that lie
   within `K/2` (**12**) bases of another proposed correction in the same read
   are both skipped (the flanks that built the pile may themselves be wrong).
   Applied corrections set the base and its quality to **0**, and statistics
   (reads corrected, corrections per read, corrections skipped,
   `frac_bases_pre_corrected` PerfStat) are logged.

Parameters live in `PC_Params` (`FindErrorsCore.h:19-37`) with the defaults
shown above.  Template-instantiated for `QualVecVec` and `QualNibbleVecVec`.

Call sites: `CorrectionSuite` (`LoadCorrectCore.cc:514`, K_PC=25),
`Preclose.cc:289`, and (with different heuristics flags)
`paths/long/PreCorrectAlt1.cc` / `paths/long/PreCorrectOldNew.cc` as
alternatives.

## 6. FindErrors: kmer-stack correction of whole reads (`find_errors_parallel`)

A more powerful (and more expensive) corrector that uses kmer-sharing *stacks*
of entire reads.  Code:

* driver: `src/paths/FindErrorsCore.{h,cc}` — `find_errors_parallel()` /
  single-threaded `find_errors()`
* kernel: `src/kmers/naif_kmer/KernelErrorFinder.h` — class `ErrorFinder`
* column statistics: `ColumnQualityScores` (same file)

Algorithm (`find_errors_parallel`, `FindErrorsCore.cc:404-498`, run for
`NUM_CYCLES` iterations):

1. **Quality capping**: `cap_quality_scores` replaces every quality score by
   the minimum in a window of radius `qual_ceil_radius` (**2**)
   (`FindErrorsCore.cc:209-267`), every cycle (`qcr_always` default true).
   This stops one high-quality base from "blessing" a bad neighborhood.
2. All K-mers (K up to 124; template dispatch on `Kmer29/Kmer60/Kmer124`)
   are piled up with `naif_kmerize` exactly as in PreCorrect, but with the
   *full* kmer as the key.
3. `ErrorFinder::summarize` filters stacks by depth (`min_readstack_depth`
   **5** .. `max_readstack_depth` 0 = infinity, `min_basestack_depth`
   **1**) and skips palindromic stacks (unless `do_palindromes`), then:
   * `build_column_stats` (`KernelErrorFinder.h:437-467`) lays each whole read
     of the stack into a common coordinate system (the kmer occupies columns
     `[nb_max-K, nb_max)`), building per-column `ColumnQualityScores`
     (per-base QSsum, count, and count of bases with `Q >= min_Q_to_support`
     = **20**).  Both fw and rc rows are handled by complementation.
   * `process_column_stats` walks columns left and right from the kmer
     (under-kmer columns are skipped by default via `skip_under`), and stops
     at the first **branch** (column with an alternative supported base,
     `found_alt_support`) unless `do_branches`.
   * `ColumnQualityScores::process` (`KernelErrorFinder.h:55-165`):
     - *recommendation* (replace base X by winner W) iff W has the highest
       QSsum >= `min_QSsum_to_win` (**40**), X is unsupported (QSsum <
       `min_QSsum_to_support` **60**), `QSsum(X)/QSsum(W) <=
       max_QSsum_ratio_to_correct` (**0.25**), and the column overall is
       unambiguous (`full_max_QSsum_ratio <= max_QSsum_ratio_to_correct2`
       **0.25**).
     - *confirmation* (lock the base; immune to later correction) iff no
       alternative support exists, the winner is supported or the column is
       outside the kmer, ratios are within `max_QSsum_ratio_to_confirm`
       (**1.0**) / `max_QSsum_ratio_to_correct`, and the winner's evidence is
       `min_QSsum_to_confirm` (**90**) or `min_n_bases_to_confirm` (**3**)
       high-quality bases.  `auto_confirm` (default off) confirms any base
       with QSsum above a threshold.
   * `confirm_and_recommend` emits per-read, per-column records; the flush
   step (`flush_recommends_confirms`, `KernelErrorFinder.h:325-424`) applies
   **confirmations first** (lock the base in `base_locked`, revert any
   pending correction), then **recommendations**, discarding positions that
   receive two different recommendations (the position is locked; both are
   counted as conflicts at the end).
4. `apply_corrections` (`FindErrorsCore.cc:344-397`) commits the recommended
   bases into the reads: unlocked changed positions are rewritten with
   quality 0 (`n_corrections`); locked-and-changed positions are counted as
   `n_conflicts`; locked-and-unchanged as `n_confirmed`.  Statistics are
   logged per cycle via PerfStat.

All thresholds live in `EF_Params` with defaults in `efp_default`
(`FindErrorsCore.h:44-130`).

In this snapshot, `find_errors`/`find_errors_parallel` are library code (the
header notes: "useful when correcting small local read sets in modules other
than FindErrors") and are **not** invoked by the `DiscovarExp` pipeline; the
PreCorrect kernel (Section 5) is what runs.

## 7. Friend-stack correction (`Correct1Pre`, `readstack`)

`src/paths/long/Correct1Pre.cc` corrects reads using alignments to
"friends" — other reads sharing a kmer:

1. `FriendAligner` (`src/paths/long/FriendAligns.h`; implementation
   `FriendAlignFinderQ.{h,cc}`) indexes kmers (respecting kmer
   *context*, so reads sharing a unique kmer are found) and returns, for each
   read, the list of friends with offsets (`FriendAlignFinderQ.h:63-71`).
2. For each editable read (`to_edit`), `Correct1Pre.cc:100-163` builds a
   `readstack` (`src/paths/long/ReadStack.{h,cc}`): rows = the read plus its
   aligned friends, columns = alignment positions.
3. Suspect friends are pruned with `stack.HighQualDiff(30, 1, suspect)`
   (friends disagreeing with high-quality consensus) and optionally
   `HighQualDiffWindow`.
4. Correction itself, `readstack::CorrectAll` (`ReadStack.cc:1017-1079`):
   * per column, quality-sum per base, counting Q<=2 bases as ~0.2 of a
     vote; the best competitor's single strongest quality score is deducted
     from its sum;
   * the column is "undisputed" iff winner sum >= 50, >= 10x runner-up, and
     runner-up <= 100;
   * the read's base is replaced by the winner (quality set to 0); the first
     disputed column defines a safe trim length `trim_to` (used to trim the
     read before local assembly).
   * the EM alternative `readstack::CorrectAllEM3` (`ReadStack.cc:1000`,
     model in `src/paths/long/EMEC3.{h,cc}`) computes probabilistic
     consensus instead.
5. Corrected reads/quals are written back to `bases`/`quals`
   (`Correct1Pre.cc:165-172`).

Heuristics: `FF_MIN_FREQ`/`FF_MAX_FREQ`/`FF_MIN_QUAL`/`FF_COVERAGE` control
friend finding; `MAX_STACK` skips reads with too many friends;
`USE_EM` selects the EM corrector.  All are set in
`src/paths/long/Heuristics.h` (`long_heuristics`).

## 8. Pair-based correction and the local CorrectionSuite

`CorrectionSuite` (`src/paths/long/LoadCorrectCore.cc:442-...`) orchestrates
the complete local read-correction stack used before each local (gap)
assembly:

1. **PreCorrect** — `pre_correct_parallel(K=25)` as in Section 5
   (`LoadCorrectCore.cc:511-515`), unless `heur.PRECORRECT_ALT1` or
   `heur.PRECORRECT_OLD_NEW` select the alternate/legacy precorrectors.
2. **`ZeroCorrectedQuals`** (`LoadCorrectCore.cc:415-427`) — zeroes the quality
   scores of every base that differs from the original reads, so downstream
   stages know which bases are trustworthy.
3. **`FillPairs`** (`src/paths/long/FillPairs.{h,cc}`) — finds perfect kmer
   matches between the two reads of a pair (MIN_FREQ=5) to *fill* the gap
   between them, creating synthetic "filled" reads that span the ~450 bp
   fragment; these get flat Q40 and are marked `done` (excluded from further
   correction) — see `LoadCorrectCore.cc:546-585`.
4. **`CapQualityScores`** (`LoadCorrectCore.cc:429-440`) — sliding minimum
   with radius **4** (a stronger version of Section 6.1).
5. **`Correct1Pre`** (Section 7) over reads not yet `done`
   (`LoadCorrectCore.cc:600-604`); trimmed reads are written to
   `frag_reads_mod0.fastb`.
6. Precorrected reads are themselves assembled at K=80
   (`LongReadsToPaths`) and consistent pairs are *closed*
   (`LoadCorrectCore.cc:643-747`): pairs whose two paths through the local
   graph agree are replaced by the reconstructed fragment sequence (Q40).
7. **`CorrectPairs1`** (`src/paths/long/CorrectPairs1.cc`, entry
   `LoadCorrectCore.cc:749-757`) — the final, most aggressive step: builds
   read stacks around placed pairs and produces **corrected reads as
   efasta** (`VecEFasta& corrected`), i.e. where the evidence is split the
   output keeps ambiguity codes like `{A,C}` rather than guessing.
   Optionally a second pass (`heur.CP2`) corrects the remaining "special"
   reads.
8. `DefinePairingInfo` re-establishes the pairing of the corrected reads,
   which then feed the local assembly (`LongHyper`, `GapToyTools.cc:1153`).

The same suite in a lighter form (steps 1, 3, 4, 5 plus a K=60 local
unipath rebuild) is used by `Preclose` at the K=60 stage
(`src/paths/long/large/Preclose.cc:278-343`) when `PRECLOSE=True`.

## 9. Graph-level hygiene and diagnostics

* **`Clean200`/`Clean200x`** (`src/paths/long/large/Clean200.cc`): after
  repathing to K=200, walks the graph and deletes *weak branches* —
  extensions whose read-path support is too low — i.e. it removes error
  artifacts that survived into the graph (`GapToyCore.cc:876-881`).
* **`SelectSpecials`** (`src/paths/long/large/GapToyTools3.cc:31`): a
  diagnostic (flag `SELECT_SPECIAL`) that identifies reads that could not be
  placed on the K=200 graph and "ought to be error corrected", by looking for
  unplaced reads with plausible quality-sum alignments to graph edges
  (`GapToyCore.cc:980-984`).
* **`EvalCorrected`** (`src/paths/long/EvalCorrected.{h,cc}`) evaluates the
  accuracy of corrected reads against a reference (development tool).
* **`PerfStat` reporting**: `frac_bases_pre_corrected`,
  `frac_bases_confirmed_<cycle>`, `frac_bases_corrected_<cycle>`,
  `frac_conflicts_<cycle>` are emitted to the perf-stats log
  (`FindErrorsCore.cc:163-169, 473-484`).

## 10. Related code that is *not* correction

* `src/ReadError.{h,cc}` — a *record type* describing an individual read error
  (substitution / insertion / deletion, with read and reference base) that is
  used by the read simulator (`src/simulation/ReadTemplate.h`,
  `src/simulation/ReadSimulatorSimpleCore.h`) to inject errors into reads for
  testing assemblers; it does not fix anything.
* `src/paths/long/ultra/*` — the "ultra" correction machinery
  (`ThreadedBlocks`, `Prefab`, `CorrectSomeMoreReads`) belongs to the
  long-read (HALC-style) pipeline also present in this tree, not to the
  short-read `DiscovarExp` path.

## 11. Quick reference

| What | Symbol | File | Default |
|---|---|---|---|
| Precorrect kmer size | `K_PC` | `LoadCorrectCore.cc:512`, `Preclose.cc:287` | 25 |
| Min pile size | `PC_Params::n_kmers_min` | `paths/FindErrorsCore.h` | 6 |
| Min winner quality-sum | `PC_Params::q_sum_min_to_win` | " | 60 |
| Loser ratio | `PC_Params::ratio_loser` | " | 0.25 |
| Max loser high-Q calls | `PC_Params::n_q_high_max_to_lose` | " | 1 |
| High-Q threshold | `PC_Params::q_high` | " | 20 |
| Correction isolation window | `K/2` | `FindErrorsCore.cc:120` | 12 |
| Trusted kmer size | `K` | `BuildReadQGraph.cc:51` | 60 |
| Trusted kmer min freq | `MIN_FREQ` | `GapToyCore.cc:277` | 3 |
| Good-tail qual floor | `MIN_QUAL` | `GapToyCore.cc:278` | 7 |
| Friend-stack min win | `readstack::CorrectAll` | `ReadStack.cc:1031` | 50, 10x, 100 |
| Qual capping radius (local) | `CapQualityScores` | `LoadCorrectCore.cc:430` | 4 |
| Qual capping radius (FindErrors) | `EF_Params::qual_ceil_radius` | `FindErrorsCore.h` | 2 |

To experiment with correction heuristics, the two parameter structs to look
at are `PC_Params` and `EF_Params` (both in `src/paths/FindErrorsCore.h`),
and `long_heuristics` (src/paths/long/Heuristics.h) for the suite-level
switches (`PRECORRECT_ALT1`, `PRECORRECT_OLD_NEW`, `CORRECT_PAIRS`, `CP2`,
`PRECORRECT_SEQ`, `FF_*`, `MAX_STACK`, `USE_EM`).
