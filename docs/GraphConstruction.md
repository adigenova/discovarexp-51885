# Unitig graph construction (K=60 → K=200) and memory optimization

This document describes how the DiscovarDenovo assembly graph is built, with
particular focus on the **K=200 unitig graph** (`a.200`), and documents an
analysis of its memory footprint together with a ranked plan for reducing it.
All file and line references are to this source tree (`discovarexp-51885`).

Intended audience: developers modifying the assembly pipeline.

Companion documents: `docs/ErrorCorrection.md` covers the same K=60 stage
from the *error-suppression* perspective (quality trimming, trusted-kmer
filtering, local read correction); this document covers the *graph
construction* mechanics and the memory cost of the K=200 stage.
`docs/Modernization.md` is the development roadmap (shallow long-read
hybrid assembly and phasing) that builds on the memory plan in §7.

---

## 1. Pipeline overview

Entry point: `DiscovarExp` (`src/paths/long/large/DiscovarExp.cc`) calls
`GapToyCore()` (`src/paths/long/large/GapToyCore.cc`), which drives the whole
assembly. The two kmer sizes are hard-coded there:

```cpp
// GapToyCore.cc:444-445
const int K0 = 60;    // error-correcting unipath graph
const int K  = 200;  // assembly graph
```

The design (stated in the header comment of `GapToyCore.cc`, lines 9–17) is a
two-stage de Bruijn scheme:

```
 reads (BAM)
    │  ExtractReads                         GapToyCore.cc:738
    ▼
 bases + quals
    │  buildReadQGraph  (K=60)               GapToyCore.cc:799
    ▼
 hbv60: HyperBasevector + ReadPathVec (read paths on K=60 unitigs)
    │  FixPaths / [JoinPaths / Preclose]     GapToyCore.cc:837-854
    │  quals.unload()                        GapToyCore.cc:859
    ▼
 Repath  (K=60 → K=200)                      GapToyCore.cc:860
    │   places ─► basevectors ─► LongReadsToPaths(K=200)
    │   = re-unipath at K=200 + path translation
    ▼
 hb200: HyperBasevector + inv + paths        (a.200.hbv / .inv / .paths)
    │  Clean200 (error correction)           GapToyCore.cc:876-881
    ▼
 a.200/  (a.hbv, a.hbx, a.fastb, a.inv, a.paths, a.paths.inv, a.k)
    │  local gap assembly + patching → a.patched
    │  Simplify / lines → a.fin
    │  MakeGaps (scaffolding) → a.final
    ▼
 final assembly
```

Why two kmer sizes: 200-mers cannot be counted directly in noisy reads. The
K=60 graph — built only from quality-trimmed, sufficiently frequent kmers —
serves as the error-corrected representation of the trusted sequence. Reads
are mapped onto it as paths; the **distinct observed paths** are then
re-assembled at K=200 into long, unambiguous unitigs.

### File map

| Component | File | Role |
|---|---|---|
| Driver | `src/paths/long/large/GapToyCore.cc` | pipeline control, stage orchestration |
| K=60 graph | `src/paths/long/BuildReadQGraph.cc` | kmer dict, unitigs, read placement at K=60 |
| K=60 → K=200 | `src/paths/long/large/Repath.cc` | places, place sequences, path translation |
| K dispatcher | `src/paths/long/LongReadsToPaths.cc` | chooses builder; K=200 → BigK pather |
| K=200 unitigs | `src/kmers/BigKPather.cc` | 200-mer dict, unitig construction, re-pathing |
| 200-mer record | `src/kmers/BigKMer.h` | kmer as (sequence view + context + hash) |
| Kmer branching | `src/kmers/KMerContext.h` | 8-bit predecessor/successor bitmask |
| Graph assembly | `src/paths/long/HBVFromEdges.cc` | edges → vertices → HyperBasevector, HyperKmerPath |
| K templates | `src/paths/long/LargeKDispatcher.h` | compile-time K dispatch table |
| K=200 cleanup | `src/paths/long/large/Clean200.cc` | weak-branch error correction |
| Path fixup | `src/paths/long/large/GapToyTools.cc:699` | `FixPaths` |
| Manual | `docs/DiscovarExp` | user-facing description |

---

## 2. Stage A — the K=60 unipath graph (`buildReadQGraph`)

`src/paths/long/BuildReadQGraph.cc` (K hard-coded to 60 at line 51), called
from `GapToyCore.cc:799` with `MIN_QUAL=7`, `MIN_FREQ=3` (defaults from
`GapToyCore.cc:277-278`).

### 2.1 Kmer dictionary (`createDict`, lines 148–189)

1. `GoodLenTailFinder` (lines 59–83) finds, per read, the longest suffix in
   which every base has quality ≥ `MIN_QUAL`; only that suffix is kmerized.
2. A MapReduce `Kmerizer` (lines 100–146) streams all 60-mer *instances*
   through a disk-spilling sort (`MapReduceEngine`) and reduces them into a
   `KmerDict` keeping only kmers with count ≥ `MIN_FREQ`. A first MapReduce
   pass (lines 164–173) counts distinct qualifying kmers **to size the
   dictionary exactly**.
3. Each surviving kmer carries a `KMerContext` (`src/kmers/KMerContext.h`):
   8 bits — one predecessor bit and one successor bit per base A/C/G/T. This
   is the local branching structure of the de Bruijn graph at that kmer.
4. Because weak kmers were dropped, `recomputeAdjacencies()`
   (line 185–186) patches up contexts across the holes.

### 2.2 Unitig (edge) construction (`EdgeBuilder`, lines 191–401)

The classical **maximal non-branching path** algorithm, parallelized over
dictionary buckets:

- For each unassigned canonical kmer (`buildEdge`, lines 197–207):
  - palindrome → single-kmer edge (`make1KmerEdge`);
  - upstream extension possible (unique predecessor, and that predecessor has
    a unique successor) → extend upstream; likewise downstream; both possible
    → this kmer is interior to an edge someone else will build, skip;
  - neither → single-kmer edge.
- `extend()` (lines 307–326) walks while the current kmer has exactly one
  predecessor and one successor; palindromes and loss of uniqueness break
  the edge. Results whose canonical form is REV are discarded (built from
  the other strand).
- `addEdge` (lines 340–368) records `(edge id, offset)` back into every
  kmer's `KDef` (`src/kmers/ReadPather.h:104-145`), so each kmer knows its
  unitig and position; double-assignment is fatal (`ForceAssert`).
- Kmers left unassigned form **smooth circles**: `simpleCircle`
  (lines 210–234, single-threaded) walks the cycle and
  `canonicalizeCircle` (237–259) rotates it to start at the minimal kmer.

### 2.3 Optional gap fixes

- `fillGaps` (lines 765–847): re-examines captured gaps between path parts
  and inserts kmers spanning them (enabled by `FILL_JOIN=True`).
- `joinOverlaps` (lines 907–1017): stitches edge ends supported by ≥
  `minFreq2` observed overlaps (also `FILL_JOIN`).

Both rebuild the edge set afterward.

### 2.4 HyperBasevector assembly (`buildHBVFromEdges`)

`src/paths/long/HBVFromEdges.cc:243-281`:

- **Vertices** are the unique (K−1)-mer end sequences of unitigs
  (`VertexDictBuilder`, lines 135–167, a MapReduce over `EdgeEnd`s; the
  (K−1)-mer at each edge end identifies the vertex).
- `HBVBuilder` (lines 169–237) traverses connected components in a
  breadth-first queue, adding each unitig **in forward and in
  reverse-complement orientation** (`fwdXlat`/`revXlat` edge-id
  translations). The resulting `HyperBasevector` is symmetric: every edge
  has an RC twin and `hbv.Involution(inv)` (GapToyCore.cc:813) yields the
  reverse-complement map used everywhere downstream. Palindromic edges are
  self-RC.
- Edge processing order is canonicalized (length descending, then lexical,
  `EOComp` lines 122–133) so output is deterministic across runs.

### 2.5 Read placement (`HBVPather`, lines 1020–1284)

Each read is seeded on a trusted kmer and extended by exact (or
quality-weighted fuzzy, `path_more`) matching against the unitigs, producing
`PathPart`s — `(edge, rc, offset, length)` runs separated by gaps
(`Pather::path`, lines 562–604). `algorithmTwo` (the default "new aligner",
`GapToyCore.cc:806`) additionally drops seeds on short hanging ends, rejects
non-conforming captured gaps, and stretches placements using quals
(`ExtendReadPath`). `pathPartsToReadPath` (lines 1218–1253) converts to a
`ReadPath` (`src/paths/long/ReadPath.h`): the edge-id list plus the read's
offset on the first edge (negative = read extends past the edge start).
Unpathed reads may get a second chance via `ShortKmerReadPather`
(`SHORT_KMER_READ_PATHER=True`, GapToyCore.cc:327).

Back in `GapToyCore.cc:837-854`, `FixPaths` (`GapToyTools.cc:699-712`)
truncates any read path at its first adjacency violation, and optional
`JoinPaths`/`Preclose` variants can merge or pre-close paths.

---

## 3. Stage B — K=200 re-unipathing (`Repath`)

`Repath` (`src/paths/long/large/Repath.cc:23-26`), called at
`GapToyCore.cc:860`:

```cpp
Repath( hbv, edges, inv, paths, hbv.K() /* 60 */, 200,
        run_head+".200", /*REPATH_TRANSLATE=*/True, /*INVERT_PATHS=*/True,
        EXTEND_PATHS );
```

Inputs: the K=60 `HyperBasevector`, its edge sequences, the edge involution,
and the reads' paths on it. **The reads themselves are not an input.**

### 3.1 Places (lines 34–61)

"Places" are the distinct read paths on the K=60 graph:

1. each read's path is copied; paths implying < 200 bases of sequence are
   dropped;
2. the path is canonicalized to the lexicographically smaller of itself and
   its RC (`placesm.push_back(x<y?x:y)`, line 55);
3. all places are unique-sorted (`sortInPlaceParallel`, `Unique`).

Note what this dedup implies: a read that lies entirely within one K=60
edge contributes the single-element place `[e]` — **one place per edge**,
not one per read. Multi-edge places are distinct *edge lists* crossing
branch vertices. The place count is therefore governed by the K=60 graph
geometry (millions), not the read count (hundreds of millions).

`EXTEND_PATHS` (off by default, GapToyCore.cc:302) optionally extends each
place through unambiguous K=60 regions before dedup (lines 65–86).

### 3.2 Place sequences (lines 91–113)

Each place is converted to a basevector by concatenating its K=60 edges with
their (K−1)=59-base overlaps removed (`b.resize(b.isize() - (K-1)); b = Cat(...)`,
lines 100–102). For multi-edge places the **first and last edges are
truncated so each contributes at most 200 bases** (lines 103–112); the
truncation amounts (`left_trunc`, `right_trunc`) are retained for the offset
arithmetic in the translation step. Single-edge places contribute the whole
edge. The result is the `all` vecbvec, one sequence per place.

### 3.3 The K=200 unitig build (line 122)

```cpp
unsigned const COVERAGE = 2u;   // Repath.cc:121
LongReadsToPaths( all, K2=200, COVERAGE, /*logLevel=*/0,
                   /*useOldLRPMethod=*/False, &hb2, &h2, &xpaths );
```

`LongReadsToPaths` (`src/paths/long/LongReadsToPaths.cc:263`) dispatches to
`buildBigKHBVFromReads(200, all, coverage, ...)` (line 293) →
`src/kmers/BigKPather.cc:521-527`, which selects the compile-time
specialization for K=200 via `BigK::dispatch`
(`src/paths/long/LargeKDispatcher.h`: the `gK[]` table, lines 22–24; K must
be in that table).

`readsToHBV<200>` (`BigKPather.cc:449-508`) is the K=200 analogue of Stage A
(with no quality/frequency filtering — the input is already trusted):

1. **Kmerize** (`BigKMerizer<200>::kmerize`, lines 39–54): every place is
   walked with a rolling `KMerHasher<200>` FNV hash; each 200-mer is
   canonicalized (RC if needed) and inserted into a single parallel
   `HashSet` (`BigDict<200>` = `HashSet<BigKMer<200>>`), OR-merging its
   `KMerContext` at every occurrence.
   A `BigKMer` (`src/kmers/BigKMer.h`) does **not** store 200 bases: it is a
   24-byte view `{bvec const* mpBV, unsigned mOffset, bool mRC,
   bool mAssigned, KMerContext mContext, size_t mHash}` into the place
   sequence (bases resolved on demand through `operator[]`), plus the FNV
   hash (used for bucketing via `BigKMer::hasher`, lines 29–34).
2. **Unitigs** (`BigKEdgeBuilder<200>::buildEdges`, lines 110–137): the same
   maximal-non-branching-walk algorithm as §2.2 — `buildEdge`
   (143–158), `upstream/downstreamExtensionPossible` (214–236),
   `extend` (243–265), palindromes (206–212, 238–241), circles
   (`simpleCircle` 161–187 + `canonicalizeCircle` 189–204). The invariant
   `AssertEq(edges.getKmerCount(200), bigDict.size())` (line 482) checks
   that every distinct 200-mer ends up on exactly one edge.
3. **Dict re-pointing** (`updateDict`, lines 85–93, invoked at 493–495):
   every kmer entry is re-located to point into the *edge* sequence it
   occurs in, so placement is O(1) from kmer → (edge, offset).
4. **Graph assembly** — `buildHBVFromEdges(edges, 200, ...)`:
   `HBVFromEdges.cc:243` (as §2.4, vertices = shared 199-mers at edge ends,
   forward+RC edge insertion, canonical order), yielding `hb2`; then
   `buildHKPFromHBV` (`HBVFromEdges.cc:284-337`) builds the `HyperKmerPath`
   mirror `h2`, assigning each edge a consecutive run of kmer ids
   (palindrome ids from a reserved region).
5. **Re-pathing the places** (`Pather<200>`, lines 308–392): each place
   sequence is walked across the new K=200 edges via the dict, producing
   `xpaths` — one `KmerPath` (global-kmer-id intervals) per place — which
   is returned for the translation step. (The edge-id list
   `mTmpReadPath` is built internally but discarded; see §7 item M5.)

`all` is destroyed immediately after (`Repath.cc:123`).

### 3.4 Path translation ("very ugly", lines 139–255)

The original reads' K=60 paths must be restated on the K=200 graph. Two
lookup structures are built:

- `hpaths`: the `KmerPath` of each K=200 edge (`h2.EdgeObject(e)`);
- `hpathsdb`: a kmer-interval database over `hpaths` (`CreateDatabase`,
  line 147).

For each place, its `xpaths` kmer path is matched interval-by-interval
against `hpathsdb` (`Contains`, lines 158–176) to recover
`ipaths2[place]` = the K=200 edge sequence the place traverses, plus
`starts[place]` / `stops[place]` = how the place's offset maps into the
first/last K=200 edge (lines 178–209; gaps in coverage or partial matches
mark the place `bad` and it is skipped).

Then each original read (lines 216–248): recompute its canonical path,
binary-search it in `places` (`BinPosition`, line 233), copy
`ipaths2[place]` (reversed + RC'd via `inv2` if the read was canonicalized),
and set the offset from `paths[id].getOffset()` plus `starts/stops` minus
the §3.2 truncations. Results are written as `a.200.paths`; with
`INVERT_PATHS`, the edge→reads index `a.200.paths.inv` is also built.

Known bug (cosmetic today, fatal if touched carelessly): the path-length
filter indexes `edges` by path *position* instead of edge id —
`nkmers += edges[j].isize( ) - ( (int) K - 1 );` at **Repath.cc:51 and
:227** — should read `edges[ x[j] ]`. It only mis-sizes the "too short to
contribute" filter.

### 3.5 Outputs

Written by Repath at `run_head+".200"`: `.hbv` (HyperBasevector, K=200),
`.hbx` (compact variant), `.fastb` (edge sequences), `.inv` (edge
involution), `.paths`, `.paths.inv`. `GapToyCore.cc:902-906` relabels them
into `a.200/` with `a.k = 200`.

---

## 4. Stage C — cleanup and artifacts

`Clean200` (`src/paths/long/large/Clean200.cc:20-`, selected via
`CLEAN_200V`, GapToyCore.cc:876-881; version 3 uses `Clean200x`) runs **two
passes** (line 31) of branch-resolving error correction using the reads:

for each branch vertex it enumerates alternative extensions up to read
length (250), converts them to basevectors, aligns the read paths against
each alternative (offset-adjusted, quality-weighted), and deletes
weakly-supported edges. `REMOVE_TINY` optionally drops tiny standalone
edges. This is where the reads (`bases`, `quals` — reloaded at
GapToyCore.cc:877) are needed again.

Artifacts of `a.200/` after this stage: `a.hbv`, `a.hbx`, `a.fastb`,
`a.inv`, `a.paths`, `a.paths.inv`, `a.k`.

Everything downstream — local gap assembly (`AssembleGaps`), patching
(`a.patched`), simplification, lines (`a.fin`), scaffolding
(`MakeGaps` → `a.final`) — operates on this K=200 graph. Per
`docs/DiscovarExp`, circular sequences in the final output overlap by
exactly K−1 = 199 bases.

---

## 5. Key data structures and invariants

| Structure | Definition | Notes |
|---|---|---|
| `HyperBasevector` | `src/paths/HyperBasevector.h` | directed graph, edges are basevectors, `K()` = 200 |
| `HyperKmerPath` | `src/paths/HyperKmerPath.h` | same graph, edges are consecutive kmer-id runs (`buildHKPFromHBV`) |
| `ReadPath` / `ReadPathVec` | `src/paths/long/ReadPath.h` | edge-id list + signed offset on first edge; feudal, pooled storage |
| `KMerContext` | `src/kmers/KMerContext.h` | 8-bit pred/succ branching record |
| `BigKMer<200>` | `src/kmers/BigKMer.h` | 24-byte (sequence view, context, hash); equality = hash match **then** base-compare (`operator==`, lines 111–114) |
| `KmerDict`/`KDef` | `src/kmers/ReadPather.h` | K=60 dict; entry records (edge, offset) once assigned |
| `vecbvec` | `src/Basevector.h:19-21` | `MasterVec<BaseVec>` — feudal, pooled element storage |
| `vecKmerPath` | `src/paths/KmerPath.h:325` | `MasterVec<KmerPath>` — feudal, pooled |
| `vec<T>` | `src/Vec.h` | wraps `std::vector` — **no pooling** (relevant for `vec<vec<int>>` below) |

Invariants worth preserving when modifying the code:

1. **Involution symmetry** — both graph builders insert every unitig in
   forward and RC orientation; downstream code (translation, scaffolding,
   `TestInvolution`) assumes an edge-level RC map exists.
2. **Kmer→edge totality** — every distinct kmer is assigned to exactly one
   edge (asserted at `BigKPather.cc:482` and enforced by `ForceAssert`s in
   `addEdge`).
3. **Determinism** — edges and vertices are processed in canonical order
   (`EOComp` in `HBVFromEdges.cc`, tag-sort in `LongReadsToPaths.cc`);
   keep any new output canonicalized if bit-stable results matter.
4. **K dispatch table** — changing K requires editing the `gK[]` table
   *and* the matching `case` list in `LargeKDispatcher.h` (see warnings
   there), then recompiling; the code is fully templated per K.
5. **The K=200 `HyperBasevector` must have `K()==200`** (e.g. asserted in
   `FixInversion`, `GapToyTools.cc:750-751`).

---

## 6. Memory profile of the a.200 construction

Numbers below are **estimates for a diploid human at ~60× coverage with
250-base reads** (trusted sequence ≈ 3.3–3.5 Gbp, ~720M reads, ~10–40M
places); exact values depend on K=60 graph geometry. The pipeline itself
warns at ~2 bytes per read base (`GapToyCore.cc:772-787`).

### 6.1 Verified allocation mechanics

1. **The 200-mer dictionary is sized by *instances*, not distinct kmers.**
   `Repath.cc:121` hard-codes `COVERAGE = 2`, and
   `BigKPather.cc:473` allocates `BigDict(nKmers/coverage)` where
   `nKmers = reads.getKmerCount(200)` counts every 200-mer **instance**
   across all places (`OuterVec.h:188-193`). Places overlap heavily
   (neighboring fragments traverse the same K=60 edges), so instances ≫
   distinct — only ~trusted-genome-sized distinct 200-mers (~3.3–3.5 G)
   are ever inserted.
2. **The `HashSet` eagerly allocates the entire bucket array at
   construction.** `HashSet::init` (`src/feudal/HashSet.h:851-864`) creates
   *every* inner `HopscotchHashSet` upfront, each with its full bucket array
   (`HopscotchHashSet` ctor, lines 184–189: `new T[capacity]` + 1 byte of
   `BktStatus` per slot). With `sizeof(BigKMer<200>) = 24`, a slot costs
   **25 bytes**. Reserved ≈ `(nKmers/2)/0.75` slots ≈ **16.7 bytes per
   kmer instance**, versus the **~33 bytes per distinct kmer**
   (25 B/slot at 0.75 load) actually required. At typical place redundancy
   (3–10×) the dictionary is 1.5–5× oversized — pure, immediate waste.
   (Under-sizing is *safe* — the set grows by splitting, `HashSet.h:906-937`
   — at the cost of rehash CPU and a warning at teardown,
   `HashSet.h:953-954`.)
3. **The reads stay resident although `Repath` never touches them.**
   `Repath`'s inputs are `hbv, edges, inv, paths` only; `bases` (~45 GB) is
   needed again just for `Clean200` (GapToyCore.cc:877). The code already
   unloads/reloads the quality vectors around exactly this window
   (`quals.unload()` at line 859, `quals.load()` at 877) — but not the
   bases.
4. **Feudal vs. plain containers.** `vecbvec`/`ReadPathVec`/`vecKmerPath`
   are `MasterVec`s with pooled element storage (no per-element mallocs).
   Plain `vec<vec<int>>` (used for `places`, `ipaths2` in Repath.cc) wraps
   `std::vector` (`Vec.h`) — each inner vector is a separate heap
   allocation with ~40–60 B effective overhead. With 10–40M places this is
   ~1–3 GB (not dominant, but also the data is scattered, which hurts the
   parallel sort/unique).
5. **The spill-sort buffer pool is sized by the machine, not the work.**
   `MapReduceEngine::run` (`MapReduceEngine.h:403-424`) budgets its
   per-thread key buffers at up to 90% of the *allowed* memory
   (`MemAvailable(.9)`, i.e. of `MAX_MEM_GB`) and then scales batches up
   by another ×2 — regardless of how small the input is. This dominated
   small-genome runs entirely: the E. coli 50× test peaked at 6.8 GB,
   nearly all of it kmer-instance buffers sized from the 8 GB allowance.
   (Addressed by M8 below.)

### 6.2 Approximate peak budget (stage B, `Repath` window)

| Component | Approx. size | Notes |
|---|---|---|
| `BigDict<200>` | **60–375 GB** | 25 B/slot × (instances/1.5); the dominant term, see §6.1(1)(2) |
| `bases` (reads) | ~45 GB | unused during Repath; reloadable |
| K=60 `paths` (ReadPathVec) | ~15–20 GB | needed (places + translation) |
| `all` (place sequences) | ~2–8 GB | single-edge places = full K=60 edges once + boundary paths |
| `xpaths` + `hpathsdb` + `ipaths2` | ~2–5 GB | translation machinery |
| `hbv60` + K=60 `edges` copy | ~1.5–3 GB | `GapToyCore.cc:811` duplicates `hbv.Edges()` |
| K=200 edges/HBV/HKP | ~2–3 GB | stored fwd+RC by design |

The dictionary and the reads dominate; everything else is ≤ ~30 GB
combined.

---

## 7. Memory improvement plan (ranked)

Effort ⓔ easy / ⓜ moderate / ⓗ hard. Savings are for the §6.2 assumptions.

### M1 — Right-size the 200-mer dictionary ⓔ (tens to ~250 GB) — **IMPLEMENTED**

Size the hash set by an estimate of **distinct** 200-mers instead of
`instances/2`. The best estimate is already available at the call site: the
K=60 `edges` sum to ≈ the trusted sequence, i.e.
`edges.SizeSum() ≈ distinct 200-mer count` (pad ×1.1–1.2, since under-
sizing costs a split/rehash pass while over-sizing costs RAM linearly —
§6.1(2)).

Concrete edits:

- `Repath.cc:121`: replace the fixed `COVERAGE = 2u` with
  `coverage = nKmers / max(1, distinctEstimate)` — or, cleaner, extend the
  signatures `LongReadsToPaths` (`LongReadsToPaths.h:25-32`) /
  `buildBigKHBVFromReads` (`BigKPather.h:29-33`) / `readsToHBV`
  (`BigKPather.cc:450`) to accept an explicit `dictSize` and use it at
  `BigKPather.cc:473`.
- Optionally pass a higher `maxLoadFactor` (5th ctor arg of `HashSet`,
  `HashSet.h:626-632`, default 0.75; hopscotch tolerates ~0.9) → another
  ~15% fewer buckets.

Implementation notes (as landed): `Repath` computes the estimate as
`Σ over K=60 edges (len−K2+1)`, counting each edge/RC-twin **involution
pair once** via `inv` (the `edges` array holds both orientations, so a
naive sum would double the estimate), plus a 25% junction margin; a zero
estimate falls back to the legacy `COVERAGE=2` sizing. The hint is threaded
as a trailing `dictSize=0` parameter through `LongReadsToPaths`,
`buildBigKHBVFromReads` and `readsToHBV` (all other callers are
unaffected), and capped at the total kmer-instance count. The legacy path
converts a hint into the equivalent coverage value.

Expected: dictionary drops to ~110–120 GB (25 B/slot at 0.75 LF,
~3.5 G distinct) regardless of place redundancy.

### M2 — Unload the reads during `Repath` ⓔ (~45 GB) — **IMPLEMENTED**

Mirror the existing `quals.unload()/load()` pattern for `bases` around
`GapToyCore.cc:859-861` and `:874`.

Implementation notes (as landed): instead of writing a dedicated copy, the
reads are re-read from `data/frag_reads_orig.fastb` — the same file the
`CACHE=True` path already reloads from (`GapToyCore.cc:913`), which
establishes it as the canonical on-disk copy of the selected reads (every
`ExtractReads` flow, including subsampling, rewrites it before returning).
`bases` is freed with a swap against an empty `vecbvec` before `Repath`
and reloaded immediately after, before `Clean200`. No extra disk write.

```cpp
quals.unload();
{    vecbvec empty; bases.swap(empty);    }
Repath( ... );
bases.ReadAll(tmp_dir1 + "/frag_reads_orig.fastb");
```

### M3 — Shrink `BigKMer<200>` from 24 to 16 bytes ⓜ (−35% of the dict)

Current layout (`BigKMer.h:133-138`): pointer(8) + offset(4) + 2 bools +
context(1) + **64-bit hash**(8) → 24 B padded. Proposed:

- `{uint32 placeId, uint32 offset}` (8 B) indexing a small segment table —
  replaces the raw `bvec*` (place count < 2³²; offsets can exceed 2²⁰ for
  single-edge places, so keep 32 bits);
- 32-bit stored hash (4 B);
- fold `mRC`/`mAssigned` into spare bits; `KMerContext` (1 B).

→ 16 B padded ⇒ 17 B/slot: dictionary ~110 GB → ~70–75 GB.

**Caveat (verified):** `operator==` (`BigKMer.h:111-114`) gates the
base-by-base comparison behind a full 64-bit hash match. With a 32-bit
stored hash, equality must become "stored-hash match, then *always*
base-compare"; otherwise distinct kmers sharing a 32-bit hash would
falsely match. Bucketing uses only the stored hash (`BigKMer::hasher`), so
32-bit hashes are fine there; collisions just cost a 200-base walk. All
base access already funnels through `BigKMer::operator[]`/iterators, so the
change is contained in `BigKMer.h` + `BigKPather.cc`.

### M4 — Virtualize place sequences ⓜ (~2–8 GB)

Places are concatenations of K=60 edges (already in memory as `Repath`'s
`edges` argument). With M3's segment-table indirection in place, a natural
follow-on is to make `BigKMerizer` walk places **virtually** (resolve bases
on demand from the K=60 edge sequences) instead of materializing `all`
(`Repath.cc:91-113`): eliminates the `all` copy and its co-residency with
the dictionary. Do M3 first; this shares its machinery.

### M5 — Emit K=200 edge lists directly; drop the "very ugly" translation ⓜ (few GB + CPU + removes a bug)

`Pather<200>` (`BigKPather.cc:322-352`) already builds the K=200 edge-id
list (`mTmpReadPath`) and the offset/stop arithmetic per place, then
discards the edge list and keeps only kmer-id paths (`xpaths`). Meanwhile
`Repath.cc:139-209` re-derives those same edge lists via
`xpaths`/`hpaths`/`hpathsdb` interval matching. Extending `Pather` to also
emit `(edge-id list, start, stop)` per place lets `Repath`'s translation
collapse to: re-path each original read through the place index
(`BinPosition`, `Repath.cc:233`) + RC-reversal via `inv2`. This deletes
`xpaths`, `hpaths`, `hpathsdb` and ~70 lines of interval algebra, removes
the code containing the `edges[j]`-vs-`edges[x[j]]` indexing bug
(§3.4), and speeds the stage up.

While in there, consider flattening `places`/`ipaths2` from
`vec<vec<int>>` (plain `std::vector`s per `Vec.h`) to CSR
(flat `vec<int>` + offsets — the pattern `VecULongVec` already uses):
saves the per-element heap overhead and makes `sortInPlaceParallel`/
`Unique` cache-friendly.

### M6 — Minor items ⓔ

- `vecbvec edges( hbv.Edges()... )` at `GapToyCore.cc:811` duplicates the
  K=60 edges (~0.75–1.5 GB); pass `hbv` and index `hbv.Edges()` where
  possible.
- The `paths_index` inversion is built twice at ~10 GB each: inside
  `Clean200` (`Clean200.cc:45-46`) and for `a.200.paths.inv`
  (`GapToyCore.cc:893-897`); build once and reuse.
- `MALLOC_PER_THREAD=1` (see `docs/DiscovarExp`) reduces allocator arena
  fragmentation with these huge transient allocations.

### M7 — Disk-backed kmerization ⓗ (last resort, for very large genomes)

Trade the ~70–120 GB dictionary floor for disk I/O: build the 200-mer dict
via the existing spill-sort `MapReduceEngine`, exactly as the K=60
`Kmerizer` does (`BuildReadQGraph.cc:146-189`) — map compact 16–20 B
instance records (hash + placeId/offset + context) to sorted runs, reduce
to distinct kmers, then insert only distinct entries into the in-RAM dict.
Unnecessary once M1–M3 are in place.

### M8 — Cap the MapReduce buffer pool by need ⓔ (small-input sanity) — **IMPLEMENTED**

`MapReduceEngine::run` sized its key-buffer pool from the memory
allowance (90% of `MAX_MEM_GB`, then ×2), not from the work: a bacterial
genome on a big machine grabbed the whole allowance. The engine gains an
optional `maxMemBytes` budget on `run()`, and a helper,

```cpp
mrePoolBytes( nKeys, keySize, divisor = 8., minBytes = 512 MB )
```

sizes it by the work: one eighth of the raw key-instance stream, with a
512 MB floor, never exceeding the old allowance. Small inputs then run
in a few extra passes over the (already in-memory) input — cheap — while
large inputs keep today's behavior because the work-based estimate
exceeds the allowance and the old cap wins.

Call sites updated: the two `createDict` kmerize runs
(`BuildReadQGraph.cc:170/179`) and the vertex-dict run
(`HBVFromEdges.cc:261`, both K=60 and K=200 stages).

Measured on the E. coli 50× test: **peak 6.8 GB → 1.1 GB**, runtime
52 s → 65 s, and a bit-identical result (same final checksum), confirming
that the pool size affects buffering only, not the dictionary or the
assembly.

### M9 — hash-partitioned disk bins for the K=60 kmerize ⓔⓜ — **IMPLEMENTED**

The KMC3/DSK lesson: the kmerize never needed a giant in-RAM sort pool —
it needed *partitioning*. `createDict` now streams each canonical
k-mer instance to one of N **disk bins** (record = raw kmer words + 1
context byte, 17 B at K=60; bin = a cheap word-hash of the canonical
kmer % N). Per-thread 16 KB append buffers with lock-free `pwrite` at
per-thread offsets keep the writers parallel; N is chosen for ~256 MB
bins, capped at 992 (file-descriptor limit). Bins are then processed a
few at a time (concurrency capped by a memory budget): each bin is
read, sorted in RAM as packed records, and reduced with the exact
`summarizeEntries` semantics (group size = count, OR of context bytes)
— qualifying kmers are streamed to a compact 25 B *qualified* record
file (words, context, count, bucket id). The dictionary is then built
directly from that stream (below). This eliminates the allowance-bound
MapReduce pool at every scale, at the cost of ~17 B/instance of
temporary disk (≈ the read set in 2-bit + slack), and *no* Bloom filter
is needed (the error-k-mer explosion makes Bloom sizing
O(distinct-including-errors) ≈ 30–40× the real kmer space — see the
exploration notes in the repo history).

As landed (`createDict` in `BuildReadQGraph.cc`): pass 1 kmerizes reads
in ~10-per-thread work chunks, each with chunk-local bin buffers and
atomically reserved `pwrite` byte ranges per bin; pass 2 sorts/reduces
bins under a concurrency cap (`MemAvailable(.25)`-based); qualifying
kmers stream to a compact record file; a sequential sweep places them
into the M10 arena by bucket.  The old `Kmerizer`/`MapReduceEngine`
counting passes are gone (legacy `rePath` still uses the feudal
`KmerDict` via `DictionaryKmerEater`).  Two progress lines report the
instance count and the qualifying kmer count.  Engineering notes from
bring-up: per-chunk file offsets must be *shared per bin* (chunks
overwrite each other otherwise — lock-free `fetch_add` reservation);
`omp_get_thread_num` is invalid inside `parallelFor` (pthreads-based
`WorklistN`), so all per-thread state is chunk-local.

### M10 — bucketed sorted-array K=60 dictionary ⓔⓜ — **IMPLEMENTED**

The qualifying stream from M9 is already deduplicated, so the feudal
`HashSet` backing of the K=60 `KmerDict` is replaced (for the
DiscovarExp path) by `KmerDictPacked<K>` (`kmers/KmerDictPacked.h`): a
contiguous arena of 24 B `KmerDictEntry`s, bucketed by the same
word-hash into ~4 K entries per bucket (directory = prefix offsets),
sorted by kmer within each bucket, with the full `KmerDict` method
surface ported (`findEntry`/`findEntryCanonical` = bucket + binary
search over an ~80 KB range; `apply`/`applyCanonical`; bulk + legacy
`insertEntry`; `parallelForEachHHS` over disjoint array chunks;
`nullEntries`; `recomputeAdjacencies`; HHS-style iteration; a
single-threaded legacy `process`). Savings versus the HashSet:
no 1 B/slot status, no 0.75 load-factor slack, no split transients
(a hash split briefly doubles its table), and *exact* sizing from the
qualifying count — ~24 B/kmer contiguous versus ~33 B effective, i.e.
~1.4× on the largest stage-A structure, plus better cache behavior for
the EdgeBuilder walk. Legacy code keeps the old `KmerDict` (`rePath`
pins to it explicitly). The K=200 dictionary is *not* converted here —
it has no sorted-input stream yet; that becomes **M10b**, naturally
built on a stage-B binning pass.

**Caveat (acceptance):** dictionary iteration order changes, so edge
id assignment is relabeled and the `final checksum` is *expected* to
change. The equality gate is: identical stats block (edges/kb, N50s,
recovery, perfect stretch, misassemblies), identical `a.fasta`
edge-sequence multiset, and bit-identical output across reruns.

As landed, on the E. coli 50× test: **stats block byte-identical,
edge-sequence multiset fingerprint identical, and the final checksum
itself identical (617468569301748) across 18 consecutive runs** — the
edge walk reproduces the same assignment order in practice.  Peak
memory 1.2–1.3 GB (vs 1.1 GB under M8; the arena plus bin sorts
displace the pool at this scale — the structural win is at genome
scale, where the allowance-bound pool disappears entirely).  Bring-up
notes: the context byte round-trip must preserve the KMerContext
nibble layout (preds high, succs low — a swapped decode produces
overlapping edge walks, caught by `addEdge`'s "preoccupied kmers"
check); the whole-dict iterator must not re-enter `settle()` after
arena exhaustion (it parks back at the arena start, looping forever).

**M10b — compact K=200 dictionary** (`BigKMer`/`BigKPather`): the
K=200 entries shrink from 24 B pointer-views to 16 B index-views —
`{seqIdx:u32, offset:u32, hash32:u32, rc+assigned+context}` — resolved
through a **thread-local sequence registry** (places first, edges after
the edge set is final) set on entry to every parallel phase, with
rolling-scratch temporaries (`mEdgeSeq`) resolving through a
thread-local scratch pointer.  A truncated 32-bit FNV is kept as a fast
inequality filter and doubles as the in-bucket discriminator; the
hopscotch load factor rises to 0.9.  Disk binning was *rejected* for
stage B: there is no count threshold to evaluate (the input is already
trusted), and instance materialization would add a transient the
hash-OR insert does not have.  Human-scale dictionary: ~33 B/entry
effective → ~19 B (≈117 GB → ~67 GB).  Critical bring-up lesson: the
local gap assemblies run several `readsToHV` calls concurrently inside
an OpenMP region, so a *static* registry is clobbered across blobs —
thread-local storage is mandatory.  Acceptance: same gates as above,
all identical, 10 consecutive runs after the fix below.

**Pre-existing feudal bug found by the M10b hunt (fixed):**
`PQVecEncoder::decode` (`feudal/PQVec.cc`) reads 64-bit words aligned
down from each quality vector's start, so decoding the last vector in
an exactly-sized `Mempool` pre-allocation reads up to 7 bytes past
the pool end — an AddressSanitizer-confirmed over-read that had been
silently absorbed by allocator slack since 2010, occasionally
segfaulting when the pool abuts an unmapped page (the ~7% flake seen
during M10b bring-up).  `Mempool::preAllocate` now keeps 8 bytes of
unbooked tail slack.  After the fix: a full ASan run is clean
end-to-end, and the E. coli battery runs 10/10 clean.

### M11 — minimizer-sparse dictionaries ⓗ (planned)

Store only anchors: the 1-in-w minimizer subset of kmers (plus unitig
boundary handles), à la sparse de Bruijn graphs (Chikhi et al.) and
BLight. Theory says the graph is *identical* if the scheme is
deterministic, strand-symmetric, and handles are kept — and Discovar
has two lucky properties: K=60 and K=200 are both even (no palindromic
kmers), and canonical minimizers preserve the involution. Persistent
dictionaries drop ~10–20× (anchors-as-views for K=200 ≈ 16 B each).
The one honest semantic change is at K=60: trust becomes per-anchor
rather than per-kmer (interior counts never individually observed);
practically negligible, but it must be gated by the M10 acceptance
procedure on real data. Requires a bounded-density minimizer scheme
(bd-anchors/mod-minimizers, not naive lexicographic) to avoid anchor
deserts in tandem repeats.

### M12 — retire post-construction dictionaries via an edge-minimizer
index ⓗ (planned)

`KDef` exists only to path reads; the Pather is already seed-then-extend.
Seeding on minimizers *of edge sequences* (minimap2-style) needs only
the trusted sequence plus a ~12 B/anchor index — no per-kmer structure
at all after construction. Human-scale peak plausibly lands at
~15–25 GB, and the same index is exactly the LR-1 substrate the ONT
plan needs (projection of long reads onto the graph) — dual-use.

### Expected result and order of attack

| Layer | Stage A transient | Persistent dicts | Peak estimate |
|---|---|---|---|
| today (M1/2/8) | allowance-bound pool | 113 + 115 GB | ~150–160 GB |
| + M8 (pool by need) | same at human scale; 1.1 GB on E. coli | 113 + 115 GB | ~150–160 GB |
| + M9 (disk bins) | ~1–2 GB + disk ≈ readset×2.2 B (E. coli: 1.2–1.3 GB peak) | 113 + 115 GB | ~200 GB (dict-bound) |
| + M10 (bucketed array, K=60) | — | ~82 + 115 GB (E. coli dict: 4.6M×24 B = 110 MB) | ~160 GB |
| + M10b (K=200, needs stage-B bins) | — | ~82 + ~58 GB | ~120 GB |
| + M11 (minimizer-sparse) | — | ~6 + ~7 GB | ~35 GB (paths60 co-dominant → flatten) |
| + M12 (edge-minimizer pathing) | — | none | **~15–25 GB** |

E. coli sanity at every layer: trivial memory; its role is the
acceptance gate (stats block + edge-sequence multiset + rerun
determinism).

Order: **M2 → M1 → M8 → M9 → M10 → [M10b] → M11 → M12.**
M2/M1/M8/M9/M10 are landed; M11+ are planned.

---

## 8. Profiling and verification guidance

- Existing checkpoints: `MemUsageBytes()` is printed after
  `buildReadQGraph` (`GapToyCore.cc:808-810`) and at `Clean200` start
  (`Clean200.cc:51`); stage timings appear as "… used in ReadQGrapher"
  (`GapToyCore.cc:907`) and in `clock.log` (`LogTime`, GapToyCore.cc:453).
  When touching `Repath`/`BigKPather`, add `MemUsageBytes()` prints
  around: places construction, `all` construction, `LongReadsToPaths`,
  translation — that is where §6.2 can be confirmed or corrected.
- Cheap end-to-end runs: `X=tiny` / `X=small` / `X=medium`
  (`GapToyCore.cc:509-607`) exercise the full stage on fosmid-scale data;
  `EXIT=PATHS` (`GapToyCore.cc:849-850`) stops right after the K=60 stage;
  `CACHE=True` (GapToyCore.cc:196) resumes from prebuilt files — useful to
  isolate stage B when iterating.
- Graph health checks: `hb.CheckSum()` is printed for K=60 and K=200
  (`GapToyCore.cc:856-857, 952-953`); `TestInvolution` guards the edge RC
  map; `a.200` edge/N50 statistics are the first thing to compare before
  and after any change (see the historical tables in the
  `GapToyCore.cc` header comment for reference values).
- Feudal containers report utilization through `UtilizationReporter`
  (`src/feudal/Generic.h:17-37`; invoked from `OuterVec.h:115`,
  `SmallVec.h:87`, `FieldVec.h:256`) — a ready-made hook for tracking
  container waste.

---

*Document based on discovarexp-51885 (Broad Institute DiscovarDenovo
experimental). Line numbers refer to this exact source tree.*
