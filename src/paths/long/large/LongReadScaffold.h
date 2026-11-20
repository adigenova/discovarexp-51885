///////////////////////////////////////////////////////////////////////////////

#ifndef LONG_READ_SCAFFOLD_H_
#define LONG_READ_SCAFFOLD_H_

#include "Basevector.h"
#include "CoreTools.h"
#include "paths/HyperBasevector.h"

// P3/LR-2 (docs/Modernization.md): scaffold the K=200 assembly graph
// using long-read (ONT/PacBio) evidence, following the Liger architecture
// (transitive reduction + maximum weight matching) adapted for the
// Discovar doubled de Bruijn graph with diploid variant awareness.
//
// Pipeline:
// 1. Classify edges: repeat (coverage > threshold), short (< min_line),
//    normal.  Variant edges from het bubbles have ~0.5x coverage and are
//    NOT flagged as repeats (diploid-aware).
// 2. Seed ONT reads on edge sequences (exact k-mer matches).
// 3. Chain anchors into blocks; block transitions across non-adjacent
//    edges produce scaffold links.  Links involving repeats or short
//    edges are discarded.
// 4. Transitive reduction: for each link, search the HBV for a path
//    connecting the same pair (BFS up to max_bfs_depth).  If found,
//    the link is redundant (already explained by the graph) -> drop.
// 5. Maximum weight matching (greedy): sort links by support
//    descending; greedily accept links whose both endpoints are
//    unused.  Each edge-end participates in at most one link.
// 6. Build scaffold paths: follow HBV adjacency within components,
//    scaffold links across gaps.
// 7. Output: scaffolded FASTA, link stats, N50.

struct ScaffoldLink
{
     int left_edge;
     int right_edge;
     int gap;
     int n_support;
};

void ScaffoldWithLongReads( const HyperBasevector& hb,
     const vec<int>& inv, const vecbasevector& edges,
     const vec<double>& edge_cov,   // per-edge coverage from a.covs
     const String& reads_file,
     const unsigned K, const int min_read, const int min_support,
     const int anchor_step, const double repeat_factor,
     const int min_line, const int max_bfs_depth,
     const int max_mult,          // max k-mer occurrences used as anchors
     const int min_anchor_block,  // min anchors on each flanking block
     const double min_chain_cov,  // min anchored span / read length
     const String& mapper,        // "seed" (internal) or "mm2" (minimap2)
     const String& mm2_bin,       // path to the bundled minimap2 binary
     const int min_map_len,       // mm2: min aligned length per mapping
     const double min_map_id,     // mm2: min mapping identity
     const int min_seg,           // mm2: min per-edge segment length
     const String& debug_read,    // mm2: dump link storage for this read
     const String& out_head, const Bool VERBOSE );

#endif /* LONG_READ_SCAFFOLD_H_ */
