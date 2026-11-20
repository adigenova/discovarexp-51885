///////////////////////////////////////////////////////////////////////////////

#ifndef LONG_READ_PROJECTION_H_
#define LONG_READ_PROJECTION_H_

#include "Basevector.h"
#include "CoreTools.h"
#include "paths/HyperBasevector.h"

// P2/LR-1 (docs/Modernization.md): project long reads (ONT/PacBio) onto
// the K=200 assembly graph.
//
// The graph edges are error-corrected, so exact k-mer matches between a
// long read and the edge sequences are dense even at ONT error rates (a
// 19-mer survives intact with probability ~0.7 at 2% per-base error).
// We build a sorted k-mer index over the (doubled) edge sequences, seed
// each read by exact matches, and chain collinear anchors.  Because the
// HyperBasevector stores both orientations of every edge, one query per
// read k-mer covers both strands -- a read matching the reverse strand
// hits the mirror edges directly, and a consistent chain never mixes
// worlds (fwd and mirror edges live in separate graph components).
//
// Consecutive anchors are consistent if they are on the same edge with
// linear coordinates, or on graph-adjacent edges (to_right[e1] ==
// to_left[e2]) with coordinates matching the de Bruijn junction algebra
// (edges overlap by K-1 bases at their shared vertex).
//
// Acceptance (P2): >= 95% of reads >= MIN_READ project onto >= 2 distinct
// edges with consistent coordinates.

// One projected segment of a read: read bases [r1,r2) match edge bases
// [o1,o2) of edge e.
struct LRSeg
{
     int edge;
     int o1, o2;
     int r1, r2;
};

// Per-read result: the longest consistent chain, plus diagnostics.
struct LRProj
{
     int64_t id;
     vec<LRSeg> chain;
     int n_anchors;
     int n_breaks;
     int64_t read_len;
};
TRIVIALLY_SERIALIZABLE(LRProj);

void ProjectLongReads( const HyperBasevector& hb,
     const vecbasevector& edges, const String& reads_file,
     const unsigned K, const int min_read, const int anchor_step,
     const int max_anchors, const String& out_head, const Bool VERBOSE );

#endif /* LONG_READ_PROJECTION_H_ */
