///////////////////////////////////////////////////////////////////////////////

// ScaffoldLongReads: P3/LR-2 tool (docs/Modernization.md).  Scaffold
// the K=200 assembly graph using long-read (ONT/PacBio) evidence,
// following the Liger architecture with diploid awareness.

#include "MainTools.h"
#include "paths/long/large/LongReadScaffold.h"
#include "paths/long/large/Lines.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(HBV, "",
          "a.hbv of the assembly graph (e.g. OUT/a.fin/a.hbv)");
     CommandArgument_String_OrDefault_Doc(FASTB, "",
          "a.fastb of the same graph");
     CommandArgument_String_OrDefault_Doc(INV, "",
          "a.inv of the same graph");
     CommandArgument_String_OrDefault_Doc(COVS, "",
          "a.covs of the same graph (per-edge coverage)");
     CommandArgument_String_OrDefault_Doc(READS, "",
          "long reads: fasta or fastq, possibly gzipped");
     CommandArgument_String_OrDefault_Doc(OUT_HEAD, "scaffold",
          "prefix for output files");
     CommandArgument_UnsignedInt_OrDefault_Doc(K, 19,
          "seed length for exact k-mer anchoring");
     CommandArgument_Int_OrDefault_Doc(MIN_READ, 20000,
          "minimum read length for scaffold evidence");
     CommandArgument_Int_OrDefault_Doc(MIN_SUPPORT, 2,
          "minimum reads supporting a scaffold link");
     CommandArgument_Int_OrDefault_Doc(ANCHOR_STEP, 64,
          "approximate spacing between anchors along a read");
     CommandArgument_Double_OrDefault_Doc(REPEAT_FACTOR, 1.5,
          "edge with coverage > REPEAT_FACTOR * median is a repeat");
     CommandArgument_Int_OrDefault_Doc(MIN_LINE, 2000,
          "edges shorter than this are excluded from scaffolding");
     CommandArgument_Int_OrDefault_Doc(MAX_BFS_DEPTH, 2,
          "reserved (transitive-reduction BFS of the old link-based"
          " design; read threading no longer uses it)");
     CommandArgument_Int_OrDefault_Doc(MAX_MULT, 1,
          "seed k-mers occurring more than this many times in the graph "
          "are skipped (ambiguous anchors: repeat copies, diploid "
          "shared k-mers; reads then thread their own haplotype)");
     CommandArgument_Int_OrDefault_Doc(MIN_ANCHOR_BLOCK, 2,
          "min anchors on each flanking block for a scaffold link "
          "(single-anchor blocks are unreliable)");
     CommandArgument_Double_OrDefault_Doc(MIN_CHAIN_COV, 0.5,
          "min anchored span / read length for a read to count as "
          "high-quality scaffold evidence");
     CommandArgument_String_OrDefault_Doc(MAPPER, "mm2",
          "read-to-graph mapper: 'mm2' (bundled minimap2, one-strand "
          "reference, PAF paths in fwd/rev) or 'seed' (internal "
          "k-mer aligner)");
     CommandArgument_String_OrDefault_Doc(MINIMAP2,
          "minimap2/minimap2/minimap2",
          "path to the bundled minimap2 binary (MAPPER=mm2)");
     CommandArgument_Int_OrDefault_Doc(MIN_MAP_LEN, 300,
          "mm2: minimum aligned length per mapping");
     CommandArgument_Double_OrDefault_Doc(MIN_MAP_ID, 0.70,
          "mm2: minimum mapping identity (1 - NM/blen)");
     CommandArgument_String_OrDefault_Doc(DEBUG_READ, "",
          "mm2 (with VERBOSE=True): dump the raw PAF, the per-edge"
          " path segments, and every stored link record for this"
          " read");
     CommandArgument_String_OrDefault_Doc(K60_FASTB, "",
          "optional (LoRDEC-style K=60 polish; default empty = ONT"
          " MSA only).  a.fastb of the K=60 base graph (SAVE_60="
          "True writes OUT_DIR/a.60/a.fastb).  The consensus fills"
          " are seeded onto the K=60 unipath edges and bridged"
          " through the graph, validated against the reads.  On"
          " E. coli this helped 15 of 171 fills but could not"
          " improve the 9 emitted ones (rrn-adjacent repeats) --"
          " the MultipleAligner ONT consensus is the better"
          " default for ONT gaps");
     CommandArgument_Int_OrDefault_Doc(MIN_SEG, 3500,
          "mm2: drop per-edge path segments shorter than this before"
          " voting.  Swept on E. coli (500..10000): 3500 gives ONE"
          " circular 4.62 Mb scaffold, zero misjoins vs the"
          " reference (N50 34.7x PE-only); 1000 gives 10 paths"
          " (N50 5.4x); 0 gives fine-grained threading through"
          " variant and hub edges (26 paths).  Waypoint density"
          " scales with ONT coverage and genome structure -- re-tune"
          " for other datasets");
     CommandArgument_Bool_OrDefault_Doc(VERBOSE, False,
          "print scaffold links as they are found");
     EndCommandArguments;

     ForceAssert( !HBV.empty( ) && !FASTB.empty( ) && !READS.empty( ) );

     HyperBasevector hb;
     BinaryReader::readFile( HBV, &hb );
     vecbasevector edges(FASTB);
     vec<int> inv;
     BinaryReader::readFile( INV, &inv );

     // Load per-edge coverage from a.covs.
     vec<double> edge_cov( hb.E( ), 0.0 );
     if ( !COVS.empty( ) && IsRegularFile( COVS ) )
     {    vec<vec<covcount>> covs;
          BinaryReader::readFile( COVS, &covs );
          if ( covs.size( ) > 0 )
               for ( int e = 0; e < hb.E( ); e++ )
                    if ( covs[0][e].Def( ) )
                         edge_cov[e] = covs[0][e].Cov( );    }

     ScaffoldWithLongReads( hb, inv, edges, edge_cov, READS,
          K, MIN_READ, MIN_SUPPORT, ANCHOR_STEP,
          REPEAT_FACTOR, MIN_LINE, MAX_BFS_DEPTH,
          MAX_MULT, MIN_ANCHOR_BLOCK, MIN_CHAIN_COV,
          MAPPER, MINIMAP2, MIN_MAP_LEN, MIN_MAP_ID, MIN_SEG,
          DEBUG_READ, K60_FASTB, OUT_HEAD, VERBOSE );

     Scram(0);
}
