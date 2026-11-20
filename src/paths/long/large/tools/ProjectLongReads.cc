///////////////////////////////////////////////////////////////////////////////

// ProjectLongReads: P2/LR-1 tool (docs/Modernization.md).  Project long
// reads (ONT/PacBio) onto the K=200 assembly graph by seeding exact
// k-mer matches on the (error-corrected) edge sequences and chaining
// collinear anchors.  See paths/long/large/LongReadProjection.{h,cc}.

#include "MainTools.h"
#include "paths/long/large/LongReadProjection.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(HBV, "",
          "a.hbv of the graph to project onto (e.g. OUT/a.200/a.hbv)");
     CommandArgument_String_OrDefault_Doc(FASTB, "",
          "a.fastb of the same graph");
     CommandArgument_String_OrDefault_Doc(READS, "",
          "long reads: fasta or fastq, possibly gzipped");
     CommandArgument_String_OrDefault_Doc(OUT_HEAD, "ont",
          "prefix for .proj.txt (text) and .proj (binary) outputs");
     CommandArgument_UnsignedInt_OrDefault_Doc(K, 19,
          "seed length for exact k-mer anchoring");
     CommandArgument_Int_OrDefault_Doc(MIN_READ, 20000,
          "minimum read length for acceptance statistics");
     CommandArgument_Int_OrDefault_Doc(ANCHOR_STEP, 64,
          "approximate spacing between anchors along a read");
     CommandArgument_Int_OrDefault_Doc(MAX_ANCHORS, 1024,
          "cap on anchors per read (even subsampling above this)");
     CommandArgument_Bool_OrDefault_Doc(VERBOSE, False,
          "print example chains");
     EndCommandArguments;

     ForceAssert( !HBV.empty( ) && !FASTB.empty( ) && !READS.empty( ) );

     HyperBasevector hb;
     BinaryReader::readFile( HBV, &hb );
     vecbasevector edges(FASTB);
     ForceAssertEq( int(edges.size( )), hb.EdgeObjectCount( ) );

     ProjectLongReads( hb, edges, READS, K, MIN_READ, ANCHOR_STEP,
          MAX_ANCHORS, OUT_HEAD, VERBOSE );

     Scram(0);
}
