///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
    "DISCOVAR de novo (experimental) is a de novo genome assembler that requires "
    "only a single PCR-free paired end Illumina library containing 250 base reads.";

// MakeDepend: library JEMALLOC

#include "MainTools.h"
#include "paths/long/large/GapToyCore.h"

int main(int argc, char *argv[])
{
     RunTime();

     BeginCommandArgumentsNoHeader;

     CommandDoc(DOC);

     CommandArgument_String_Doc(READS,
          "Comma-separated list of input files, see manual for details");
     CommandArgument_String_Doc(OUT_DIR, "name of output directory");
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS,0,
          "Number of threads.  By default, the number of processors online.");
     CommandArgument_String_OrDefault_Doc(REFHEAD, "",
          "use reference sequence REFHEAD.fasta to annotate assembly, and also "
          "REFHEAD.names if it exists");
     CommandArgument_Double_OrDefault_Doc( MAX_MEM_GB, 0,
          "if specified, maximum allowed RAM use in GB; in some cases may be "
          "exceeded by our code");
     EndCommandArguments;

     GapToyCore(argc,argv);

     Scram(0);    }
