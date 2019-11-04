///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// GapToy is a prototype assembler for large genomes, set up as an experimental
// test program.  It has the following steps:
// 1. Generate K=60 unipath graph from high-quality kmers (ReadQGrapher).
// 2. Map reads to paths on this graph (ReadQGrapher).
// 3. Convert paths to basevectors and re-unipath at K=200 (ReadQGrapher).
// 4. Translate the paths to K=200 (ReadQGrapher).
// 5. Identify gaps and carry out local assemblies as in LongProto/Discovar.
// 6. Use the local assemblies to patch the global assembly graph.
// 7. Evaluate using the Fosmid reference sequences.
//
// Results for X=all EVALUATE=True
//
// * means results inferred, not actually run verbatim
// (local): comparable numbers from local assisted assemblies, see Baseline*.cc
// run times are for crd26-29
//
//                            tig   scaff
// start   rev    time  edges N50    N50  mis gaps indels subs  notes
// date           hours  /kb  kb     kb
//
// 1/13/14 48330* 23.7  23.2                  227  179    697
// 1/17/14 48381  22.2   8.9                  205  201    818
// 1/29/14 48546  26.6   8.9                  169  183    746
// 2/21/14 48747  26.9   8.9  10.7            158  193    778
// 2/22/14 48760  37.7   8.9  10.8            145  194    778   PLACE_MORE
// 3/9/14  48944  26.6   8.4  18.7            158  193    752
// ------- 48989* ----   5.0  21.9            148  188    708   MIN_RATIO2=8
// 3/14/14 49051  28.3   4.5  34.0            120  176    872
// 4/8/14  49475  19.9   4.6  41.6             80  157    541
// 4/13/14 49530  21.0   4.4  50.2             76  167    542   
// 4/15/14 49553  19.2   4.4  76.7   91.2      76  167    546
// 4/15/14 49558  23.1   4.2  87.2  108.0      79  168    573
// 4/26/14 49703  23.4   3.2  92.5  122.1      81  161    557
// 5/9/14  49817  27.9   3.4  99.9  135.5  27  61  113    334
// 5/16/14 49855  29.2   3.4 100.3  153.0  27  70  110    401
// 6/10/14 49962  28.0   3.3 100.3  153.0  29  70  110    401
// 6/17/14 50010  26.5   3.3 102.3  156.1  27  70  110    381
// 6/22/14 50043  26.3   3.3 102.2  156.0  27  70  110    381
// 7/18/14 50189  26.8   3.3 102.2  156.1  26  70  110    381
// 7/22/14 50218  28.7   3.3 102.2  156.1  26  70  110    381
// 7/25/14 50261  29.1   3.3 102.2  156.1  26  70  110    381  41069 N50 perfect
// 8/9/14  50433  28.7   3.3 102.2  156.0  26  70  110    381  41069 N50 perfect
// 8/10/14 50434  27.9   3.3 102.2  156.3  28  69  114    384  41069 N50 perfect
// 8/19/14 50564  28.9   3.3 102.2  163.5  30  71  114    402  41069 N50 perfect
// ------- 50577* ----   3.3 102.8  165.3  32  71  114    413  41069 N50 perfect
// 8/21/14 50614  28.8   3.3 102.3  163.7  32  71  114    401  41069
// 8/22/14 50621  27.2   3.3 102.3  163.8  32  71  117    436  41069
// 8/25/14 50642  29.2   3.3 102.5  164.2  33  71  117    435  41069
// 9/10/14 50858  29.8   3.3 102.5  164.3  33  71  117    434  41069
// 10/16/14 51155 29.7   3.3 102.5  164.3  33  71  117    434  41069
// 10/28/14 51366 28.4   3.3 102.5  164.3  33  65  114    425  41069
// 10/29/14 *51383 28.4  3.3 103.3  164.4  33  65  115    467  41069
// 10/31/14 51400 29.4   3.3 103.3  164.4  33  65  115    467  41069
// 11/2/14  51418 30.0   3.3 103.3  164.4  33  65  115    467  41069
// 11/10/14 51454 30.0   3.3 103.3  164.4  33  65  115    467  41069
// 11/14/14 51464 28.5   3.3 103.3  164.4  33  65  115    467  41069
// 11/17/14 51472 29.6   3.3 103.3  164.4  33  65  115    467  41069
// 11/20/14 51481 28.4   3.3 103.3  164.4  33  65  115    467  41069
// 12/16/14 51648 28.8   3.3 103.3  164.4  33  65  115    467  41036
//
// (local)               1.6                    4   69    131
//
// RHINO
// 8/7/14  50408  16.4   2.6 115.9  171.7
//
// AARDVARK
// 8/8/14  50408  38.5   6.6  41.0  46.7
//
// CANCER1
// ------- 49305* ----   ???  35.0                             49264.HCC1954
// ------- 49333* ----   ???  34.9                             49264.HCC1954/a.fin2
// 4/9/14  49509* 33.6   4.3  43.3                             49498.HCC1954
//
// NORMAL1
// 3/29/14 49308  28.4   4.6  40.8                             49308.HCC1954BL
// ------- 49333* ----   4.6  40.7                             49308.HCC1954BL/a.fin2
//
// CANCER2
// 4/6/14  49441  33.0   2.7  57.3                             49441.HCC1143
// 4/10/14 49513  23.9   2.7  67.1                             49513.HCC1143
// ------- 49580* ----   2.6 114.6  165.4             49573.HCC1143/a.fin2/a.s
// 6/11/14 49962  42.9   2.3 102.3  235.2                      49962.HCC1143
//
// NORMAL2
// 4/1/14  49355  35.0   4.8  44.3                    49355.HCC1143BL/a.fin/a.u
// ------- 48592* ----   --- 121.2  169.8             49563.HCC1143BL/a.fin2/a.s
//
// CANCER2+NORMAL2
// 5/22/14 49875  54.2   3.1 104.0  210.2   7         49875.HCC1143+BL

// MakeDepend: library JEMALLOC
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

// MakeDepend: dependency MakeLookupTable
// MakeDepend: dependency QueryLookupTable

#include <omp.h>

#include "FastIfstream.h"
#include "FetchReads.h"
#include "Intvector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "VecUtilities.h"
#include "feudal/ObjectManager.h"
#include "feudal/PQVec.h"
#include "lookup/LookAlign.h"
#include "graph/DigraphTemplate.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/UnibaseUtils.h"
#include "paths/long/BuildReadQGraph.h"
#include "paths/long/PlaceReads0.h"
#include "paths/long/ReadPath.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/large/AssembleGaps.h"
#include "paths/long/large/Clean200.h"
#include "paths/long/large/DiscoStats.h"
#include "paths/long/large/ExtractReads.h"
#include "paths/long/large/FinalFiles.h"
#include "paths/long/large/GapToyCore.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/Lines.h"
#include "paths/long/large/MakeGaps.h"
#include "paths/long/large/Preclose.h"
#include "paths/long/large/Repath.h"
#include "paths/long/large/Samples.h"
#include "paths/long/large/Simplify.h"
#include "paths/long/fosmid/Fosmids.h"
#include "system/HostName.h"

int GapToyCore( int argc, char *argv[] )
{
     double all_clock = WallClockTime( );
     String start_time = Date( );

     BeginCommandArguments;

     // PUBLIC INTERFACE (should be same as in DiscoverExp.cc, 
     //                   except some of the args are optional here)

     CommandArgument_String_OrDefault_Doc(READS, "",
          "Comma-separated list of input files, see manual for details" );
     CommandArgument_String_OrDefault_Doc(OUT_DIR, "", "name of output directory");
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS,0,
          "Number of threads.  By default, the number of processors online.");
     CommandArgument_String_OrDefault_Doc(REFHEAD, "", 
          "use reference sequence REFHEAD.fasta to annotate assembly, and also "
          "REFHEAD.names if it exists");
     CommandArgument_Double_OrDefault_Doc( MAX_MEM_GB, 0,
          "if specified, maximum allowed RAM use in GB; in some cases may be "
          "exceeded by our code");

     // OTHER KEY OPTIONS

     CommandArgument_String_OrDefault_Doc(SAMPLE, "", "see Samples.cc");
     CommandArgument_String_OrDefault_Doc(X, "", "tiny, small, medium, "
          "mediump, large, all or allf, or e.g. 3:20M-20.1M to get the given "
          "region on chr3; only all for rhino and aardvark");
     CommandArgument_String_OrDefault_Doc(DATASET, "1",
          "1, X1, X2, X3 or X4, and only meaningful for SAMPLE=NA12878");
     CommandArgument_String_OrDefault_Doc(F, "default",
          "list of Fosmid ids, otherwise determined by X");
     CommandArgument_String_OrDefault_Doc(INSTANCE, "1",
          "to allow multiple concurrent runs");
     CommandArgument_String_OrDefault_Doc(FIN, "", 
          "suffix a.fin and a.final directory names with this");
     CommandArgument_String_OrDefault_Doc(SELECT_FRAC, "", 
          "fraction of data to use; only applies in certain cases; "
          "comma-separated pair of doubles for HCC1143+BL");
     CommandArgument_Int_OrDefault_Doc(READS_TO_USE, -1, 
          "bound for number of reads to use; only applies in certain cases; "
          "alternative to SELECT_FRAC");

     // CONTROL OVER WHAT PART OF PROGRAM IS RUN

     CommandArgument_String_OrDefault_Doc(EVALUATE, "",
          "evaluate assembly versus reference; default is False");
     CommandArgument_Bool_OrDefault_Doc(EVALUATE_ONLY, False,
          "just evaluate assembly versus reference");
     CommandArgument_Bool_OrDefault_Doc(SCAFFOLD_ONLY, False,
          "just generate scaffolds");
     CommandArgument_Bool_OrDefault_Doc(ALIGN_ONLY, False,
          "just align to reference");
     CommandArgument_Bool_OrDefault_Doc(CACHE, False, "use prebuilt data");
     CommandArgument_String_OrDefault_Doc(EXIT, "", 
          "exit at given point; only LOAD, PATHS and PATCHES allowed now");
     CommandArgument_Bool_OrDefault_Doc(START_WITH_PATCHES, False, 
          "start with assembly plus patches");
     CommandArgument_Bool_OrDefault_Doc(START_PATCHED, False, 
          "start with patched assembly");
     CommandArgument_Bool_OrDefault_Doc(START_FIN, False, 
          "start with a.fin assembly");
     CommandArgument_Bool_OrDefault_Doc(START_FINAL, False, 
          "start with a.final assembly");
     CommandArgument_Bool_OrDefault_Doc(PERFSTATS, False,
          "record performance statistics");

     // NON-ALGORITHMIC OPTIONS

     CommandArgument_Bool_OrDefault_Doc(ANNOUNCE, False, 
          "announce start and stop of main parallel loop instances");
     CommandArgument_String_OrDefault_Doc(FETCH_PIDS, "",
          "start,stop -- assuming that there is just one Fosmid, find the "
          "pair ids incident upon this interval on the Fosmid reference sequence, "
          "and exit; also find the original edges");
     CommandArgument_Int_OrDefault_Doc(FETCH_PIDS_MQ, 30,
          "max qual score for aligning reads for FETCH_PIDS");
     CommandArgument_String_OrDefault_Doc(DUMP_LOCAL, "",
          "lroot,rroot -- dump local dot file associated to this");
     CommandArgument_String_OrDefault_Doc(USER, "",
          "use this instead of Getenv(USER)");
     CommandArgument_String_OrDefault_Doc(ROOT, "/wga/scr4", 
          "root for output directory location");
     CommandArgument_Int_OrDefault_Doc(PAD, 30000, "flank for X=allf");
     CommandArgument_Bool_OrDefault_Doc(EVALUATE_VERBOSE, False, 
          "verbose evaluation output");
     CommandArgument_Bool_OrDefault_Doc(PLACE_MORE_VERBOSE, False, 
          "verbose output for PLACE_MORE");
     CommandArgument_Bool_OrDefault_Doc(PLACE_MORE_DUMP, False, 
          "write aux files for PLACE_MORE");
     CommandArgument_Bool_OrDefault_Doc(ANALYZE_BRANCHES_VERBOSE, False, 
          "verbose output for analyze branches");
     CommandArgument_Bool_OrDefault_Doc(ANALYZE_BRANCHES_VERBOSE2, False, 
          "verbose output for analyze branches 2");
     CommandArgument_Bool_OrDefault_Doc(RQGRAPHER_VERBOSE, False,
          "verbose output from ReadQGrapher");
     CommandArgument_Bool_OrDefault_Doc(SELECT_SPECIAL, False, 
          "find reads that we ought to separately correct");
     CommandArgument_Bool_OrDefault_Doc(PAIRS_60, False, 
          "generate pairs file for K=60");
     CommandArgument_Bool_OrDefault_Doc(CREATE_LOCAL_DIR, False,
          "create the directory 'local' in the working directory");
     CommandArgument_Bool_OrDefault_Doc(DOT_MATCHING, False, 
          "call PrintDotMatchingGenome");
     CommandArgument_Bool_OrDefault_Doc(ALIGN_TO_GENOME, True, 
          "call AlignToGenome");
     CommandArgument_Bool_OrDefault_Doc(PULL_APART_VERBOSE, False,
          "when using PullApart print copious diagnostic information");
     CommandArgument_IntSet_OrDefault_Doc(PULL_APART_TRACE,"{}",
             "trace these edges during PullAparter");
     CommandArgument_Bool_OrDefault_Doc(SCAFFOLD_VERBOSE, False,
          "generate verbose output when scaffolding");
     CommandArgument_Int_OrDefault_Doc(CLEAN_200_VERBOSITY, 0,
          "to generate verbose output for CLEAN_200");
     CommandArgument_String_OrDefault_Doc(TRACE_SEQ, "",
          "trace the given DNA sequence; very partially implemented; "
          "puts files in fin_dir/trace");
     CommandArgument_Int_OrDefault_Doc(GAP_CAP, -1,
          "to cap number of gap assemblies, for timing experiments only");
     CommandArgument_Bool_OrDefault_Doc(KEEP_NAMES, False,
          "keep read names (partially implemented)");
     CommandArgument_Bool_OrDefault_Doc(SAVE_60, False, "create a.60 directory");

     // ALGORITHMIC OPTIONS

     CommandArgument_Bool_OrDefault_Doc(INJECT, False, 
          "inject nearby global edges into local assembly");
     CommandArgument_Bool_OrDefault_Doc(FILL_JOIN, False, 
          "turn on FILL_GAPS and JOIN_OVERLAPS in ReadQGrapher");
     CommandArgument_Bool_OrDefault_Doc(EXTRA_PATHS, False, 
          "try to path more reads");
     CommandArgument_Bool_OrDefault_Doc(CUTTING_EDGE, False, 
          "turn on a bunch of cutting-edge options");
     CommandArgument_Bool_OrDefault_Doc(EXTEND, False, "extend terminal edges");
     CommandArgument_Int_OrDefault_Doc(MIN_FREQ, 3, "passed to ReadQGrapher");
     CommandArgument_Int_OrDefault_Doc(MIN_QUAL, 7, "passed to ReadQGrapher");
     CommandArgument_Bool_OrDefault_Doc(LOCAL_LAYOUT, False, 
          "lay out reads in local assembly, not functional yet");
     CommandArgument_Bool_OrDefault_Doc(CONSERVATIVE_KEEP, False, "use more "
          "conservative definition of keep; does not seem to help globally");
     CommandArgument_Int_OrDefault_Doc(K2_FLOOR, 0, 
          "minimum value for K2 in local assemblies; setting to 60 does not seem "
          "to help with global assemblies");
     CommandArgument_Bool_OrDefault_Doc(PLACE_MORE, False, 
          "try to place more reads");
     CommandArgument_Int_OrDefault_Doc(PLACE_MORE_LEVEL, 1, 
          "1 or 2; governs cpu/memory tradeoff; 2 uses less memory and does not "
          "increase cpu a lot");
     CommandArgument_Bool_OrDefault_Doc(ANALYZE_BRANCHES, False, 
          "analyze branches, seems to increase number of edges in global assembly");
     CommandArgument_Bool_OrDefault_Doc(ANALYZE_BRANCHES_REV, False, 
          "analyze branches in reverse direction");
     CommandArgument_Bool_OrDefault_Doc(KEEP_ALL_LOCAL, False,
          "if true, then all local assemblies are given unique TMP "
               "directories such that reads and corrected reads are saved");
     CommandArgument_UnsignedInt_OrDefault_Doc(JOIN_PATHS, 0,
          "rewrite K=60 ReadPaths to merge overlapping pairs: 0=none (default), "
               "1=zero-order method, 2=new method, 3=new method+trace;"
               "2 possibly somewhat helpful");
     CommandArgument_Bool_OrDefault_Doc(EXTEND_PATHS, False,
          "for purposes of generating the K=200 graph, extend K=60 paths through "
          "unambiguous regions of the K=60 graph; does not appear to help global "
          "assembly");
     CommandArgument_IntSet_OrDefault_Doc(TRACE_PATHS, "{}",
               "paths to trace in AddNewStuff (and elsewhere?)");
     CommandArgument_Int_OrDefault_Doc(MAX_DEL2, 200, "passed to RemoveHangs");
     CommandArgument_Int_OrDefault_Doc(MIN_RATIO2, 8, "passed to AnalyzeBranches");
     CommandArgument_Bool_OrDefault_Doc(PRECLOSE, False, "preclose at K=60 stage");
     CommandArgument_Bool_OrDefault_Doc(PF_ONLY, False, 
          "use only PF reads; note that this resizes the non-PF reads to zero "
          "rather than eliminate them altogether");
     CommandArgument_Int_OrDefault_Doc(MIN_GAIN, 5, "for read extension");
     CommandArgument_Int_OrDefault_Doc(MAX_SUPP_DEL, 0, 
          "for deletion of low_support edges; increasing to 1 might help");
     CommandArgument_Bool_OrDefault_Doc(EXTEND_PAIRS_60, False, 
          "extend pairs; does not seem to help");
     CommandArgument_Bool_OrDefault_Doc(TAMP_EARLY, True, 
          "add early call to tamp down hanging ends; of ambiguous value as it "
          "trades off recovery of a small number of SNPs for a small increase "
          "in total edge count");
     CommandArgument_Int_OrDefault_Doc(EXT_MODE, 1, "path extension mode, 1 or 2; "
          "mode 2 doesn't seem to help but worth trying again");
     CommandArgument_Bool_OrDefault_Doc(PLACE_PARTNERS, False, 
          "try to place partners on a.patched");
     CommandArgument_Bool_OrDefault_Doc(SHORT_KMER_READ_PATHER, False,
          "Tries to find read paths for unplaced reads");
     CommandArgument_Bool_OrDefault_Doc(CYCLIC_SAVE, True,
          "try to recover some local assemblies that contain cycles");
     CommandArgument_Bool_OrDefault_Doc(CLEAN_200_POST, False,
          "clean some weak edges from patched K=200 assembly; "
          "seems to do some damage but retry later");
     CommandArgument_Int_OrDefault_Doc(MAX_CELL_PATHS, 50, "for finding lines");
     CommandArgument_Int_OrDefault_Doc(MAX_DEPTH, 10, "for finding lines");
     CommandArgument_Int_OrDefault_Doc(A2V, 5, "version for A2");
     CommandArgument_Int_OrDefault_Doc(CLEAN_200V, 2, 
          "version for CLEAN_200; note option of using 3");
     CommandArgument_Int_OrDefault_Doc(MIN_LINE, 5000, "for scaffolding");
     CommandArgument_Int_OrDefault_Doc(MIN_LINK_COUNT, 3, "for scaffolding; "
          "lowering MIN_LINK_COUNT to 2 increases scaffold length but also "
          "results in more misassemblies");
     CommandArgument_Int_OrDefault_Doc(MAX_PROX_LEFT, 400, "for pair selection");
     CommandArgument_Int_OrDefault_Doc(MAX_PROX_RIGHT, 400, "for pair selection");
     CommandArgument_Bool_OrDefault_Doc(REMOVE_TINY, False,
          "remove tiny standalone edges in a.200; seems slightly deleterious");
     CommandArgument_Bool_OrDefault_Doc(EXT_FINAL, True, "extend paths near end");
     CommandArgument_Int_OrDefault_Doc(EXT_FINAL_MODE, 1, "0 or 1");
     CommandArgument_Bool_OrDefault_Doc(DEGLOOP, True,
          "remove weakly supported branches");
     CommandArgument_Int_OrDefault_Doc(DEGLOOP_MODE, 1, "1 or 2");
     CommandArgument_Double_OrDefault(DEGLOOP_MIN_DIST, 2.5);
     CommandArgument_Bool_OrDefault_Doc(IMPROVE_PATHS, True, "run ImprovePath");
     CommandArgument_Bool_OrDefault(IMPROVE_PATHS_LARGE, False);
     CommandArgument_Bool_OrDefault_Doc(FINAL_TINY, True,
          "remove tiny stuff at end");
     CommandArgument_Bool_OrDefault_Doc(GAP_CLEANUP, True,
          "clean up after making scaffolds");
     CommandArgument_Bool_OrDefault_Doc(UNWIND3, True, "unwind three-edge plasmids");
     CommandArgument_Int_OrDefault_Doc(MAX_BPATHS, 100000,
          "cap for local assembly expansion");

     EndCommandArguments;

     // Set computational limits.

     SetThreads(NUM_THREADS);
     int64_t max_bytes = 0;
     if ( MAX_MEM_GB > 0 )
     {    if ( MAX_MEM_GB < 1 )
          {    cout << "\nPlease don't set MAX_MEM_GB to less than 1.  It is "
                    << "likely to make our code crash.\n" << endl;
               Scram(1);    }
          int64_t max_bytes 
               = int64_t( round( MAX_MEM_GB * 1024.0 * 1024.0 * 1024.0 ) );
          SetMaxMemory(max_bytes);    }

     // Check arguments.

     ForceAssertLe(JOIN_PATHS,3u);

     // Parse DUMP_LOCAL.

     int DUMP_LOCAL_LROOT=-1, DUMP_LOCAL_RROOT=-1;      // -1 can't match edge
     if ( DUMP_LOCAL != "" ) 
     {    DUMP_LOCAL_LROOT = DUMP_LOCAL.Before(",").Int();
          DUMP_LOCAL_RROOT = DUMP_LOCAL.After(",").Int();    }

     // Implement CUTTING_EDGE.

     if (CUTTING_EDGE)
     {    INJECT = True;            // doesn't seem to help on X=all
          FILL_JOIN = True;         // doesn't seem to help on X=all
          EXTRA_PATHS = True;       // dies on X=all
          PLACE_MORE = True;        // helps X=all but slow
          ANALYZE_BRANCHES = True;  // doesn't seem to help on X=all
               }

     // Set up directories, etc.

     String user = ( USER != "" ? USER : Getenv("USER","noname") );
     String work_dir = ROOT + "/" + user + "/GapToy/" + INSTANCE;
     if ( OUT_DIR != "" ) work_dir = OUT_DIR;
     Mkpath(work_dir);
     if ( PERFSTATS ) PerfStatLogger::init(work_dir);

     // Load reference.
     
     if ( REFHEAD != "" )
     {    if ( !IsRegularFile( REFHEAD + ".fasta" ) )
          {    cout << "Can't find REFHEAD.fasta." << endl;
               Scram(1);    }
          vecbasevector genome;
          FetchReads( genome, 0, REFHEAD + ".fasta" );
          genome.WriteAll( work_dir + "/genome.fastb" );
          if ( IsRegularFile( REFHEAD + ".names" ) )
          {    vec<String> rnames;
               fast_ifstream in( REFHEAD + ".names" );
               String line;
               while(1)
               {    getline( in, line );
                    if ( in.fail( ) ) break;
                    rnames.push_back(line);    }
               if ( genome.size( ) != rnames.size( ) )
               {    cout << "The file REFHEAD.fasta has " << genome.size( )
                         << " records but\nthe file REFHEAD.names has "
                         << rnames.size( ) << " lines.  These should be the same."
                         << endl;
                    Scram(1);    }
               Cp( REFHEAD + ".names", work_dir + "/genome.names" );    }    }

     // Handle samples and subsamples.  A given sample may be 
     // divided into multiple subsamples.

     String species;
     vec<String> subsam_names;
     Samples( species, SAMPLE, X, EVALUATE, SELECT_FRAC, 
          READS_TO_USE, DATASET, READS, subsam_names );
     BinaryWriter::writeFile( work_dir + "/subsam.names", subsam_names );
     vec<int64_t> subsam_starts( subsam_names.size( ), 0 );

     // Define heuristics.

     const int K0 = 60;
     const int K = 200;

     String fin_dir = work_dir + "/a.fin" + FIN;
     String final_dir = work_dir + "/a.final" + FIN;
     Mkdir777(fin_dir);
     Mkdir777(final_dir);
     Echo( ToString(K), fin_dir + "/a.k" );
     Echo( ToString(K), final_dir + "/a.k" );
     Remove( work_dir + "/clock.log" );
     LogTime( 0, "", work_dir );
     {    string hostname = getHostName();
          hostname = hostname.substr( 0, hostname.find('.') );
          OfstreamMode( xout, work_dir + "/the_command", ios::app );
          xout << "\n" << hostname << ": " << command.TheCommand( ) << endl;    }
     String run_head = work_dir + "/a";
     String tmp_dir1 = work_dir + "/data";
     Mkdir777(tmp_dir1);
     if( CREATE_LOCAL_DIR) Mkdir777(work_dir+"/local");
     String log_dir = work_dir + "/logs";
     Mkdir777(log_dir);
     SystemSucceed( "/bin/rm -rf " + work_dir + "/loc/*" );
     Mkdir777( work_dir + "/loc" );
     Mkdir777( work_dir + "/special" );
     if ( !CACHE && !START_WITH_PATCHES && !START_PATCHED && !START_FIN 
          && !START_FINAL && !EVALUATE_ONLY && !ALIGN_ONLY && !SCAFFOLD_ONLY
          && REFHEAD == "" )
     {    Remove( work_dir + "/genome.fastb" );    }

     // Output README.

     if ( true )
     {
         Ofstream( rout, work_dir + "/README" );
         rout << "Assembly directories:\n";
         // rout << "a.60      - initial graph from 60-mers" << endl;
         rout << "a.200     - initial graph from 200-mers" << endl;
         rout << "a.patched - patched graph" << endl;
         rout << "a.fin     - pre-scaffolding final graph" << endl;
         rout << "a.final   - final assembly" << endl;
     }

     // Scaffold only if so requested.

     if (SCAFFOLD_ONLY)
     {    HyperBasevector hbx;
          BinaryReader::readFile( fin_dir + "/a.hbv", &hbx );
          vec<int> inv2;
          BinaryReader::readFile( fin_dir + "/a.inv", &inv2 );
          ReadPathVec paths2( fin_dir + "/a.paths" );
          VecULongVec invPaths;
          invert( paths2, invPaths, hbx.EdgeObjectCount( ) );
          MakeGaps( hbx, inv2, paths2, invPaths, MIN_LINE, MIN_LINK_COUNT,
               work_dir, "fin" + FIN, SCAFFOLD_VERBOSE, GAP_CLEANUP );
          BinaryWriter::writeFile( final_dir + "/a.hbv", hbx );
          BinaryWriter::writeFile( final_dir + "/a.inv", inv2 );
          TestInvolution( hbx, inv2 );
          Scram(0);    }

     // Define fosmids, regions, and results.

     vec<String> regions;
     Bool all = False;
     vec<int> fosmids;
     map<String,GapToyResults> res;
     if ( X == "tiny" )
     {    fosmids = {1};
          regions.push_back( "1:24.8M-24.9M" );      // F1
          }
     else if ( X.size() > 3 && X.StartsWith("fos") ) // X=fosY for fosmid Y
     {    int fosno = X.After("fos").Int();
          fosmids = {fosno};
          regions.push_back( FosmidRegion(fosno, PAD) );
          cout << "assembling fosmid id = " << fosno << endl;
          }
     else if ( X == "small" ) // 3.45 Mb
     {    for ( int n = 1; n <= 39; n++ )
          {    if ( n == 16 || n == 33 || n == 35 || n == 36 || n == 38 ) continue;
               fosmids.push_back(n);    }
          res[X].rev = 49066;
          res[X].nedges = 13963, res[X].meanlen = 566.536;
          res[X].gaps = 26, res[X].indels = 43, res[X].subs = 109;
          // Notes for particular Fosmids, run with CUTTING_EDGE = True.  
          // * = good candidate for investigation
          regions.push_back( "1:24.8M-24.9M" );      // F1
          // F2. There is a small gap around 43100.  There is a gap in the same 
          // place in the unassisted LongProto assembly.
          regions.push_back( "1:54.75M-54.85M" );    // F2
          regions.push_back( "1:164.9M-165M" );      // F3
          regions.push_back( "2:106.75M-106.85M" );  // F4
          // F5.  Two gaps.  The first gap is not a real gap but the alternative
          // has a 102-base indel.  The second gap is 294 bases of low-complexity
          // sequence.  It is closed in the unassisted LongProto assembly.
          regions.push_back( "2:239.3M-239.4M" );    // F5
          regions.push_back( "3:11M-11.1M" );        // F6
          regions.push_back( "3:61.5M-61.6M" );      // F7
          regions.push_back( "5:111M-111.1M" );      // F8
          regions.push_back( "5:177.65M-177.75M" );  // F9
          // F9.  Two gaps.  One gap can be eliminated with CONSERVATIVE_KEEP=True
          // and K2_FLOOR=60.
          regions.push_back( "5:179.3M-179.4M" );    // F10
          // F11.  One gap.  The same gap is in the unassisted LongProto assembly.
          regions.push_back( "5:179.2M-179.3M" );    // F11
          // F12.  One gap, new.
          regions.push_back( "6:19.6M-19.7M" );      // F12
          regions.push_back( "6:92.55M-92.65M" );    // F13
          regions.push_back( "7:3.85M-3.95M" );      // F14
          regions.push_back( "7:38.7M-38.8M" );      // F15
          // F17.  Two gaps.
          regions.push_back( "8:23.2M-23.3M" );      // F17
          regions.push_back( "8:30.75M-30.85M" );    // F18
          // F19.  One gap, new.
          regions.push_back( "8:72.75M-72.85M" );    // F19
          // *F20.  One gap.  New to r48420.  This involves a huge homopolymer-like
          // cell that contains the true path but which DeleteLowCoverage wrecks.
          // Turning off DeleteLowCoverage also results in a gap because because
          // the cell contains a cycle.
          regions.push_back( "8:128.75M-128.85M" );  // F20
          regions.push_back( "10:30.85M-30.95M" );   // F21
          // *F22.  One gap.  New to r48420.  This gap can be closed by *deleting*
          // two pids: 6295, 6708.  This suggests that we may need a better local
          // assembly process.
          regions.push_back( "11:44.9M-45M" );       // F22
          regions.push_back( "11:45.5M-45.6M" );     // F23
          // F24.  Two large gaps.  The first one is around 10300 and the second one
          // is around 33300.  If you assemble the region using unassisted 
          // LongProto, the first gap is filled but associated with a 655 base 
          // deletion, and the second gap is a gap in the LongProto assembly.
          regions.push_back( "11:64.95M-65.05M" );   // F24
          regions.push_back( "11:67.7M-67.8M" );     // F25
          regions.push_back( "11:75.45M-75.55M" );   // F26
          regions.push_back( "11:111.8M-111.9M" );   // F27
          regions.push_back( "12:3.1M-3.2M" );       // F28
          // F29.  Two gaps.  One is a large insertion so we don't expect to get it.
          // The other gap is unclosed in unassisted LongProto, but looks like it
          // might be closable.
          regions.push_back( "12:7M-7.1M" );         // F29
          regions.push_back( "12:14.8M-14.9M" );     // F30
          // *F31.  One gap.  Unassisted LongProto assembly does not have a gap.
          // Single Fosmid assembly has the same gap.
          regions.push_back( "12:57.55M-57.7M" );    // F31
          regions.push_back( "12:113.95M-114.05M" ); // F32
          regions.push_back( "14:104M-104.1M" );     // F34
          // F37.  Three gaps.  Not examined yet.
          regions.push_back( "15:30.4M-30.5M" );     // F37
          regions.push_back( "15:74.9M-75.0M" );     // F39
               }
     else if ( X == "medium" ) // 6 Mb
     {    fosmids = {1,2,3,4,5,6,7,8,9,10,11,12};
          res[X].rev = 48385;
          res[X].nedges = 31518, res[X].meanlen = 450.469;
          res[X].gaps = 9, res[X].indels = 14, res[X].subs = 24;
          regions.push_back( "1:24.6M-25.1M" );      // F1
          regions.push_back( "1:54.55M-55.05M" );    // F2
          regions.push_back( "1:164.7M-165.2M" );    // F3
          regions.push_back( "2:106.55M-107.05M" );  // F4
          regions.push_back( "2:239.1M-239.6M" );    // F5
          regions.push_back( "3:10.8M-11.3M" );      // F6
          regions.push_back( "3:61.3M-61.8M" );      // F7
          regions.push_back( "5:110.8M-111.3M" );    // F8
          regions.push_back( "5:177.45M-177.95M" );  // F9
          regions.push_back( "5:178.8M-179.8M" );    // F10,11
          regions.push_back( "6:19.4M-19.9M" );      // F12
               }
     else if ( X == "mediump" ) // 11.5 Mb
     {    fosmids = {1,2,3,4,5,6,7,8,9,10,11,12};
          res[X].rev = 48385;
          res[X].nedges = 62015, res[X].meanlen = 439.684;
          res[X].gaps = 11, res[X].indels = 17, res[X].subs = 36;
          regions.push_back( "1:24.4M-25.4M" );      // F1
          regions.push_back( "1:54.35M-55.35M" );    // F2
          regions.push_back( "1:164.5M-165.5M" );    // F3
          regions.push_back( "2:106.35M-107.35M" );  // F4
          regions.push_back( "2:238.9M-239.9M" );    // F5
          regions.push_back( "3:10.6M-11.6M" );      // F6
          regions.push_back( "3:61.1M-62.1M" );      // F7
          regions.push_back( "5:110.6M-111.6M" );    // F8
          regions.push_back( "5:177.25M-178.25M" );  // F9
          regions.push_back( "5:178.6M-180.1M" );    // F10,11
          regions.push_back( "6:19.2M-20.2M" );      // F12
               }
     else if ( X == "chr11" )                        // chromosome 11
     {
          fosmids = {22,23,24,25,26,27,88,89,90,91};
          regions.push_back("11");
     }
     else if ( X == "large" ) // the entire aligned genome
     {    fosmids = {1,2,3,4,5,6,7,8,9,10,11,12};
          for ( int c = 1; c <= 22; c++ )
               regions.push_back( ToString(c) );
          regions.push_back( "X" );    }
     else if ( X == "allf" ) 
     {    fosmids = AllFosmids( );
	  regions = AllFosmidRegions(PAD);    } // all fosmids with padding
     else if ( X == "all" ) // the entire genome
     {    if ( SAMPLE == "NA12878" ) 
          {    fosmids = AllFosmids( );
               if ( EVALUATE == "" ) EVALUATE = "True";    }
          all = True;    }
     else if ( X.Contains( ":" ) && X.After( ":" ).Contains( "-" ) )
     {    regions.push_back(X);    }
     else
     {    cout << "Unknown X argument." << endl;
          Scram(1);    }
     if ( F != "default" ) 
     {    vec<int> ftmp;
          ParseIntSet( F, ftmp );
          if ( !Subset(ftmp, fosmids) ) 
          {    cout << "WARNING: your eval fosmids may not intersect with the "
                    << "loaded data" << endl;    }
          fosmids = ftmp;
     }

     // Fetch pids.

     if ( FETCH_PIDS != "" ) FetchPids(FETCH_PIDS, fosmids, work_dir, FETCH_PIDS_MQ);

     // Load fosmid reference sequences.

     vecbasevector G; 
     vecString Gnames; 
     if ( SAMPLE == "NA12878" )
     {    G.reserve(fosmids.size());
          Gnames.reserve(fosmids.size());
          for ( int i = 0; i < fosmids.isize( ); i++ )
          {    vecbasevector g;
               int f = fosmids[i];
               FetchReads( g, 0, "/wga/dev/references/Homo_sapiens/"
                    "NA12878_Fosmid_Pool.regions.fin/fos." + ToString(f) 
                    + ".fasta" );
               if ( g.size() == 1 )
               { G.push_back(g[0]); Gnames.push_back( "F" + ToString(f) ); }
               else
               { String baseName = "F" + ToString(f) +'.';
                 for ( size_t idx = 0; idx != g.size(); ++idx )
                 { G.push_back(g[idx]);
                   Gnames.push_back(baseName+ToString(idx)); } } }    }
     if ( SAMPLE == "HCC1954" || SAMPLE == "HCC1954BL" )
     {    FetchReads( 
               G, 0, "/wga/dev/references/Homo_sapiens/HCC1954/genome.fasta" );
          for ( int g = 0; g < (int) G.size( ); g++ )
               Gnames.push_back( "B" + ToString(g+1) );    }

     // Run START_FINAL.

     if (START_FINAL)
     {    HyperBasevector hb;
          BinaryReader::readFile( final_dir + "/a.hbv", &hb );
          vec<int> inv2;
          BinaryReader::readFile( final_dir + "/a.inv", &inv2 );
          ReadPathVec paths2( final_dir + "/a.paths" );
          FinalFiles( hb, inv2, paths2, subsam_names, subsam_starts, work_dir, 
               final_dir, MAX_CELL_PATHS, MAX_DEPTH, ALIGN_TO_GENOME, EVALUATE, 
               EVALUATE_VERBOSE, X, res, SAMPLE, species, fosmids, G );
          cout << endl;
          Scram(0);    }

     // Evaluate if that's all that's asked for.

     if (EVALUATE_ONLY)
     {    HyperBasevector hb;
          BinaryReader::readFile( fin_dir + "/a.hbv", &hb );
          GapToyEvaluate( SAMPLE, species, hb, G, fosmids, work_dir, "fin" + FIN, 
               501, res[X], EVALUATE_VERBOSE );
          Scram(0);     }
     if (ALIGN_ONLY)
     {    HyperBasevector hb;
          BinaryReader::readFile( fin_dir + "/a.hbv", &hb );
          vec<int> inv2;
          BinaryReader::readFile( fin_dir + "/a.inv", &inv2 );
          vec< vec< pair<int,int> > > hits, hits_alt;
          vecbasevector genome( work_dir + "/genome.fastb" );
          AlignToGenome( hb, inv2, genome, hits );
          BinaryWriter::writeFile( fin_dir + "/a.aligns", hits );
          vecbasevector genome_alt;
          if ( IsRegularFile( work_dir + "/genome.fastb_alt" ) )
               genome_alt.ReadAll( work_dir + "/genome.fastb_alt" );
          HyperBasevector hbx;
          BinaryReader::readFile( final_dir + "/a.hbv", &hbx );
          BinaryReader::readFile( final_dir + "/a.inv", &inv2 );
          vec< vec< pair<int,int> > > hitsx, hitsx_alt;
          AlignToGenome( hbx, inv2, genome, hitsx );
          BinaryWriter::writeFile( final_dir + "/a.aligns", hitsx );
          if ( genome_alt.size( ) > 0 )
          {    AlignToGenome( hbx, inv2, genome_alt, hitsx_alt );
               BinaryWriter::writeFile( 
                    final_dir + "/a.aligns_alt", hitsx_alt );    }
          Scram(0);    }

     // Extract reads, convert reads to HyperBasevector and paths.

     vecbvec bases;
     ObjectManager<VecPQVec> quals(tmp_dir1 + "/frag_reads_orig.qualp");
     if ( !CACHE && !START_WITH_PATCHES && !START_PATCHED && !START_FIN )
     {    ExtractReads( SAMPLE, species, READS, SELECT_FRAC, READS_TO_USE, 
               regions, tmp_dir1, work_dir, all, PF_ONLY, KEEP_NAMES, 
               subsam_names, subsam_starts, &bases, quals );
          BinaryWriter::writeFile( work_dir + "/subsam.starts", subsam_starts );
          disco_stats stats;
          stats.Compute(work_dir);
          BinaryWriter::writeFile( work_dir + "/disco_stats", stats );
          if ( EXIT == "LOAD" ) Scram(0);

          // Report memory.

          {    uint64_t phys_mem = physicalMemory( );
               cout << Date( ) << ": see total physical memory of "
                    << ToStringAddCommas(phys_mem) << " bytes" << endl;
               rlimit max_data;
               if ( getrlimit( RLIMIT_DATA, &max_data ) == 0 )
               {    uint64_t max_mem = max_data.rlim_max;
                    if ( max_mem < phys_mem )
                    {    cout << "\nWell this is very interesting.  Apparently "
                              << "your memory usage is capped at\n" 
                              << ToStringAddCommas(max_mem)
                              << ".  This is less than the physical memory on\n" 
                              << "your machine, and may result in a crash.  "
                              << "Please let us know\nif this happens, as we can "
                              << "make our code respect the "
                              << "memory cap.\n\n";    }    }
               uint64_t allowed_mem = GetMaxMemory( );
               if ( allowed_mem > 0 && allowed_mem < phys_mem )
               {    cout << Date( ) << ": see user-imposed limit on memory of "
                         << ToStringAddCommas(allowed_mem) << " bytes" 
                         << endl;    }    }
     
          // Report bytes per base and issue warning if appropriate.

          int64_t total_bytes = GetMaxMemory( );
          double bytes_per_base = total_bytes / stats.total_bases;
          cout << Date( ) << ": " << setiosflags(ios::fixed) 
               << setprecision(2) << bytes_per_base << resetiosflags(ios::fixed)
               << " bytes per read base, " 
               << "assuming max memory available" << endl;
          // 2.25 OK in 51452.YRI
          // 1.85 failed in 51454.F3
          if ( bytes_per_base < 1.8 )
          {    cout << "\nWARNING: generally 2.0 bytes per read base of memory are "
                    << "needed.  You have\nsubstantially less than this, so the "
                    << "odds of your assembly completing are low.\n" << endl;    }
          if ( bytes_per_base < 2.0 )
          {    cout << "\nWARNING: generally 2.0 bytes per read base of memory are "
                    << "needed.  You have\nsomewhat less than this, so it is "
                    << "possible that your assembly will not finish.\n" << endl;    }

          if ( X != "all" )
          {    double mclock = WallClockTime( );
               SystemSucceed( "MakeLookupTable SOURCE=" + work_dir + "/genome.fastb "
                    "OUT_HEAD=" + work_dir + "/genome LO=True > " + log_dir 
                    + "/MakeLookupTable.out" );
               cout << TimeSince(mclock) << " used in MakeLookupTable" << endl;    }
          double rclock = WallClockTime( );
          if ( true )
          {    HyperBasevector hbv;
               ReadPathVec paths;
               buildReadQGraph( bases, quals,
                    FILL_JOIN, //doFillGaps
                    FILL_JOIN, //doJoinOverlaps
                    MIN_QUAL, MIN_FREQ,
                    .75, //minFreq2Fract
                    0, //maxGapSize
                    "", //refFasta
                    True, // new aligner
 		    SHORT_KMER_READ_PATHER, &hbv, &paths, RQGRAPHER_VERBOSE );
               cout << Date( ) << ": back from buildReadQGraph" << endl;
               cout << "memory in use = " 
                    << ToStringAddCommas( MemUsageBytes( ) ) << endl;
               vecbvec edges( hbv.Edges().begin(),hbv.Edges().end() );
               vec<int> inv;
               hbv.Involution(inv);

	       if (PRECLOSE) 
               {    Bool preclose_verbose = False;
                    Preclose( bases, quals.load(), work_dir, hbv, inv, paths,
                         preclose_verbose );
                    edges.resize( hbv.EdgeObjectCount( ) );
                    for ( int e = 0; e < hbv.EdgeObjectCount( ); e++ )
                         edges[e] = hbv.EdgeObject(e);    }

	       if (EXTEND_PAIRS_60)
               {    ExtendPairs60( hbv, inv, paths );
                    edges.resize( hbv.EdgeObjectCount( ) );
                    for ( int e = 0; e < hbv.EdgeObjectCount( ); e++ )
                         edges[e] = hbv.EdgeObject(e);    }
               String outHead = work_dir + "/a.60/a";
               if (SAVE_60)
               {    Mkdir777( work_dir + "/a.60" );
                    Echo( ToString(K0), work_dir + "/a.60/a.k" );
                    BinaryWriter::writeFile( outHead + ".hbv", hbv );
                    BinaryWriter::writeFile( outHead + ".hbx", 
                         HyperBasevectorX(hbv) );
                    edges.WriteAll( outHead+".fastb" );
                    BinaryWriter::writeFile( outHead + ".inv", inv );    }
               FixPaths( hbv, paths );
               if (SAVE_60)
               {    paths.WriteAll( outHead+".paths" );
                    VecULongVec invPaths;
                    invert( paths, invPaths, hbv.EdgeObjectCount( ) );
                    invPaths.WriteAll( outHead+".paths.inv");    }
               // if (PAIRS_60)
               // {    vec< pair<vec<int>,vec<int>> > pairs;
               //      vec<int64_t> pairs_pid;
               //      DefinePairs( 
               //           paths, inv, pairs, pairs_pid, work_dir + "/a.60" );    }

               if ( EXIT == "PATHS" )
               {    Scram(0);    }

               if ( JOIN_PATHS == 1 ) JoinPaths0( inv, paths );
               else if ( JOIN_PATHS == 2 ) JoinPaths( inv, paths, hbv );
               else if ( JOIN_PATHS == 3 ) JoinPaths( inv, paths, hbv, true );

               int64_t checksum_60 = hbv.CheckSum( );
               PRINT(checksum_60);

               quals.unload();
               Repath( hbv, edges, inv, paths, hbv.K(), 200, 
                    run_head+".200", True, True, EXTEND_PATHS );    }

          // Clean.

          HyperBasevector hb;
          vec<int> inv;
          double clock1 = WallClockTime( );
          BinaryReader::readFile( work_dir + "/a.200.hbv", &hb );
          if ( hb.E( ) == 0 )
          {    cout << "\nThe initial assembly is empty, so it is impossible "
                    << "to proceed.\nYou might try providing more data." << endl;
               Scram(0);    }
          BinaryReader::readFile( work_dir + "/a.200.inv", &inv );
          ReadPathVec paths( work_dir + "/a.200.paths" );
          cout << TimeSince(clock1) << " used reloading assembly" << endl;
          if ( CLEAN_200V <= 2 )
          {    Clean200( hb, inv, paths, bases, quals.load(), CLEAN_200_VERBOSITY,
                    CLEAN_200V, REMOVE_TINY );    }
          else
          {    Clean200x( hb, inv, paths, bases, quals.load(), CLEAN_200_VERBOSITY,
                    CLEAN_200V, REMOVE_TINY );    }

          // Write files.

          BinaryWriter::writeFile( work_dir + "/a.200.hbv", hb );
          BinaryWriter::writeFile( work_dir + "/a.200.hbx", HyperBasevectorX(hb) );
          BinaryWriter::writeFile( work_dir + "/a.200.inv", inv );
          paths.WriteAll( work_dir + "/a.200.paths" );
          vecbasevector edges( hb.EdgeObjectCount( ) );
          for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
               edges[e] = hb.EdgeObject(e);
          edges.WriteAll( work_dir + "/a.200.fastb" );
          if ( true ) // scope paths_index
          {
              VecULongVec paths_index;
              invert( paths, paths_index, hb.EdgeObjectCount( ) );
              paths_index.WriteAll( work_dir + "/a.200.paths.inv" );
          }

          // Relabel files.

          Mkdir777( work_dir + "/a.200" );
          Echo( ToString(K), work_dir + "/a.200/a.k" );
          for ( String s : { "fastb", "hbv", "hbx", "inv", "paths" } )
               Mv( work_dir + "/a.200." + s, work_dir + "/a.200/a." + s );
          Mv( work_dir+"/a.200.paths.inv", work_dir+"/a.200/a.paths.inv" );
          cout << TimeSince(rclock) << " used in ReadQGrapher" << endl;    }
     else
     {
         if ( !START_FIN )
         {
             double bqclock = WallClockTime( );
             bases.ReadAll(tmp_dir1 + "/frag_reads_orig.fastb");
             cout << TimeSince(bqclock) << " used loading bases" << endl;
         }
         BinaryReader::readFile( work_dir + "/subsam.starts", &subsam_starts );
     }

     // Proceed.

     HyperBasevector hb;
     vec<int> inv2;
     ReadPathVec paths2;
     if (START_PATCHED)
     {    double clock = WallClockTime( );
          BinaryReader::readFile( work_dir + "/a.patched/a.hbv", &hb );
          BinaryReader::readFile( work_dir + "/a.patched/a.inv", &inv2 );
          paths2.ReadAll( work_dir + "/a.patched/a.paths" );
          cout << TimeSince(clock) << " used loading assembly" << endl;    }
     if (START_FIN)
     {    double clock = WallClockTime( );
          BinaryReader::readFile( work_dir + "/a.fin/a.hbv", &hb );
          BinaryReader::readFile( work_dir + "/a.fin/a.inv", &inv2 );
          paths2.ReadAll( work_dir + "/a.fin/a.paths" );
          cout << TimeSince(clock) << " used loading assembly" << endl;    }
     if ( SAMPLE != "NA12878" && SAMPLE != "HCC1954" 
          && IsRegularFile( work_dir + "/genome.fastb" ) )
     {    G.ReadAll( work_dir + "/genome.fastb" );
          for ( int g = 0; g < (int) G.size( ); g++ )
               Gnames.push_back( ToString(g+1) );    }
     if ( !START_PATCHED && !START_FIN )
     {    double lclock = WallClockTime( );
          if (CACHE)
          {    cout << Date( ) << ": using " << ToStringAddCommas( bases.size( ) )
                    << " reads" << endl;   }
          cout << TimeSince(lclock) << " used reloading reads" << endl;

          // Load K=200 HBV and read paths.

          double xclock = WallClockTime( );
          BinaryReader::readFile( work_dir + "/a.200/a.hbv", &hb );
          int64_t checksum_200 = hb.CheckSum( );
          PRINT(checksum_200);
          int nedges = hb.EdgeObjectCount( );
          BinaryReader::readFile( work_dir + "/a.200/a.inv", &inv2 );
          paths2.ReadAll( work_dir + "/a.200/a.paths" );
          cout << "1 peak mem usage = " << PeakMemUsageGB( ) << " GB\n";
          cout << TimeSince(xclock) << " used loading stuff" << endl;

          // Dump initial assembly.

          if (DOT_MATCHING)
               PrintDotMatchingGenome( hb, G, Gnames, work_dir+"/a.200" );

          // Try to place more reads.

          if (EXTRA_PATHS) ExtraPaths( hb, bases, quals.load(), paths2 );
          if (PLACE_MORE) {
			vec<int64_t> placed;
			PlaceMore(hb, bases, quals.load(), paths2, placed,
                            PLACE_MORE_LEVEL, PLACE_MORE_VERBOSE);
			if (PLACE_MORE_DUMP) { 
			    paths2.WriteAll(work_dir + "/a.200/a.place_more.paths");
			    Ofstream(out, work_dir + "/a.200/a.place_more.ids");
			    for (size_t i = 0; i < placed.size(); i++)
				out << placed[i] << endl;
			}
	  }

          // Look for special reads that ought to be error corrected.

          if (SELECT_SPECIAL)
          {    SelectSpecials( hb, bases, quals.load(), paths2, work_dir );
               Scram(0);    }

          // Analyze branches.

        if ( ANALYZE_BRANCHES )
        {
            vec<int> to_right;
            hb.ToRight(to_right);
            AnalyzeBranches(hb, to_right, inv2, paths2, ANALYZE_BRANCHES_REV,
                            MIN_RATIO2, ANALYZE_BRANCHES_VERBOSE);
        }
          cout << "2 peak mem usage = " << PeakMemUsageGB( ) << " GB\n";

          // Assemble gaps.

          double nclock;
          if ( true ) // scope new_stuff
          {
              vecbvec new_stuff;
              if (START_WITH_PATCHES)
              {    BinaryReader::readFile( work_dir + "/new_stuff", &new_stuff );

                   // For experiments, should delete later.

                   if ( IsRegularFile( work_dir + "/new_stuff.exp" ) )
                   {    vecbvec new_stuff2;
                        new_stuff2.ReadAll( work_dir + "/new_stuff.exp" );
                        new_stuff.Append(new_stuff2);    }    }
              else
              {    VecULongVec paths2_index;
                   invert( paths2, paths2_index, hb.EdgeObjectCount( ) );
                   cout << "launching gap assemblies, mem usage = "
                        << ToStringAddCommas( MemUsageBytes( ) ) << endl;
                   AssembleGaps2( hb, inv2, paths2, paths2_index, bases, quals.load(),
                        work_dir, EXTEND, ANNOUNCE, KEEP_ALL_LOCAL, 
                        CONSERVATIVE_KEEP, INJECT, LOCAL_LAYOUT, DUMP_LOCAL, 
                        K2_FLOOR, DUMP_LOCAL_LROOT, DUMP_LOCAL_RROOT, new_stuff, 
                        CYCLIC_SAVE, A2V, GAP_CAP, MAX_PROX_LEFT, MAX_PROX_RIGHT,
                        MAX_BPATHS );
                   BinaryWriter::writeFile( 
                        work_dir + "/new_stuff", new_stuff );    }

              if ( EXIT == "PATCHES" )
              {    cout << "total time = " << TimeSince( all_clock, 1.0, "hours" )
                        << endl;
                   Scram(0);    }

              PRINT( new_stuff.size( ) );

              // Add in new stuff and rebuild graph.

              nclock = WallClockTime( );
              AddNewStuff( new_stuff, hb, inv2, paths2, bases, quals.load(),
                                  MIN_GAIN, TRACE_PATHS, work_dir, EXT_MODE );
          }

          // Stick partners on ends.

          PartnersToEnds( hb, paths2, bases, quals.load() );

          cout << TimeSince(nclock) << " used in new phase" << endl;
          PRINT2( hb.N( ), hb.EdgeObjectCount( ) ); // XXXXXXXXXXXXXXXXXXXXXXXX

          // Write modified assembly.

          TestInvolution( hb, inv2 );
          Mkdir777( work_dir + "/a.patched" );
          Echo( ToString(K), work_dir + "/a.patched/a.k" );
          double tclock = WallClockTime( );
          if (DOT_MATCHING)
               PrintDotMatchingGenome( hb, G, Gnames, work_dir + "/a.patched" );
          vecbvec(hb.Edges().begin(),hb.Edges().end()).WriteAll( 
               work_dir + "/a.patched/a.fastb" );
          BinaryWriter::writeFile( work_dir + "/a.patched/a.hbv", hb );
          BinaryWriter::writeFile( work_dir + "/a.patched/a.hbx", 
               HyperBasevectorX(hb) );
          BinaryWriter::writeFile( work_dir + "/a.patched/a.inv", inv2 );
          paths2.WriteAll( work_dir + "/a.patched/a.paths" );
          if ( true ) // scope invPaths
          {
              VecULongVec invPaths;
              invert( paths2, invPaths, hb.EdgeObjectCount( ) );
              invPaths.WriteAll( work_dir + "/a.patched/a.paths.inv" );
          }
          LogTime( tclock, "in tail 1" );    }
     Validate( hb, inv2, paths2 );

     // Simplify the assembly.

     if ( !START_FIN )
     {    if (CLEAN_200_POST)
               Clean200( hb, inv2, paths2, bases, quals.load(), CLEAN_200_VERBOSITY );
          Simplify( fin_dir, hb, inv2, paths2, bases, quals.load(), MAX_SUPP_DEL,
               TAMP_EARLY, MIN_RATIO2, MAX_DEL2, PLACE_PARTNERS, 
               ANALYZE_BRANCHES_VERBOSE2, TRACE_SEQ, DEGLOOP, EXT_FINAL,
               EXT_FINAL_MODE, PULL_APART_VERBOSE, PULL_APART_TRACE,
               DEGLOOP_MODE, DEGLOOP_MIN_DIST, IMPROVE_PATHS, 
               IMPROVE_PATHS_LARGE, FINAL_TINY, UNWIND3 );
          TestInvolution( hb, inv2 );
          if ( hb.E( ) == 0 )
          {    cout << "\nSimplified assembly is empty, so it is impossible to "
                    << "proceed.\nYou might try supplying more data.\n" << endl;
               Scram(0);    }
          TestInvolution( hb, inv2 );    }

     quals.unload();

     // Write pre-scaffolded final assembly.

     if ( !START_FIN )
     {    double pclock = WallClockTime( );
          bases.destroy( );
          BinaryWriter::writeFile( fin_dir + "/a.hbv", hb );
          BinaryWriter::writeFile( fin_dir + "/a.hbx", HyperBasevectorX(hb) );
          vecbasevector afin;
          for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
               afin.push_back( hb.EdgeObject(e) );
          afin.WriteAll( fin_dir + "/a.fastb" );
          if (DOT_MATCHING) PrintDotMatchingGenome( hb, G, Gnames, fin_dir );

          // For now, fix paths and write them and their inverse.

          for ( int i = 0; i < (int) paths2.size( ); i++ )
          {    Bool bad = False;
               for ( int j = 0; j < (int) paths2[i].size( ); j++ )
                    if ( paths2[i][j] < 0 ) bad = True;
               if (bad) paths2[i].resize(0);    }
          paths2.WriteAll( fin_dir + "/a.paths" );
          VecULongVec invPaths;
          invert( paths2, invPaths, hb.EdgeObjectCount( ) );
          invPaths.WriteAll( fin_dir + "/a.paths.inv" );
          LogTime( pclock, "writing final assembly" );
          BinaryWriter::writeFile( fin_dir + "/a.inv", inv2 );
          hb.DumpFasta( fin_dir + "/a.fasta", False );    }

     // Align to genome.

     if ( ALIGN_TO_GENOME && IsRegularFile( work_dir + "/genome.fastb" ) )
     {    vec< vec< pair<int,int> > > hits;
          vecbasevector genome( work_dir + "/genome.fastb" );
          AlignToGenome( hb, inv2, genome, hits );
          BinaryWriter::writeFile( fin_dir + "/a.aligns", hits );    }
     else Remove( fin_dir + "/a.aligns" );

     // Find lines and write files.

     vec<vec<vec<vec<int>>>> lines;
     FindLines( hb, inv2, lines, MAX_CELL_PATHS, MAX_DEPTH );
     BinaryWriter::writeFile( fin_dir + "/a.lines", lines );
     {  
	 vec<int> llens, npairs;
	 GetLineLengths( hb, lines, llens );
	 GetLineNpairs( hb, inv2, paths2, lines, npairs );
	 BinaryWriter::writeFile( fin_dir + "/a.lines.npairs", npairs );
	 vec<vec<covcount>> covs;
	 ComputeCoverage( hb, inv2, paths2, lines, subsam_starts, covs );
	 BinaryWriter::writeFile( fin_dir + "/a.covs", covs );
	 WriteLineStats( fin_dir + "/a", lines, llens, npairs, covs ); 
	 // Report CN stats
	 double cn_frac_good = CNIntegerFraction(hb, covs);
         cout << "CN fraction good = " << cn_frac_good << endl;
	 PerfStatLogger::log("cn_frac_good",ToString(cn_frac_good,2), "fraction of edges with CN near integer" );
     }
     // TestLineSymmetry( lines, inv2 );

     // Compute fragment distribution.

     FragDist( hb, inv2, paths2, work_dir + "/frags.dist" );

     // Scaffold.

     {    VecULongVec invPaths;
          invert( paths2, invPaths, hb.EdgeObjectCount( ) );
          MakeGaps( hb, inv2, paths2, invPaths, MIN_LINE, MIN_LINK_COUNT, work_dir,
               "fin" + FIN, SCAFFOLD_VERBOSE, GAP_CLEANUP );    }

     // Carry out final analyses and write final assembly files.

     FinalFiles( hb, inv2, paths2, subsam_names, subsam_starts, work_dir, final_dir,
          MAX_CELL_PATHS, MAX_DEPTH, ALIGN_TO_GENOME, EVALUATE, EVALUATE_VERBOSE, 
          X, res, SAMPLE, species, fosmids, G );

     // Done.

     cout << "\nrun started " << start_time << endl;
     cout << "peak mem usage = " << ToString( PeakMemUsageGB( ), 1 ) << " GB, ";
     cout << "total time = " << TimeSince( all_clock, 1.0, "hours" ) << endl;
     cout << "final checksum = " << hb.CheckSum( ) << endl;
     cout << "\n" << command.TheCommand( ) << endl << endl;

     double hours = double( WallClockTime( ) - all_clock ) / 3600.0;
     PerfStatLogger::log("etime_h", hours, "elapsed time in hours");

     // Note funny exit here, probably shouldn't be doing this.
     Scram(0);    }
