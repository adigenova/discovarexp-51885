///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Test for hardcoded Broad-specific samples.

#include "CoreTools.h"
#include "FastIfstream.h"
#include "TokenizeString.h"
#include "paths/long/large/Samples.h"

String PicardBamList( const String& fc, const int n );
String PicardBams( const String& fc, const String& lanes, const String& lib );

void Samples( String& species, const String& SAMPLE, String& X, String& EVALUATE,
     String& SELECT_FRAC, int& READS_TO_USE, const String& DATASET, String& BAM,
     vec<String>& subsam_names )
{
     // Define datasets.

     vec<String> NA12878_datasets 
          = { "1", "X1", "X2", "X3", "X4", "newchem", "big", "bigger" };
     if ( SAMPLE == "NA12878" && !Member( NA12878_datasets, DATASET ) )
     {    cout << "Illegal DATASET." << endl;
          Scram(1);    }
     vec<String> NA12877_datasets = { "1", "2" };
     if ( SAMPLE == "NA12877" && !Member( NA12877_datasets, DATASET ) )
     {    cout << "Illegal DATASET." << endl;
          Scram(1);    }
     vec<String> rhody_datasets = { "1", "2" };
     if ( SAMPLE == "rhody" && !Member( rhody_datasets, DATASET ) )
     {    cout << "Illegal DATASET." << endl;
          Scram(1);    }

     // Define samples.

     vec<String> normal_human_cell_lines = { "NA12877", "NA12878", "HG00096", 
          "HG00268", "HG00419", "HG00759", "HG01051", "HG01112", "HG01500", 
          "HG01565", "HG01583", "HG01595", "HG01879", "HG02568", "HG02922", 
          "HG03006", "HG03052", "HG03642", "HG03742", "NA18525", "NA18939", 
          "NA19017", "NA19625", "NA19648", "NA20502", "NA20845",
          "YRI.1", "YRI.2", "YRI.3" };
     vec<String> normal_human_non_cell_lines
          = { "16E_MD", "15E_DD", "17E_PD", "F2.1", "F2.2", "F2.3", "25H_JM", 
               "F3.1", "F3.2", "F3.3", "F3.4", "F3.5" };
     vec<String> normal_human_trios = { "F1", "F2", "F3", "CEPH", "YRI" };
     vec<String> cancer_human_samples = { "HCC1143", "HCC1143BL", "HCC1143+BL", 
          "HCC1954", "HCC1954BL", "HCC1954+BL" };
     vec<String> nonhuman_samples = { "rhino", "aardvark", "maize", "tb148", 
          "tbHaarlem", "ecoli11", "ecoli12", "ecoli_scs", "rhody", "plasmo",
          "lolium", "arabidopsis" };
     for ( int s = 1; s <= 18; s++ )
          nonhuman_samples.push_back( "helico" + ToString(s) );
     for ( int s = 1; s <= 10; s++ )
          nonhuman_samples.push_back( "arabidopsis" + ToString(s) );

     // Handle subsamples.

     if ( SAMPLE == "HCC1143+BL" || SAMPLE == "HCC1954+BL" ) 
          subsam_names = { "T", "N" };
     else if ( SAMPLE == "F1" || SAMPLE == "F2" || SAMPLE == "CEPH" )
     {    subsam_names = { "C", "M", "P" };    }
     else if ( SAMPLE == "YRI" )
     {    subsam_names = { "P", "M", "C" };    }
     else if ( SAMPLE == "F3" )
     {    subsam_names = { "C1", "C2", "M", "P1", "P2" };    }
     else subsam_names = { "C" };

     // Continue.

     vec<String> human_samples = normal_human_cell_lines;
     human_samples.append(normal_human_non_cell_lines);
     human_samples.append(normal_human_trios);
     human_samples.append(cancer_human_samples);

     species = ( Member( human_samples, SAMPLE ) ? "human" : SAMPLE );
     vec<String> samples = human_samples;
     samples.append(nonhuman_samples);
     if ( !Member( samples, SAMPLE ) && SAMPLE != "" )
     {    cout << "Unknown sample." << endl;
          Scram(1);    }
     if ( species == "human" && SAMPLE != "NA12878" )
     {    if ( X == "" ) X = "all";
          EVALUATE = "False";    }
     Bool cov_unset = ( SELECT_FRAC == "" && READS_TO_USE < 0 );

     if ( SAMPLE == "HCC1143" && cov_unset ) SELECT_FRAC = "0.75";
     if ( SAMPLE == "HCC1143+BL" && cov_unset ) SELECT_FRAC = "0.85+0.65";
     if ( SAMPLE == "HCC1954BL" && cov_unset ) SELECT_FRAC = "0.5";
     if ( SAMPLE == "HCC1954+BL" && cov_unset ) SELECT_FRAC = "0.85+0.65";

     vec<String> sx = { "HG00096:H7AGFADXX", "HG00268:H7A1KADXX",
          "HG00419:H7A26ADXX", "HG00759:H7A29ADXX", "HG01051:H7A22ADXX",
          "HG01112:H7A2BADXX", "HG01500:H7A83ADXX", "HG01565:H7A63ADXX",
          "HG01583:H7A2NADXX", "HG01595:H7A7BADXX", "HG01879:H7A2JADXX",
          "HG02568:H7A2DADXX", "HG02922:H7A2VADXX", "HG03006:H7CLPADXX",
          "HG03052:H7B6EADXX", "HG03642:H7A2EADXX", "HG03742:H7A81ADXX",
          "NA18525:H7A33ADXX", "NA18939:H781JADXX", "NA19017:H7A7YADXX",
          "NA19625:H7A1NADXX", "NA19648:H7AGDADXX", "NA20502:H7A1WADXX",
          "NA20845:H7A30ADXX" };
     for ( int i = 0; i < sx.isize( ); i++ )
     {    String sample = sx[i].Before( ":" ), fc = sx[i].After( ":" );
          if ( SAMPLE == sample ) BAM = PicardBamList( fc, 2 );    }

     if ( SAMPLE.Contains( "YRI.", 0 ) ) // should switch to fastb/qualp
     {    int id = SAMPLE.After( "YRI." ).Int( );
          String sid;
          if ( id == 1 ) sid = "NA19238"; // father
          if ( id == 2 ) sid = "NA19239"; // mother
          if ( id == 3 ) sid = "NA19240"; // child
          BAM = "/humgen/1kg/processing/cnv/phase3_highcov/bams_reheadered/"
               + sid + ".sorted.clean.recal.reheader.bam";    }

     // Trios.

     if ( SAMPLE == "F1" && cov_unset ) SELECT_FRAC = "0.4+0.6+0.55";
     if ( SAMPLE == "F2" && cov_unset ) SELECT_FRAC = "0.68+0.68+0.68";
     if ( SAMPLE == "F3" && cov_unset ) SELECT_FRAC = "0.285+0.38+0.57+0.57+0.475";
     if ( SAMPLE == "CEPH" && cov_unset ) SELECT_FRAC = "0.65+0.55+0.65";
     if ( SAMPLE == "YRI" && cov_unset ) SELECT_FRAC = "0.7+0.7+0.7";

     if ( SAMPLE == "NA12877" ) // datasets that Illumina sent us
     {    String root = "/wga/scr4/vendor/illumina/2014-05-28/Conversion";
          if ( DATASET == "1" )
          {    BAM = root + "/140517_700518R_0287_Bh9fm3adxx/1/L1_NA12878/"
                    + "H9FM3ADXX.1.unmapped.bam" + ","
                    + root + "/140517_700518R_0287_Bh9fm3adxx/2/L2_NA12878/"
                    + "H9FM3ADXX.2.unmapped.bam";    }
          if ( DATASET == "2" )
          {    BAM = root + "/140517_D00277R_0126_Bh9flhadxx/1/L1_NA12878/"
                    + "H9FLHADXX.1.aligned.bam" + ","
                    + root + "/140517_D00277R_0126_Bh9flhadxx/2/L2_NA12878/"
                    + "H9FLHADXX.2.aligned.bam";    }    }

     // NA12878.big: This is from some data that Illumina sent us, based on newer
     // chemistry.  This dataset uses two lanes from H9TEU (from a standard 
     // library), and one lane from H9T96 (a wide library).

     if ( SAMPLE == "NA12878" && DATASET == "big" )
     {    String iroot = "/wga/scr4/vendor/illumina/2014-07-04/Conversion";
          String abam = ".aligned.bam";
          BAM = iroot + 
               "/140613_HSQ1185_0757_BH9TEUADXX/1/L1_NA12878/H9TEUADXX.1" + abam
               + "," + iroot + 
               "/140613_HSQ1185_0757_BH9TEUADXX/2/L2_NA12878/H9TEUADXX.2" + abam
               + "," + iroot + 
               "/140620_HSQ1185_0761_BH9T96ADXX/1/L1_NA12878/H9T96ADXX.1" 
               + abam;    }

     if ( SAMPLE == "NA12878" && DATASET == "newchem" )
     {    String iroot = "/wga/scr4/vendor/illumina/2014-07-04/Conversion";
          String abam = ".aligned.bam";
          BAM = iroot + 
               "/140613_HSQ1185_0757_BH9TEUADXX/1/L1_NA12878/H9TEUADXX.1" + abam
               + "," + iroot + 
               "/140613_HSQ1185_0757_BH9TEUADXX/2/L2_NA12878/H9TEUADXX.2" 
               + abam;    }

     if ( SAMPLE == "NA12878" && DATASET == "bigger" )
     {    String iroot = "/wga/scr4/vendor/illumina/2014-07-04/Conversion";
          String abam = ".aligned.bam";
          BAM = "/wga/scr4/picard/H01UJADXX/C1-508_2012-11-01_2012-11-04/1/"
               "Solexa-125532/H01UJADXX.1.aligned.duplicates_marked.bam," 
               + iroot + 
               "/140613_HSQ1185_0757_BH9TEUADXX/1/L1_NA12878/H9TEUADXX.1" + abam
               + "," + iroot + 
               "/140613_HSQ1185_0757_BH9TEUADXX/2/L2_NA12878/H9TEUADXX.2" + abam
               + "," + iroot + 
               "/140620_HSQ1185_0761_BH9T96ADXX/1/L1_NA12878/H9T96ADXX.1" 
               + abam;    }
     
     // F2
     // F2.1 = 23H_LM
     // F2.2 = 24H_CM
     // F2.3 = 25H_JM

     // F3
     // F3.1 = 65T_CR
     // F3.2 = 66T_NG
     // F3.3 = 67T_SR
     // F3.4 = 68T_DR
     // F3.5 = 69T_GG

     if ( SAMPLE == "NA12878" && DATASET == "X1" )
     {    BAM = "/seq/tng/datasets/20140317-hiseqx-wave2/genomes/formulation1_fca/"
               "NA12878/NA12878.bam";    }
     if ( SAMPLE == "NA12878" && DATASET == "X2" )
     {    BAM = "/seq/tng/datasets/20140317-hiseqx-wave2/genomes/formulation1_fcb/"
               "NA12878/NA12878.bam";    }
     if ( SAMPLE == "NA12878" && DATASET == "X3" )
     {    BAM = "/seq/tng/datasets/20140317-hiseqx-wave2/genomes/formulation2_fca/"
               "NA12878/NA12878.bam";    }
     if ( SAMPLE == "NA12878" && DATASET == "X4" )
     {    BAM = "/seq/tng/datasets/20140317-hiseqx-wave2/genomes/formulation2_fcb/"
               "NA12878/NA12878.bam";    }

     // Arabidopsis.

     if ( SAMPLE.Contains( "arabidopsis", 0 ) )
     {    int n = SAMPLE.After( "arabidopsis" ).Int( );
          BAM = "/seq/picard/H75EWADXX/*/*/Solexa-" 
               + ToString( 183228 - 1 + n ) + "/*.bam";    }

     // Butterflies.

     String hroot = "/seq/picard_aggregation/G77375";
     if ( SAMPLE.Contains( "helico", 0 ) )
     {    int n = SAMPLE.After( "helico" ).Int( );
          BAM = hroot;
          if ( n == 1 ) BAM += "/DAS_09-93/current/DAS_09-93";
          if ( n == 2 ) BAM += "/DAS_09-132/current/DAS_09-132";
          if ( n == 3 ) BAM += "/DAS_09-229/current/DAS_09-229";
          if ( n == 4 ) BAM += "/DAS_09-282/current/DAS_09-282";
          if ( n == 5 ) BAM += "/DAS_09-323/current/DAS_09-323";
          if ( n == 6 ) BAM += "/DAS_110-111/current/DAS_110-111";
          if ( n == 7 ) BAM += "/DAS_GELE/current/DAS_GELE";
          if ( n == 8 ) BAM += "/DAS_REL13_139/current/DAS_REL13_139";
          if ( n == 9 ) BAM += "/DAV_1_timareta/current/DAV_1_timareta";
          if ( n == 10 ) BAM += "/DAV_2_cydno/current/DAV_2_cydno";
          if ( n == 11 ) BAM += "/DAV_3_hecale/current/DAV_3_hecale";
          if ( n == 12 ) BAM += "/DAV_4_numata/current/DAV_4_numata";
          if ( n == 13 ) BAM += "/DAV_5_melpomene/current/DAV_5_melpomene";
          if ( n == 14 ) BAM += "/MCM_doris_6897/current/MCM_doris_6897";
          if ( n == 15 ) 
               BAM += "/MCM_eratoxhimera_F1_6700/current/MCM_eratoxhimera_F1_6700";
          if ( n == 16 ) 
               BAM += "/MCM_himera_father_6363/current/MCM_himera_father_6363";
          if ( n == 17 ) 
               BAM += "/MCM_himera_mother_6394/current/MCM_himera_mother_6394";
          if ( n == 18 ) BAM += "/REE_vanillae/current/REE_vanillae";
          BAM += ".bam";    }

     // Bacteria.

     if ( SAMPLE == "ecoli11" )
     {    BAM = "/seq/picard/A42CA/C1-516_2013-05-28_2013-05-30/1/"
               "Ec_11_9941_PCR_Free/A42CA.1.aligned.duplicates_marked.bam";    }
     if ( SAMPLE == "ecoli12" ) // reads = 2893480
     {    BAM = "/seq/picard/A6FWA/C1-516_2013-12-20_2013-12-22/1/Solexa-208789/"
               "A6FWA.1.aligned.duplicates_marked.bam";
          EVALUATE = "False"; // runs too long, takes too much memory
               }
     if ( SAMPLE == "ecoli_scs" )
     {    BAM = PicardBams( "A6FWA", "1", "Solexa-208786" );
          // Could add A6VP5 and could add Solexa-208787.
          // BAM = PicardBams( "A6VP5", "1", "Solexa-208786" );
          SELECT_FRAC = "0.7";    }
     if ( SAMPLE == "plasmo" )
     {    BAM = "frac:0.25::" + PicardBams( "A2RD9", "1", "Solexa-160369" );    }
     if ( SAMPLE == "rhody" ) // reads = 3981302
     {    BAM = "frac:0.5::" + PicardBams( "A6FWA", "1", "Solexa-208790" );
          EVALUATE = "False"; // evaluation crashed on 512 GB box
          if ( DATASET == "2" )
               BAM += "," + PicardBams( "A6VP5", "1", "Solexa-208790" );    }
     if ( SAMPLE == "tb148" )
     {    BAM = "/seq/picard/A6FWA/C1-516_2013-12-20_2013-12-22/1/"
               "Solexa-208791/A6FWA.1.aligned.duplicates_marked.bam";
          if (cov_unset) SELECT_FRAC = "0.4";    }
     if ( SAMPLE == "tbHaarlem" )
     {    BAM = "/seq/picard/A42CA/C1-516_2013-05-28_2013-05-30/1/Haarlem_PCR_Free/"
               "A42CA.1.aligned.duplicates_marked.bam";
          if (cov_unset) READS_TO_USE = 5000000;
          EVALUATE = "False";    }

     if ( ( SAMPLE == "rhino" || SAMPLE == "aardvark" ) && X != "all" )
     {    cout << "Rhino and aardvark only do all." << endl;
          Scram(1);    }
     if ( ( SAMPLE == "rhino" || SAMPLE == "aardvark" ) && EVALUATE == "True" )
     {    cout << "Rhino and aardvark can't be evaluated." << endl;
          Scram(1);    }
     if ( BAM != "" && SAMPLE != "NA12878" && X == "" ) X = "all";
     if ( SELECT_FRAC == "" ) SELECT_FRAC = "1.0";    }

String PicardBamList( const String& fc, const int n )
{    fast_pipe_ifstream in( "find /seq/picard/" + fc + " -name \"*.bam\" -print" );
     vec<String> lines;
     String line;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          lines.push_back(line);    }
     ForceAssertEq( lines.isize( ), n );
     String x;
     for ( int j = 0; j < n; j++ )
     {    if ( j > 0 ) x += ",";
          x += lines[j];    }
     return x;    }

String PicardBams( const String& fc, const String& lanes, const String& lib )
{    String answer;
     String path = "/seq/picard/" + fc;
     if ( !IsDirectory(path) )
     {    cout << "Can't find directory " << path << "." << endl;
          cout << "Abort." << endl;
          Scram(1);    }
     vec<String> dir1 = AllFiles(path);
     if ( !dir1.solo( ) )
     {    cout << "Directory " << path << " does not contain a unique entry." 
               << endl;
          cout << "Abort." << endl;
          Scram(1);    }
     path += "/" + dir1[0];
     vec<String> lanesv;
     Tokenize( lanes, ',', lanesv );
     for ( int l = 0; l < lanesv.isize( ); l++ )
     {    String path2 = path + "/" + lanesv[l] + "/" + lib;
          if ( !IsDirectory(path2) )
          {    cout << "Can't find directory " << path2 << "." << endl;
               cout << "Abort." << endl;
               Scram(1);    }
          String bam = fc + "." + lanesv[l] + ".aligned.duplicates_marked.bam";
          if ( !IsRegularFile( path2 + "/" + bam ) )
          {    cout << "Can't find file " << path2 + "/" + bam << "." << endl;
               cout << "Abort." << endl;
               Scram(1);    }
          if ( answer != "" ) answer += ",";
          answer += path2 + "/" + bam;    }
     return answer;    }
