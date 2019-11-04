///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// AffineAlign.  Produce visual affine alignment of two sequences.  Sadly, the two
// sequences have to match up exactly on the ends, or something bad happens.

#include "FetchReads.h"
#include "MainTools.h"
#include "PrintAlignment.h"
#include "pairwise_aligners/SmithWatAffine.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(F1, "first fasta file; uses only first record");
     CommandArgument_String_Doc(F2, "second fasta file; uses only first record");
     CommandArgument_Bool_OrDefault_Doc(RC1, False, 
          "reverse complement first entry");
     CommandArgument_Bool_OrDefault_Doc(RC2, False, 
          "reverse complement second entry");
     CommandArgument_Int_OrDefault(MISMATCH, 3);
     CommandArgument_Int_OrDefault(GAP_OPEN, 12);
     CommandArgument_Int_OrDefault(GAP_EXTEND, 1);
     CommandArgument_Int_OrDefault_Doc(ID1, 0, "use this record on F1");
     CommandArgument_Int_OrDefault_Doc(ID2, 0, "use this record on F2");
     CommandArgument_Int_OrDefault_Doc(START2, 0, 
          "start at this zero-based position on ID2");
     CommandArgument_Int_OrDefault_Doc(STOP2, -1, 
          "stop at this zero-based position on ID2");
     EndCommandArguments;

     vecbasevector f1, f2;
     FetchReads( f1, 0, F1 );
     FetchReads( f2, 0, F2 );

     if( RC1){
         for( auto& entry: f1){
             entry.ReverseComplement();
         }
     }
     if( RC2){
         for( auto& entry: f2){
             entry.ReverseComplement();
         }
     }

     if ( STOP2 < 0 ) STOP2 = f2[ID2].size( );
     basevector X2( f2[ID2], START2, STOP2 - START2 );

     alignment a;
     SmithWatAffineParallel( f1[ID1], X2, a, true, true,
          MISMATCH, GAP_OPEN, GAP_EXTEND );

     PrintVisualAlignment( True, cout, f1[ID1], X2, a );    }
