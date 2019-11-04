///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef MAKE_LOOKUP_H
#define MAKE_LOOKUP_H

#include "Basevector.h"
#include "CoreTools.h"
#include "kmers/KmerRecord.h"

template<int K> void MakeKmerLookup0Pre( const vecbasevector& unibases,
     const String& prefix, vec< triple<kmer<K>,int,int> >& kmers_plus );

#endif
