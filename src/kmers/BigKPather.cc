///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * BigKPather.cc
 *
 *  Created on: Dec 13, 2013
 *      Author: tsharpe
 */

#include "kmers/BigKPather.h"
#include "feudal/HashSet.h"
#include "kmers/BigKMer.h"
#include "kmers/KMerHasher.h"
#include "paths/long/HBVFromEdges.h"
#include "paths/long/LargeKDispatcher.h"
#include "system/SpinLockedData.h"
#include "system/WorklistN.h"
#include <iostream>

namespace
{

template <unsigned BIGK>
using BigDict = HashSet<BigKMer<BIGK>,typename BigKMer<BIGK>::hasher>;


template <unsigned BIGK>
class BigKMerizer
{
public:
    typedef BigDict<BIGK> BigKDict;
    BigKMerizer( BigKDict* pDict ) : mDict(*pDict) {}

    void kmerize( bvec const& bv, uint32_t seqIdx )
    { if ( bv.size() < BIGK ) return;
      KMerHasher<BIGK> hasher;
      auto itr = bv.begin();
      size_t hash = hasher.hash(itr);
      if ( bv.size() == BIGK ) { add(BigKMer<BIGK>(bv,hash,seqIdx)); return; }
      BigKMer<BIGK> kmer(bv,hash,seqIdx,
                         KMerContext::initialContext(bv[BIGK]));
      add(kmer);
      auto end = bv.end()-BIGK;
      while ( ++itr != end )
      { hash = hasher.stepF(itr);
        kmer.successor(hash,KMerContext(itr[-1],itr[BIGK]));
        add(kmer); }
      hash = hasher.stepF(itr);
      kmer.successor(hash,KMerContext::finalContext(itr[-1]));
      add(kmer); }

    // We build the dictionary with BigKMers pointing to reads.  This updates
    // all the BigKMers in the dictionary to point to edges instead.
    // seqIdx must be the edge's registry index (places first, then edges).
    void updateDict( bvec const& edge, uint32_t seqIdx )
    { KMerHasher<BIGK> hasher;
      auto itr = edge.begin();
      BigKMer<BIGK> kmer(edge,hasher.hash(itr),seqIdx);
      update(kmer.isRev()?kmer.rc():kmer);
      auto end = edge.end()-BIGK+1;
      while ( ++itr != end )
      { kmer.successor(hasher.stepF(itr));
        update(kmer.isRev()?kmer.rc():kmer); } }

private:
    void add( BigKMer<BIGK> const& kmer )
    { canonicalAdd( kmer.isRev() ? kmer.rc() : kmer ); }

    void update( BigKMer<BIGK> const& kmer )
    { const_cast<BigKMer<BIGK>*>(mDict.lookup(kmer))->updateLocation(kmer); }

    void canonicalAdd( BigKMer<BIGK> const& kmer ) const
    { mDict.apply(kmer,
        [&kmer]( BigKMer<BIGK> const& entry )
        { const_cast<BigKMer<BIGK>&>(entry).addContext(kmer.getContext()); }); }

    BigKDict& mDict;
};

template <unsigned BIGK>
class BigKEdgeBuilder
{
public:
    typedef BigKMer<BIGK> BigKmer;
    typedef BigKMerizer<BIGK> BigKmerizer;
    typedef BigDict<BIGK> BigKDict;

    static void buildEdges( BigKDict const& dict, vecbvec* pEdges,
                            std::vector<bvec const*> const* pReg )
    {
        BigKEdgeBuilder eb(dict,pEdges);
        dict.parallelForEachHHS(
                [eb,pReg]( typename BigKDict::HHS const& hhs ) mutable
                { BigKMer<BIGK>::setRegistry(pReg);
                  BigKMer<BIGK>::setScratch(&eb.mEdgeSeq);
                  for ( BigKmer const& entry : hhs )
                    if ( entry.isUnassigned() )
                      eb.buildEdge(entry); });

        size_t nRegularEdges = pEdges->size();
        size_t nTotalLength = pEdges->SizeSum();

        // if a kmer isn't marked as being on an edge, it must be a part of a
        // smooth circle: add those edges, too.  simpleCircle method isn't
        // thread-safe, so this part is single-threaded.
        BigKMer<BIGK>::setRegistry(pReg);
        BigKMer<BIGK>::setScratch(&eb.mEdgeSeq);
        for ( auto const& hhs : dict )
            for ( auto const& entry : hhs )
                if ( entry.isUnassigned() )
                    eb.simpleCircle(entry);
    }

private:
    BigKEdgeBuilder( BigKDict const& dict, vecbvec* pEdges )
    : mDict(dict), mEdges(*pEdges) {}

    void buildEdge( BigKmer const& entry )
    { if ( isPalindrome(entry) )
        make1KmerEdge(entry);
      else if ( upstreamExtensionPossible(entry) )
      { if ( downstreamExtensionPossible(entry) )
          return;
        BigKmer rc = entry.rc();
        mEdgeSeq.assign(rc.begin(),rc.end());
        mEdgeEntries.push_back(&entry);
        extend(rc.getContext()); }
      else if ( downstreamExtensionPossible(entry) )
      { mEdgeSeq.assign(entry.begin(),entry.end());
        mEdgeEntries.push_back(&entry);
        extend(entry.getContext()); }
      else
        make1KmerEdge(entry); }

    // not thread-safe
    void simpleCircle( BigKmer const& entry )
    { BigKmer const* pFirstEntry = &entry;
      mEdgeSeq.assign(entry.begin(),entry.end());
      mEdgeEntries.push_back(pFirstEntry);
      auto itr = mEdgeSeq.begin();
      KMerHasher<BIGK> hasher;
      BigKmer kmer(mEdgeSeq,hasher.hash(itr));
      KMerContext context = entry.getContext();
      while ( true )
      { ForceAssertEq(context.getPredecessorCount(),1u);
        ForceAssertEq(context.getSuccessorCount(),1u);
        mEdgeSeq.push_back(context.getSingleSuccessor());
        kmer.successor(hasher.stepF(++itr));
        BigKmer const* pEntry = lookup(kmer,&context);
        if ( pEntry == pFirstEntry )
        { mEdgeSeq.pop_back();
          break; }
        if ( !pEntry->isUnassigned() )
        { std::cout << "Failed to close circle.\n";
          for ( auto beg=mEdgeEntries.begin(),end=mEdgeEntries.end();
                  beg != end; ++beg )
            std::cout << *beg << ' ' << **beg;
          std::cout << pEntry << ' ' << *pEntry << std::endl;
          CRD::exit(1); }
        mEdgeEntries.push_back(pEntry); }
      canonicalizeCircle();
      addEdge(); }

    void canonicalizeCircle()
    { auto itr = std::min_element(mEdgeEntries.begin(),mEdgeEntries.end(),
                            []( BigKmer const* pEnt1, BigKmer const* pEnt2 )
                            { return *pEnt1 < *pEnt2; });
      BigKmer const& minKmer = **itr;
      size_t idx = itr - mEdgeEntries.begin();
      if ( !std::equal(minKmer.begin(),minKmer.end(),mEdgeSeq.begin(idx)) )
      { mEdgeSeq.ReverseComplement();
        idx = mEdgeSeq.size()-idx-BIGK;
        Assert(std::equal(minKmer.begin(),minKmer.end(),mEdgeSeq.begin(idx))); }
      bvec bv; bv.reserve(mEdgeSeq.size());
      bv.assign(mEdgeSeq.begin(idx),mEdgeSeq.end());
      bv.append(mEdgeSeq.begin(BIGK-1),mEdgeSeq.begin(BIGK-1+idx));
      Assert(std::equal(mEdgeSeq.begin(),mEdgeSeq.begin(BIGK-1),
                          mEdgeSeq.end()-(BIGK-1)));
      mEdgeSeq = bv; }

    bool isPalindrome( BigKmer const& kmer )
    { auto itr = kmer.begin();
      if ( !(BIGK&1) )
        return CF<BIGK>::isPalindrome(itr);
      if ( CF<BIGK-1>::isPalindrome(itr) )
          return true;
      return CF<BIGK-1>::isPalindrome(++itr); }

    bool upstreamExtensionPossible( BigKmer const& entry )
    { KMerContext context = entry.getContext();
      if ( context.getPredecessorCount() != 1 )
        return false;
      mEdgeSeq.clear().push_back(context.getSinglePredecessor());
      mEdgeSeq.append(entry.begin(),entry.end()-1);
      if ( CF<BIGK>::isPalindrome(mEdgeSeq.begin()) )
        return false;
      BigKmer kmer(mEdgeSeq,KMerHasher<BIGK>()(mEdgeSeq.begin()));
      lookup(kmer,&context);
      return context.getSuccessorCount() == 1; }

    bool downstreamExtensionPossible( BigKmer const& entry )
    { KMerContext context = entry.getContext();
      if ( context.getSuccessorCount() != 1 )
        return false;
      mEdgeSeq.assign(entry.begin()+1,entry.end());
      mEdgeSeq.push_back(context.getSingleSuccessor());
      if ( CF<BIGK>::isPalindrome(mEdgeSeq.begin()) )
        return false;
      BigKmer kmer(mEdgeSeq,KMerHasher<BIGK>()(mEdgeSeq.begin()));
      lookup(kmer,&context);
      return context.getPredecessorCount() == 1; }

    void make1KmerEdge( BigKmer const& entry )
    { mEdgeSeq.assign(entry.begin(),entry.end());
      mEdgeEntries.push_back(&entry);
      addEdge(); }

    void extend( KMerContext context )
    { auto itr = mEdgeSeq.begin();
      KMerHasher<BIGK> hasher;
      BigKmer kmer(mEdgeSeq,hasher.hash(itr));
      while ( context.getSuccessorCount() == 1 )
      { mEdgeSeq.push_back(context.getSingleSuccessor());
        kmer.successor(hasher.stepF(++itr));
        if ( isPalindrome(kmer) )
        { mEdgeSeq.pop_back();
          break; }
        BigKmer const* pEntry = lookup(kmer,&context);
        if ( context.getPredecessorCount() != 1 )
        { mEdgeSeq.pop_back();
          break; }
        mEdgeEntries.push_back(pEntry); }
      switch (mEdgeSeq.getCanonicalForm())
      {case CanonicalForm::PALINDROME:
         ForceAssertEq(mEdgeSeq.size(),BIGK);
        // allow flow-through to FWD case
       case CanonicalForm::FWD:
         addEdge(); break;
       case CanonicalForm::REV:
         mEdgeSeq.clear(); mEdgeEntries.clear(); break; } }

    BigKmer const* lookup( BigKmer const& kmer, KMerContext* pContext )
    { BigKmer const* result;
      if ( kmer.isRev() )
      { result = mDict.lookup(kmer.rc());
        ForceAssert(result);
        *pContext = result->getContext().rc(); }
      else
      { result = mDict.lookup(kmer);
        ForceAssert(result);
        *pContext = result->getContext(); }
      return result; }

    void addEdge()
    { static SpinLockedData gLock;
      size_t edgeId;
      if ( mEdgeSeq.getCanonicalForm() == CanonicalForm::REV )
          mEdgeSeq.ReverseComplement();
      if ( true )
      { SpinLocker lock(gLock);
        edgeId = mEdges.size();
        mEdges.push_back(mEdgeSeq); }
      unsigned offset = 0;
      bool err = false;
      for ( BigKmer const* pEnt : mEdgeEntries )
      { if ( pEnt->isUnassigned() )
          const_cast<BigKmer*>(pEnt)->setAssigned();
        else
        { std::cout << edgeId << ':' << offset++ << ' ' << pEnt
                    << " Already occupied." << std::endl;
          err = true; } }
        if ( err )
            FatalErr("Having trouble with preoccupied kmers.");
        mEdgeSeq.clear();
        mEdgeEntries.clear(); }

    BigKDict const& mDict;
    vecbvec& mEdges;
    std::vector<BigKmer const*> mEdgeEntries;
    bvec mEdgeSeq;
};

template <unsigned BIGK>
class Pather
{
public:
    typedef BigDict<BIGK> BigKDict;

    Pather( vecbvec const& reads, BigKDict const& dict, vecbvec const& edges,
                vec<int> const& fwdXlat, vec<int> const& revXlat,
                ReadPathVec* pReadPaths,
                HyperKmerPath const* pHKP, vecKmerPath* pKmerPaths,
                std::vector<bvec const*> const* pReg )
    : mReads(reads), mDict(dict), mEdges(edges),
      mFwdXlat(fwdXlat), mRevXlat(revXlat), mpReadPaths(pReadPaths),
      mpHKP(pHKP), mpKmerPaths(pKmerPaths), mpRegistry(pReg) {}

    void operator()( size_t readId )
    { BigKMer<BIGK>::setRegistry(mpRegistry);
      bvec const& read = mReads[readId];
      if ( read.size() < BIGK ) return;
      KMerHasher<BIGK> hasher;
      BigKMer<BIGK> kmer(read,hasher(read.begin()),uint32_t(readId));
      BigKMer<BIGK> entry = lookup(kmer);
      mTmpReadPath.setFirstSkip(entry.getOffset());
      size_t idx = &entry.getBV()-&mEdges[0];
      size_t edgeId = entry.isRC()?mRevXlat[idx]:mFwdXlat[idx];
      mTmpReadPath.push_back(edgeId);
      size_t readLenRemaining = read.size();
      size_t edgeLenRemaining = entry.getBV().size()-entry.getOffset();
      while ( readLenRemaining > edgeLenRemaining )
      { readLenRemaining = readLenRemaining-edgeLenRemaining+BIGK-1;
        unsigned readOffset = read.size()-readLenRemaining;
        BigKMer<BIGK> nextKmer(read,hasher(read.begin(readOffset)),
                                uint32_t(readId),KMerContext(),readOffset);
        BigKMer<BIGK> nextEntry = lookup(nextKmer);
        ForceAssertEq(nextEntry.getOffset(),0u);
        idx = &nextEntry.getBV()-&mEdges[0];
        edgeId = nextEntry.isRC()?mRevXlat[idx]:mFwdXlat[idx];
        mTmpReadPath.push_back(edgeId);
        edgeLenRemaining = nextEntry.getBV().size(); }
      edgeLenRemaining -= readLenRemaining;
      /* mTmpReadPath.setLastSkip(edgeLenRemaining); */
      if ( mpReadPaths )
        (*mpReadPaths)[readId] = mTmpReadPath;
      if ( mpKmerPaths )
      { buildKmerPath(edgeLenRemaining);
        (*mpKmerPaths)[readId] = mTmpKmerPath; mTmpKmerPath.clear(); }
      mTmpReadPath.clear(); }

private:
    BigKMer<BIGK> lookup( BigKMer<BIGK> const& kmer )
    { if ( kmer.isRev() )
      { BigKMer<BIGK> const* pEnt = mDict.lookup(kmer.rc());
        ForceAssert(pEnt);
        return pEnt->rc(); }
      BigKMer<BIGK> const* pEnt = mDict.lookup(kmer);
      ForceAssert(pEnt);
      return *pEnt; }

    void buildKmerPath( size_t edgeLenRemaining )
    { if ( mTmpReadPath.size() == 1 )
      { KmerPath const& edge = mpHKP->EdgeObject(mTmpReadPath.front());
        mTmpKmerPath.Assign(edge.front().Start()+mTmpReadPath.getFirstSkip(),
                            edge.back().Stop()-edgeLenRemaining); }
      else
      { KmerPath const& edge1 = mpHKP->EdgeObject(mTmpReadPath.front());
        mTmpKmerPath.Assign(edge1.front().Start()+mTmpReadPath.getFirstSkip(),
                            edge1.back().Stop());
        auto end = mTmpReadPath.end()-1;
        for ( auto itr=mTmpReadPath.begin()+1; itr != end; ++itr )
        { KmerPath const& edge = mpHKP->EdgeObject(*itr);
          SerfVec<KmerPathInterval>& tmp = mTmpKmerPath;
          tmp.append(edge.begin(),edge.end()); }
        KmerPath const& edgeN = mpHKP->EdgeObject(*end);
        kmer_id_t stop = edgeN.back().Stop()-edgeLenRemaining;
        mTmpKmerPath.Append(edgeN.front().Start(),stop); } }

    vecbvec const& mReads;
    BigKDict const& mDict;
    vecbvec const& mEdges;
    vec<int> const& mFwdXlat;
    vec<int> const& mRevXlat;
    ReadPathVec* mpReadPaths;
    HyperKmerPath const* mpHKP;
    vecKmerPath* mpKmerPaths;
    std::vector<bvec const*> const* mpRegistry;
    ReadPath mTmpReadPath;
    KmerPath mTmpKmerPath;
};

#if 0
void ValidatePaths( vecbvec const& reads, ReadPathVec const& readPaths,
                        HyperBasevector const& hbv )
{
    auto iPath = readPaths.begin();
    for ( bvec const& read : reads )
    {
        ReadPath const& path = *iPath;
        ++iPath;
        if ( path.empty() )
        {
            std::cout << "Empty path.\n";
            continue;
        }
        auto eBeg = hbv.EdgeObject(path.front()).begin()+path.getFirstSkip();
        auto eEnd = hbv.EdgeObject(path.front()).end();
        auto rBeg = read.begin();
        auto rEnd = read.end();
        size_t dist = std::min(eEnd-eBeg,rEnd-rBeg);
        if ( !std::equal(eBeg,eBeg+dist,rBeg) )
        {
            std::cout << "Bad path compare.\n";
            continue;
        }
        if ( eBeg+dist == eEnd )
        {
            rBeg = rBeg+dist-hbv.K()+1;
        }
        else if ( path.size() > 1 )
        {
            std::cout << "Extra path info.\n";
            continue;
        }
        for ( auto itr=path.begin()+1,end=path.end(); itr != end; ++itr )
        {
            eBeg = hbv.EdgeObject(*itr).begin();
            eEnd = hbv.EdgeObject(*itr).end();
            dist = std::min(eEnd-eBeg,rEnd-rBeg);
            if ( !std::equal(eBeg,eBeg+dist,rBeg) )
            {
                std::cout << "Bad path compare.\n";
                break;
            }
            if ( eBeg+dist == eEnd )
                rBeg = rBeg+dist-hbv.K()+1;
            else if ( itr+1 != end )
            {
                std::cout << "Extra path info.\n";
                break;
            }
        }
    }
}
#endif

template <unsigned BIGK>
void readsToHBV( vecbvec const& reads, unsigned coverage, size_t dictSize,
                    HyperBasevector* pHBV, ReadPathVec* pReadPaths,
                    HyperKmerPath* pHKP, vecKmerPath* pKmerPaths )
{
    if ( pReadPaths )
    {
        pReadPaths->clear();
        pReadPaths->resize(reads.size());
    }
    if ( pKmerPaths )
    {
        pKmerPaths->clear();
        pKmerPaths->resize(reads.size());
    }

    size_t nKmers = reads.getKmerCount(BIGK);
    if ( !nKmers )
    {
        pHBV->Clear();
        if ( pHKP ) pHKP->Clear();
        return;
    }

    // M1 (docs/GraphConstruction.md, section 7): if a dictionary-size hint
    // was supplied, use it -- the hint estimates the number of *distinct*
    // kmers, which is all the dictionary can ever hold.  The legacy
    // sizing nKmers/coverage is based on the total kmer-instance count
    // and overshoots by the redundancy of the input.  An undersized hint
    // is safe (the hash set grows by splitting) but costs rehashing.
    size_t dictSizeToUse = dictSize ?
        ( dictSize < nKmers ? dictSize : nKmers ) : nKmers/coverage;
    if ( !dictSizeToUse ) dictSizeToUse = 1;
    // M10b: compact dictionary entries refer to sequences by registry
    // index: the places first, the edges once the edge set is final.
    std::vector<bvec const*> registry;
    registry.reserve(reads.size()+1000);
    for ( size_t i = 0; i < reads.size(); ++i )
        registry.push_back( &reads[i] );
    // Thread-local registry: every parallel phase sets it on entry (the
    // local gap assemblies run several readsToHV calls concurrently).

    BigDict<BIGK> bigDict( dictSizeToUse,
        typename BigKMer<BIGK>::hasher(), std::equal_to<BigKMer<BIGK>>(),
        TFactory<BigKMer<BIGK>>(), .9 );
    BigKMerizer<BIGK> kmerizer(&bigDict);
    parallelFor(0ul,reads.size(),
            [kmerizer,&reads,&registry]( size_t readId ) mutable
            { BigKMer<BIGK>::setRegistry(&registry);
              kmerizer.kmerize(reads[readId],uint32_t(readId)); });

    vecbvec edges;
    edges.reserve(bigDict.size()/100);
    BigKEdgeBuilder<BIGK>::buildEdges(bigDict,&edges,&registry);
    AssertEq(edges.getKmerCount(BIGK),bigDict.size());
    vec<int> fwdXlat;
    vec<int> revXlat;
    if ( !pReadPaths && !pKmerPaths )
    {
        buildHBVFromEdges(edges,BIGK,pHBV,&fwdXlat,&revXlat);
        if ( pHKP )
            buildHKPFromHBV(*pHBV,fwdXlat,revXlat,pHKP);
    }
    else
    {
        size_t const nPlaces = registry.size();
        for ( size_t i = 0; i < edges.size(); ++i )
            registry.push_back( &edges[i] );
        parallelForBatch(0ul,edges.size(),100,
                [&edges,&kmerizer,nPlaces,&registry]( size_t edgeId )
                { BigKMer<BIGK>::setRegistry(&registry);
                  kmerizer.updateDict(edges[edgeId],
                     uint32_t(nPlaces+edgeId)); });

        buildHBVFromEdges(edges,BIGK,pHBV,&fwdXlat,&revXlat);
        HyperKmerPath hkp;
        if ( pKmerPaths && !pHKP )
            pHKP = &hkp;
        if ( pHKP )
            buildHKPFromHBV(*pHBV,fwdXlat,revXlat,pHKP);

        Pather<BIGK> pather(reads,bigDict,edges,fwdXlat,revXlat,
                                pReadPaths,pHKP,pKmerPaths,&registry);
        parallelForBatch(0ul,reads.size(),10000,pather);
    }
}

template <int BIGK>
struct SillyFunctor
{
    void operator()( vecbvec const& reads, unsigned coverage, size_t dictSize,
            HyperBasevector* pHBV, ReadPathVec* pReadPaths,
            HyperKmerPath* pHKV, vecKmerPath* pKmerPaths )
    { readsToHBV<BIGK>(reads,coverage,dictSize,pHBV,pReadPaths,pHKV,pKmerPaths); }
};

}

void buildBigKHBVFromReads( unsigned K, vecbvec const& reads, unsigned coverage,
                            HyperBasevector* pHBV, ReadPathVec* pReadPaths,
                            HyperKmerPath* pHKP, vecKmerPath* pKmerPaths,
                            size_t dictSize )
{
    BigK::dispatch<SillyFunctor>(K,reads,coverage,dictSize,
                                    pHBV,pReadPaths,pHKP,pKmerPaths);
}
