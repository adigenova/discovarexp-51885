///////////////////////////////////////////////////////////////////////////////

#ifndef KMER_DICT_PACKED_H_
#define KMER_DICT_PACKED_H_

// KmerDictPacked: a memory-compact replacement for the feudal
// KmerDict<K> on the DiscovarExp kmerize path
// (docs/GraphConstruction.md, items M9/M10).  Storage is one contiguous
// arena of KmerDictEntry<K>, partitioned into power-of-two buckets by a
// cheap hash of the kmer words and sorted by kmer within each bucket,
// with a prefix-offset directory.  It provides the method surface that
// the kmerize / unitig-walk / read-pathing code uses on KmerDict, plus a
// bulk-build interface fed by the M9 disk-bin qualified stream.  Legacy
// code paths keep the original feudal KmerDict.
//
// Lookup cost: a bucket holds ~4K entries (~80 KB at K=60), so a probe is
// a dozen or so cache-friendly comparisons of an in-L2 range.
//
// Thread safety: findEntry/findEntryCanonical are safe for concurrent
// readers.  applyCanonical/insertEntry take a mutex (they are only used
// by the legacy fill-joins/process paths); the bulk-build interface must
// be driven from a single thread, and finalizeBulk parallelizes the
// per-bucket sorts internally.

#include "Basevector.h"
#include "kmers/KMer.h"
#include "kmers/KMerContext.h"
#include "kmers/ReadPather.h"
#include "system/WorklistN.h"
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <mutex>
#include <vector>

template <unsigned K>
class KmerDictPacked
{
public:
    typedef KmerDictEntry<K> Entry;
    typedef unsigned long word_t;
    static unsigned const WORDS = (K+31)/32;   // storage units of KMer<K>

    explicit KmerDictPacked( unsigned nBuckets = 1u )
    : mpE(0), mSize(0), mBuilt(false)
    {    setBuckets(nBuckets);    }

    ~KmerDictPacked() { delete[] mpE; }

    KmerDictPacked( KmerDictPacked const& ) = delete;
    KmerDictPacked& operator=( KmerDictPacked const& ) = delete;

    size_t size() const { return mSize; }
    unsigned nBuckets() const { return mBuckets; }

    // Cheap, deterministic hash of the kmer words.  The M9 binning code
    // uses the same hash (modulo the number of disk bins).
    static word_t hashWords( word_t const* words )
    {    word_t h = 14695981039346656037ul;
         for ( unsigned i = 0; i < WORDS; i++ )
         {    h ^= words[i];
              h *= 1099511628211ul;    }
         return h;    }

    static void wordsOf( KMer<K> const& kmer, word_t* words )
    {    std::memcpy( words, &kmer, WORDS*sizeof(word_t) );    }

    // ---------------- bulk build (M9 qualified stream) ----------------

    // bucketCounts must hold the exact per-bucket totals of the stream.
    void beginBulk( size_t nEntries, std::vector<size_t> const& bucketCounts )
    {    ForceAssert( !mBuilt && !mpE );
         ForceAssertEq( bucketCounts.size(), mBuckets );
         mSize = nEntries;
         mpE = mSize ? new Entry[mSize] : 0;
         mStart.resize(mBuckets+1);
         mStart[0] = 0;
         for ( unsigned b = 0; b < mBuckets; b++ )
              mStart[b+1] = mStart[b] + bucketCounts[b];
         ForceAssertEq( mStart[mBuckets], mSize );
         mCursor.assign( mStart.begin(), mStart.end()-1 );   }

    void placeBulk( unsigned bucket, word_t const* words,
                    unsigned char ctx, unsigned count )
    {    ForceAssert( !mBuilt );
         size_t& cur = mCursor[bucket];
         ForceAssertLt( cur, mStart[bucket+1] );
         Entry& e = mpE[cur++];
         KMer<K> kmer;
         std::memcpy( &kmer, words, WORDS*sizeof(word_t) );
         e = Entry( kmer, contextFromByte(ctx) );
         e.getKDef().setCount(count);   }

    void finalizeBulk()
    {    ForceAssert( !mBuilt );
         mCursor.clear();
         #pragma omp parallel for schedule(dynamic,64)
         for ( int64_t b = 0; b < int64_t(mBuckets); b++ )
              std::sort( mpE+mStart[b], mpE+mStart[b+1], EntryCmp() );
         // (loop var is int64_t: needs <cstdint>, included above)
         mBuilt = true;   }

    // ---------------- lookups ----------------

    // kmer must be in canonical (stored) form.
    Entry const* findEntryCanonical( KMer<K> const& kmer ) const
    {    word_t w[WORDS];
         wordsOf( kmer, w );
         return findWords(w);    }

    // canonicalizes REV forms; FWD/PALINDROME are used as-is.
    Entry const* findEntry( KMer<K> const& kmer ) const
    {    if ( kmer.getCanonicalForm() == CanonicalForm::REV )
              return findEntryCanonical( KMer<K>(kmer).rc() );
         return findEntryCanonical(kmer);    }

    // ---------------- context updates ----------------

    // HashSet::apply semantics: the kmer is expected in canonical form;
    // if absent it is inserted (with a default KDef), then func applied.
    // Serialized by a mutex for the legacy parallel paths.
    template <class Func>
    void applyCanonical( KMer<K> const& kmer, Func const& func )
    {    std::lock_guard<std::mutex> lock(mMutex);
         word_t w[WORDS];
         wordsOf( kmer, w );
         Entry const* pEnt = findWords(w);
         if ( !pEnt )
         {    unsigned b = hashWords(w) & (mBuckets-1);
              mSpill[b].push_back( Entry(kmer, KMerContext()) );
              pEnt = &mSpill[b].back();    }
         func(*pEnt);    }

    // ---------------- legacy inserts ----------------

    // Old insertEntry assumed novelty; we merge defensively (nobody can
    // legitimately rely on duplicate entries).
    void insertEntry( Entry const& entry )
    {    std::lock_guard<std::mutex> lock(mMutex);
         word_t w[WORDS];
         std::memcpy( w, &entry, WORDS*sizeof(word_t) );
         Entry const* pEnt = findWords(w);
         if ( pEnt )
         {    KDef& kDef = const_cast<KDef&>(pEnt->getKDef());
              kDef.orContext( entry.getKDef().getContext() );    }
         else
         {    unsigned b = hashWords(w) & (mBuckets-1);
              mSpill[b].push_back(entry);    }    }

    // ---------------- iteration ----------------

    // Range type passed to parallelForEachHHS procs.
    class ChunkRange
    {    public:
          ChunkRange( Entry const* beg, Entry const* end ) : mB(beg), mE(end) {}
          Entry const* begin() const { return mB; }
          Entry const* end() const { return mE; }
         private:
          Entry const* mB;
          Entry const* mE;    };
    typedef ChunkRange HHS;

    template <class Proc>
    void parallelForEachHHS( Proc proc ) const
    {    std::vector<ChunkRange> ranges;
         if ( mSize )
         {    size_t nChunks = 1 + mSize/262144;
              size_t chunk = 1 + mSize/nChunks;
              for ( size_t i = 0; i < mSize; i += chunk )
                   ranges.push_back( ChunkRange( mpE+i,
                        mpE+std::min(i+chunk,mSize) ) );    }
         for ( unsigned b = 0; b < mBuckets; b++ )
              if ( !mSpill[b].empty() )
                   ranges.push_back( ChunkRange( mSpill[b].data(),
                        mSpill[b].data()+mSpill[b].size() ) );
         if ( ranges.empty() ) return;
         parallelFor( 0ul, ranges.size(),
              [&ranges,proc]( size_t i ) mutable
              {    proc( ranges[i] );    } );    }

    // HHS-style whole-dict iteration (the legacy loops do a nested
    // for over "hash sets" then entries; we present one pseudo-set).
    class AllIter
    {    public:
          AllIter( KmerDictPacked const* dict, bool atEnd );
          Entry const& operator*() const;
          AllIter& operator++();
          bool operator!=( AllIter const& that ) const;
         private:
          void settle();        // initial placement
          void advanceSpill();  // next spill entry, or exhaust
          KmerDictPacked const* mpDict;
          Entry const* mpPos;          // position in arena
          unsigned mSpillBucket;
          size_t mSpillPos;            };

    class AllRange
    {    public:
          AllRange( KmerDictPacked const* dict ) : mpDict(dict) {}
          AllIter begin() const { return AllIter(mpDict,false); }
          AllIter end() const { return AllIter(mpDict,true); }
         private:
          KmerDictPacked const* mpDict;    };

    class HHSIter
    {    public:
          explicit HHSIter( KmerDictPacked const* dict ) : mpDict(dict) {}
          AllRange operator*() const { return AllRange(mpDict); }
          HHSIter& operator++() { mpDict = 0; return *this; }
          bool operator!=( HHSIter const& that ) const
          {    return mpDict != that.mpDict;    }
         private:
          KmerDictPacked const* mpDict;    };

    HHSIter begin() const { return HHSIter(this); }
    HHSIter end() const { return HHSIter(0); }

    // ---------------- maintenance ----------------

    void nullEntries()
    {    parallelForEachHHS(
              []( ChunkRange const& hhs )
              {    for ( Entry const& entry : hhs )
                       const_cast<Entry&>(entry).getKDef().setNull(); } );    }

    void recomputeAdjacencies()
    {    parallelForEachHHS( AdjProc(*this) );    }

    // ---------------- legacy: kmerize a read set ----------------
    // Single-threaded, safe, and only used by legacy paths.

    void process( vecbvec const& reads, bool verbose=false,
                   bool validate=false, unsigned nThreads=0,
                   size_t batchSize=10000 )
    {    for ( auto itr = reads.begin(), end = reads.end(); itr != end; ++itr )
         {    bvec const& read = *itr;
              if ( read.size() < K+1 ) continue;
              auto beg = read.begin(), itr2 = beg+K, last = beg+(read.size()-1);
              KMer<K> kkk(beg);
              KMerContext kc = KMerContext::initialContext(*itr2);
              feed( kkk, kc );
              while ( itr2 != last )
              {    unsigned char pred = kkk.front();
                   kkk.toSuccessor(*itr2); ++itr2;
                   kc = KMerContext(pred,*itr2);
                   feed( kkk, kc );    }
              kc = KMerContext::finalContext(kkk.front());
              kkk.toSuccessor(*last);
              feed( kkk, kc );    }    }

private:
    class EntryCmp
    {    public:
          bool operator()( Entry const& a, Entry const& b ) const
          {    return compare( a, b ) < 0;    }    };

    class AdjProc
    {    public:
          AdjProc( KmerDictPacked const& dict ) : mDict(dict) {}
          void operator()( ChunkRange const& hhs ) const
          {    for ( Entry const& entry : hhs )
               {    KDef& kDef = const_cast<KDef&>(entry.getKDef());
                    KMerContext context = kDef.getContext();
                    if ( context.getSuccessors() )
                    {    KMer<K> kmer(entry);
                         kmer.toSuccessor(0);
                         for ( unsigned succCode = 0; succCode < 4u; ++succCode )
                         {    if ( context.isSuccessor(succCode) )
                              {    kmer.setBack(succCode);
                                   if ( !mDict.findEntry(kmer) )
                                        context.removeSuccessor(succCode); } } }
                    if ( context.getPredecessors() )
                    {    KMer<K> kmer(entry);
                         kmer.toPredecessor(0);
                         for ( unsigned predCode = 0; predCode < 4u; ++predCode )
                         {    if ( context.isPredecessor(predCode) )
                              {    kmer.setFront(predCode);
                                   if ( !mDict.findEntry(kmer) )
                                        context.removePredecessor(predCode); } } }
                    kDef.setContext(context);    }    }
         private:
          KmerDictPacked const& mDict;    };

    void feed( KMer<K> const& kmer, KMerContext context )
    {    switch ( kmer.getCanonicalForm() )
         {    case CanonicalForm::FWD:
                  applyCanonical( kmer, Functor(context) );
                  break;
              case CanonicalForm::REV:
                  applyCanonical( KMer<K>(kmer).rc(), Functor(context.rc()) );
                  break;
              case CanonicalForm::PALINDROME:
                  applyCanonical( kmer, Functor(context|=context.rc()) );
                  break;    }    }

    class Functor
    {    public:
          Functor( KMerContext kc ) : mKC(kc) {}
          void operator()( Entry const& entry ) const
          {    KDef& kDef = const_cast<KDef&>(entry.getKDef());
               kDef.setContext( mKC | kDef.getContext() );
               kDef.incrementCount();    }
     private:
          KMerContext mKC;    };

    // Byte layout (shared with the M9 binning code): predecessors occupy
    // the high nibble, successors the low nibble -- exactly the KMerContext
    // bit layout (pred bit 4+code, succ bit code).
    static KMerContext contextFromByte( unsigned char b )
    {    KMerContext kc;
         for ( unsigned code = 0; code < 4u; code++ )
         {    if ( (b >> (4+code)) & 1 ) kc.setPredecessor(code);
              if ( (b >> code) & 1 ) kc.setSuccessor(code);    }
         return kc;    }

    void setBuckets( unsigned nBuckets )
    {    unsigned b = 1;
         while ( b < nBuckets ) b <<= 1;
         mBuckets = b;
         mStart.resize(mBuckets+1,0);
         mSpill.resize(mBuckets);    }

    Entry const* findWords( word_t const* w ) const
    {    unsigned b = hashWords(w) & (mBuckets-1);
         Entry const* beg = mpE + mStart[b];
         Entry const* end = mpE + mStart[b+1];
         KMer<K> kmer;
         std::memcpy( &kmer, w, WORDS*sizeof(word_t) );
         while ( beg < end )
         {    Entry const* mid = beg + (end-beg)/2;
              if ( compare( *mid, kmer ) < 0 ) beg = mid+1;
              else end = mid;    }
         if ( beg < mpE + mStart[b+1] && !compare( *beg, kmer ) )
              return beg;
         for ( Entry const& e : mSpill[b] )
              if ( !compare( e, kmer ) ) return &e;
         return 0;    }

    Entry* mpE;
    size_t mSize;
    unsigned mBuckets;
    std::vector<size_t> mStart;      // mBuckets+1 prefix offsets
    std::vector<size_t> mCursor;     // bulk placement cursors
    std::vector<std::vector<Entry>> mSpill;   // post-build inserts
    std::mutex mMutex;
    bool mBuilt;
};

template <unsigned K>
KmerDictPacked<K>::AllIter::AllIter( KmerDictPacked<K> const* dict, bool atEnd )
: mpDict( atEnd ? 0 : dict ), mpPos(0), mSpillBucket(0), mSpillPos(0)
{    if ( !atEnd ) settle();    }

template <unsigned K>
void KmerDictPacked<K>::AllIter::settle()
{    if ( mpDict->mSize )
     {    mpPos = mpDict->mpE;
          return;    }
     advanceSpill();    }

template <unsigned K>
void KmerDictPacked<K>::AllIter::advanceSpill()
{    while ( mSpillBucket < mpDict->mBuckets )
     {    if ( mSpillPos < mpDict->mSpill[mSpillBucket].size() ) return;
          mSpillBucket++;
          mSpillPos = 0;    }
     mpDict = 0;    }

template <unsigned K>
typename KmerDictPacked<K>::Entry const&
KmerDictPacked<K>::AllIter::operator*() const
{    if ( mpPos ) return *mpPos;
     return mpDict->mSpill[mSpillBucket][mSpillPos];    }

template <unsigned K>
typename KmerDictPacked<K>::AllIter&
KmerDictPacked<K>::AllIter::operator++()
{    if ( mpPos )
         {    if ( ++mpPos == mpDict->mpE + mpDict->mSize )
              {    mpPos = 0;
                   mSpillBucket = 0;
                   mSpillPos = 0;
                   advanceSpill();    }
              return *this;    }
     ++mSpillPos;
     if ( mSpillPos == mpDict->mSpill[mSpillBucket].size() )
     {    ++mSpillBucket;
          mSpillPos = 0;
          advanceSpill();    }
     return *this;    }

template <unsigned K>
bool KmerDictPacked<K>::AllIter::operator!=( AllIter const& that ) const
{    // The end iterator and an exhausted iterator both carry a null dict
    // pointer; a live iterator never does.  (An exhausted iterator may
    // hold a stale spill-bucket index, so only mpDict is compared.)
    return mpDict != that.mpDict;    }

#endif /* KMER_DICT_PACKED_H_ */
