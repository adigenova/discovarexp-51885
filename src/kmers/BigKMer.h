/////////////////////////////////////////////////////////////////////////////

#ifndef BIGKMER_H_
#define BIGKMER_H_

#include "Basevector.h"
#include "dna/CanonicalForm.h"
#include "feudal/Iterator.h"
#include "kmers/KMerContext.h"
#include <cstdint>
#include <ostream>
#include <vector>

// BigKMer: a view of a K-mer in a set of sequences, in a compact 16-byte
// layout (docs/GraphConstruction.md, item M10b).  The entry no longer
// stores a bvec pointer; it stores a 32-bit index into a sequence
// *registry* (a vector of bvec pointers covering the place sequences,
// then the edge sequences, set up by readsToHV in BigKPather.cc).  A
// truncated 32-bit FNV hash is kept as a fast inequality filter, and
// doubles as the in-bucket discriminator for the hash set.
//
// Temporary kmers built from rolling scratch buffers (BigKEdgeBuilder's
// mEdgeSeq) are not registry members; they carry the DIRECT_SEQ sentinel
// and resolve through a thread-local scratch pointer set by the edge
// builder before each walk.

template <unsigned K>
class BigKMer
{
public:
    typedef unsigned char Call;
    struct hasher
    {
        typedef BigKMer argument_type;
        size_t operator()( BigKMer const& bigKMer ) const
        { return bigKMer.mHash; }
    };

    class iterator
    :  public std::iterator<std::random_access_iterator_tag,Call,int,void,Call>,
       public IteratorBase<iterator,unsigned,int>
    {
    public:
        iterator( BigKMer const* pBigKMer, unsigned offset )
        : IteratorBase<iterator,unsigned,int>(offset), mpBigKMer(pBigKMer) {}

        Call operator*() const
        { return (*mpBigKMer)[this->pos()]; }

        Call operator[]( int idx ) const
        { return (*mpBigKMer)[this->pos()+idx]; }

    private:
        BigKMer const* mpBigKMer;
    };

    BigKMer() : mSeqIdx(DIRECT_SEQ), mOffset(0), mHash(0), mRC(false),
                mAssigned(false), mContext()
    {}

    // View of a registry sequence: bv must be registry[seqIdx].
    BigKMer( bvec const& bv, size_t hash, uint32_t seqIdx,
                KMerContext context=KMerContext(), unsigned offset = 0 )
    : mSeqIdx(seqIdx), mOffset(offset), mHash(uint32_t(hash)),
      mRC(false), mAssigned(false), mContext(context)
    {}

    // View of the thread-local scratch sequence (see setScratch).
    BigKMer( bvec const& bv, size_t hash, KMerContext context=KMerContext() )
    : mSeqIdx(DIRECT_SEQ), mOffset(0), mHash(uint32_t(hash)),
      mRC(false), mAssigned(false), mContext(context)
    {}

    bvec const& getBV() const { return *getSeq(); }
    unsigned getOffset() const { return mOffset; }
    bool isRC() const { return mRC; }

    KMerContext getContext() const { return mContext; }
    void addContext( KMerContext context ) { mContext |= context; }

    bool isUnassigned() const { return !mAssigned; }
    void setAssigned() { mAssigned = true; }
    void setUnassigned() { mAssigned = false; }

    iterator begin() const { return iterator(this,0); }
    iterator end() const { return iterator(this,K); }

    Call operator[]( unsigned offset ) const
    { offset += mOffset;
      return mRC ? (*getSeq())[getSeq()->size()-offset-1]^3
                 : (*getSeq())[offset]; }

    Call front() const { return (*this)[0]; }
    Call back() const { return (*this)[K-1]; }

    // Sequence registry: thread-local, set on entry to every parallel
    // phase of readsToHV (see BigKPather.cc).  Thread-local because the
    // local gap assemblies run several readsToHV calls concurrently.
    static void setRegistry( std::vector<bvec const*> const* pReg )
    { gpRegistry = pReg; }
    static std::vector<bvec const*> const* getRegistry() { return gpRegistry; }

    // Thread-local scratch sequence for non-registry temporaries.
    static void setScratch( bvec const* pBV ) { gpScratch = pBV; }
    static bvec const* getScratch() { return gpScratch; }

    void updateLocation( BigKMer const& kmer )
    { mSeqIdx = kmer.mSeqIdx; mOffset = kmer.mOffset; mRC = kmer.mRC;
      AssertEq(mHash,kmer.mHash); }

    BigKMer rc() const
    { unsigned offset = getSeq()->size()-mOffset-K;
      return BigKMer( mSeqIdx, offset, !mRC, mAssigned, mContext.rc(),
                     mHash ); }

    CanonicalForm getForm() const
    { return CF<K>::getForm(begin()); }

    bool isFwd() const { return getForm()==CanonicalForm::FWD; }
    bool isRev() const { return getForm()==CanonicalForm::REV; }
    bool isPalindrome() const { return getForm()==CanonicalForm::PALINDROME; }

    BigKMer& successor( size_t hash )
    { ++mOffset; mHash = uint32_t(hash); return *this; }

    BigKMer& successor( size_t hash, KMerContext context )
    { successor(hash); mContext = context; return *this; }

    friend bool identical( BigKMer const& bigKMer1, BigKMer const& bigKMer2 )
    { return bigKMer1.mSeqIdx == bigKMer2.mSeqIdx &&
             bigKMer1.mOffset == bigKMer2.mOffset &&
             bigKMer1.mRC == bigKMer2.mRC; }

    friend bool operator==( BigKMer const& bigKMer1, BigKMer const& bigKMer2 )
    { return bigKMer1.mHash==bigKMer2.mHash &&
            (identical(bigKMer1,bigKMer2) ||
             std::equal(bigKMer1.begin(),bigKMer1.end(),bigKMer2.begin())); }

    friend bool operator<( BigKMer const& bigKMer1, BigKMer const& bigKMer2 )
    { auto itr1=bigKMer1.begin(), end=bigKMer1.end();
      for ( auto itr2=bigKMer2.begin(); itr1 != end; ++itr1, ++itr2 )
      { if ( *itr1 < *itr2 ) return true;
        if ( *itr1 > *itr2 ) return false; }
      return false; }

    friend ostream& operator<<( ostream& os, BigKMer const& kmer )
    { PrintBasesIter(os,kmer.begin(),kmer.end(),K+1); return os; }

private:
    BigKMer( uint32_t seqIdx, unsigned offset, bool rc, bool assigned,
                KMerContext context, uint32_t hash )
    : mSeqIdx(seqIdx), mOffset(offset), mHash(hash), mRC(rc),
      mAssigned(assigned), mContext(context)
    {}

    bvec const* getSeq() const
    { return mSeqIdx==DIRECT_SEQ ? gpScratch : (*gpRegistry)[mSeqIdx]; }

    static uint32_t const DIRECT_SEQ = ~uint32_t(0);
    static thread_local std::vector<bvec const*> const* gpRegistry;
    static thread_local bvec const* gpScratch;

    uint32_t mSeqIdx;
    uint32_t mOffset;
    uint32_t mHash;
    bool mRC;
    bool mAssigned;
    KMerContext mContext;
};

template <unsigned K>
thread_local std::vector<bvec const*> const* BigKMer<K>::gpRegistry = 0;

template <unsigned K>
thread_local bvec const* BigKMer<K>::gpScratch = 0;

#endif /* BIGKMER_H_ */
