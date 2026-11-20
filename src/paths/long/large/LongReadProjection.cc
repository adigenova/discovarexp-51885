///////////////////////////////////////////////////////////////////////////////

#include "paths/long/large/LongReadProjection.h"

#include "FastIfstream.h"
#include "MainTools.h"
#include "dna/Bases.h"
#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

namespace {

// One index entry: k-mer code CODE, occurring at offset POS of edge EDGE.
struct KEntry
{
     uint64_t code;
     uint32_t edge;
     uint32_t pos;
};

bool operator<( KEntry const& a, KEntry const& b )
{    if ( a.code < b.code ) return true;
     if ( a.code > b.code ) return false;
     if ( a.edge < b.edge ) return true;
     if ( a.edge > b.edge ) return false;
     return a.pos < b.pos;    }

// One anchor of a read: k-mer at read offset R matches edge EDGE at
// edge offset O.
struct Anchor
{
     uint32_t r;
     uint32_t e;
     uint32_t o;
};

inline int Slop( int dist )
{    return 250 + dist/30;    }

inline bool ValidBase( char c )
{    return c=='A' || c=='C' || c=='G' || c=='T'
          || c=='a' || c=='c' || c=='g' || c=='t';    }

} // end of anonymous namespace

void ProjectLongReads( const HyperBasevector& hb,
     const vecbasevector& edges, const String& reads_file,
     const unsigned K, const int min_read, const int anchor_step,
     const int max_anchors, const String& out_head, const Bool VERBOSE )
{
     ForceAssertLe( K, 32u );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     int const OVER = hb.K( ) - 1;   // de Bruijn overlap at vertices

     // Build the k-mer index over all (doubled) edge sequences.

     cout << Date( ) << ": indexing " << K << "-mers of "
          << edges.size( ) << " edge sequences" << endl;
     size_t nkm = 0;
     for ( size_t e = 0; e < edges.size( ); e++ )
          if ( edges[e].size( ) >= K ) nkm += edges[e].size( ) - K + 1;
     std::vector<KEntry> index;
     index.reserve(nkm);
     {    uint64_t const mask = ( uint64_t(1) << (2*K) ) - 1;
          for ( size_t e = 0; e < edges.size( ); e++ )
          {    bvec const& bv = edges[e];
               if ( bv.size( ) < K ) continue;
               uint64_t code = 0;
               for ( size_t j = 0; j < K; j++ )
                    code = ( code << 2 ) | bv[j];
               index.push_back( { code, uint32_t(e), 0 } );
               for ( size_t j = K; j < bv.size( ); j++ )
               {    code = ( ( code << 2 ) | bv[j] ) & mask;
                    index.push_back(
                         { code, uint32_t(e), uint32_t(j-K+1) } );    }    }    }
     std::sort( index.begin( ), index.end( ) );
     cout << Date( ) << ": index built, " << ToStringAddCommas(index.size( ))
          << " entries" << endl;

     // Stream the reads (fasta/fastq, possibly gzipped) and project.

     String cmd = "cat " + reads_file;
     if ( reads_file.Contains( ".gz", -1 ) ) cmd = "zcat " + reads_file;
     fast_pipe_ifstream in(cmd);

     int64_t n_reads = 0, total_bp = 0, n_ge = 0;
     int64_t n_ge_2edges = 0, n_ge_3edges = 0, n_anchored_bp = 0;
     int64_t n_ge_span_bp = 0;
     vec<String> report_lines;
     String line, header, raw, seq;

     auto process = [&]( String const& header, String const& raw )
     {    int64_t id = n_reads++;
          total_bp += raw.size( );
          if ( (int64_t) raw.size( ) < min_read ) return;
          n_ge++;

          // Sanitize: 2-bit storage can only hold ACGT.

          seq.resize(0);
          seq.reserve(raw.size( ));
          for ( int i = 0; i < raw.isize( ); i++ )
               seq.push_back( ValidBase(raw[i]) ? raw[i] : 'A' );
          bvec rd(seq);

          // Seed: exact k-mer matches, roughly every anchor_step bases.

          vec<Anchor> anchors;
          int next_anchor_pos = 0;
          uint64_t const mask = ( uint64_t(1) << (2*K) ) - 1;
          uint64_t code = 0;
          for ( size_t j = 0; j < rd.size( ); j++ )
          {    code = ( ( code << 2 ) | rd[j] ) & mask;
               if ( j+1 < K ) continue;
               int pos = int(j) - int(K) + 1;
               if ( pos < next_anchor_pos ) continue;
               auto lo = std::lower_bound( index.begin(), index.end(),
                    KEntry{ code, 0, 0 },
                    []( KEntry const& a, KEntry const& b )
                    {    return a.code < b.code;    } );
               auto hi = std::upper_bound( index.begin(), index.end(),
                    KEntry{ code, 0, 0 },
                    []( KEntry const& a, KEntry const& b )
                    {    return a.code < b.code;    } );
               if ( lo == hi ) continue;
               if ( hi - lo > 16 ) continue;   // pathological repeat
               for ( auto it = lo; it != hi; ++it )
                    anchors.push_back( { uint32_t(pos), it->edge,
                         it->pos } );
               next_anchor_pos = pos + anchor_step;    }

          // Cap the anchors per read (even subsample).

          if ( anchors.isize( ) > max_anchors )
          {    int stride = 1 + anchors.isize( ) / max_anchors;
               vec<Anchor> keep;
               for ( int i = 0; i < anchors.isize( ); i += stride )
                    keep.push_back(anchors[i]);
               anchors = keep;    }

          // Longest consistent chain (DP over anchors in read order).

          int n = anchors.isize( );
          vec<int> best( n, 1 ), prev( n, -1 );
          auto compatible = [&]( Anchor const& a, Anchor const& b )
          {    int dist = int(b.r) - int(a.r);
               if ( dist <= 0 ) return false;
               if ( a.e == b.e )
               {    int diff = int(b.o) - int(a.o);
                    return diff > 0 && abs( diff - dist ) <= Slop(dist);    }
               if ( to_right[ a.e ] != to_left[ b.e ] ) return false;
               // Junction algebra: edges overlap by K-1 at the vertex.
               int len_a = edges[ a.e ].isize( );
               int expected = int(a.o) + dist - len_a + OVER;
               if ( expected < 0 ) return false;
               if ( expected + int(K) > edges[ b.e ].isize( ) ) return false;
               return abs( int(b.o) - expected ) <= Slop(dist);    };
          int best_end = -1, best_len = 0;
          for ( int i = 0; i < n; i++ )
          {    for ( int j = 0; j < i; j++ )
                    if ( compatible( anchors[j], anchors[i] )
                         && best[j] + 1 > best[i] )
                    {    best[i] = best[j] + 1;
                         prev[i] = j;    }
               if ( best[i] > best_len )
               {    best_len = best[i];
                    best_end = i;    }    }

          // Reconstruct, merge collinear anchors into segments.

          vec<Anchor> chain;
          for ( int i = best_end; i >= 0; i = prev[i] )
               chain.push_back(anchors[i]);
          Reverse(chain);
          std::sort( chain.begin( ), chain.end( ),
               []( Anchor const& a, Anchor const& b )
               {    return a.r < b.r;    } );
          vec<LRSeg> segs;
          for ( int i = 0; i < chain.isize( ); i++ )
          {    if ( i > 0 && chain[i].e == uint32_t(segs.back( ).edge) )
               {    segs.back( ).o2 = int(chain[i].o);
                    segs.back( ).r2 = int(chain[i].r);    }
               else
                    segs.push_back( { int(chain[i].e), int(chain[i].o),
                         int(chain[i].o), int(chain[i].r),
                         int(chain[i].r) } );    }
          for ( auto& s : segs ) { s.r2 += int(K); s.o2 += int(K); }

          // Stats.

          int n_distinct = 0;
          {    vec<int> used;
               for ( auto const& s : segs ) used.push_back(s.edge);
               Unique(used);
               n_distinct = used.isize( );    }
          if ( n_distinct >= 2 ) n_ge_2edges++;
          if ( n_distinct >= 3 ) n_ge_3edges++;
          if ( !segs.empty( ) )
          {    int64_t span = segs.back( ).r2 - segs.front( ).r1;
               if ( span < 0 ) span = -span;
               n_anchored_bp += span;
               n_ge_span_bp += rd.size( );    }

          // Report line.

          ostringstream ss;
          ss << id << " " << rd.size( ) << " " << n << " "
             << segs.size( ) << " " << n_distinct << " : ";
          for ( auto const& s : segs )
               ss << s.edge << "[" << s.o1 << "-" << s.o2 << "@"
                  << s.r1 << "-" << s.r2 << "] ";
          report_lines.push_back(ss.str());

          if ( VERBOSE && n_distinct >= 2 && id < 20 )
               cout << "read " << id << " (" << rd.size( ) << " bp, "
                    << n << " anchors): " << ss.str( ) << endl;    };

     // Parse the stream: fasta (multi-line) or fastq.

     Bool fastq = False;
     Bool have_header = False;
     while ( true )
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( line.empty( ) ) continue;
          if ( !have_header )
          {    if ( line[0] == '>' ) fastq = False;
               else if ( line[0] == '@' ) fastq = True;
               else FatalErr( "Unrecognized read file format." );
               have_header = True;    }
          if ( fastq )
          {    if ( line[0] != '@' ) continue;
               header = line;
               getline( in, raw );
               getline( in, line );   // '+'
               getline( in, line );   // quals
               if ( in.fail( ) ) break;
               process( header, raw );    }
          else
          {    if ( line[0] == '>' )
               {    if ( !header.empty( ) ) process( header, raw );
                    header = line;
                    raw.resize(0);    }
               else
                    raw += line;    }    }
     if ( !fastq && !header.empty( ) ) process( header, raw );

     // Write projections (text, for eyeballing and for P3 to consume).

     {    Ofstream( tout, out_head + ".proj.txt" );
          for ( auto const& r : report_lines )
               tout << r << "\n";    }

     // Report acceptance statistics.

     cout << "\n===================================================="
          << "=======================\n";
     cout << "P2/LR-1 projection stats (" << reads_file << ")\n";
     cout << "  reads: " << n_reads
          << ", total bp: " << ToStringAddCommas(total_bp) << "\n";
     cout << "  reads >= " << min_read << " bp: " << n_ge << "\n";
     if ( n_ge )
     {    cout << fixed << setprecision(1)
               << "  projecting onto >= 2 edges: " << n_ge_2edges
               << " (" << 100. * n_ge_2edges / n_ge << "%)"
               << "   [acceptance >= 95%]\n"
               << "  projecting onto >= 3 edges: " << n_ge_3edges
               << " (" << 100. * n_ge_3edges / n_ge << "%)\n";
          if ( n_ge_span_bp )
               cout << "  mean chain span / read length: "
                    << 100. * double(n_anchored_bp) / double(n_ge_span_bp)
                    << "%\n";    }
     cout << "=========================================================="
          << "=============\n" << resetiosflags(ios::fixed);
}
