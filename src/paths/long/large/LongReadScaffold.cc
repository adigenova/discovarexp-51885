///////////////////////////////////////////////////////////////////////////////

#include "paths/long/large/LongReadScaffold.h"

#include "FastIfstream.h"
#include "MainTools.h"
#include "dna/Bases.h"
#include "paths/long/ultra/MultipleAligner.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <cstdint>
#include <deque>
#include <map>
#include <set>
#include <vector>

namespace {

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

// Edge classification for scaffolding.
enum EdgeClass { EDGE_NORMAL, EDGE_REPEAT, EDGE_SHORT };

} // end of anonymous namespace

void ScaffoldWithLongReads( const HyperBasevector& hb,
     const vec<int>& inv, const vecbasevector& edges,
     const vec<double>& edge_cov,
     const String& reads_file,
     const unsigned K, const int min_read, const int min_support,
     const int anchor_step, const double repeat_factor,
     const int min_line, const int max_bfs_depth,
     const int max_mult, const int min_anchor_block,
     const double min_chain_cov,
     const String& mapper, const String& mm2_bin,
     const int min_map_len, const double min_map_id,
     const int min_seg, const String& debug_read,
     const String& out_head, const Bool VERBOSE )
{
     ForceAssertLe( K, 32u );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     int const OVER = hb.K( ) - 1;

     // ---- 1. Classify edges ----

     // Compute median coverage of canonical edges (for repeat threshold).
     vec<double> canon_cov;
     for ( int e = 0; e < hb.E( ); e++ )
          if ( e <= inv[e] )   // canonical (or palindromic)
               if ( edge_cov[e] > 0 ) canon_cov.push_back(edge_cov[e]);
     Sort(canon_cov);
     double median_cov = canon_cov.empty( ) ? 1.0
          : canon_cov[ canon_cov.size( ) / 2 ];
     double repeat_threshold = median_cov * repeat_factor;

     vec<int> edge_class( hb.E( ), (int)EDGE_NORMAL );
     int n_repeat = 0, n_short = 0;
     for ( int e = 0; e < hb.E( ); e++ )
     {    if ( edges[e].isize( ) < min_line )
          {    edge_class[e] = EDGE_SHORT;    n_short++;    }
          else if ( edge_cov[e] > repeat_threshold )
          {    edge_class[e] = EDGE_REPEAT;  n_repeat++;    }    }

     cout << Date( ) << ": edge classification: " << hb.E( ) << " edges, "
          << n_repeat << " repeats (cov > " << repeat_threshold
          << "), " << n_short << " short (< " << min_line
          << " bp), " << hb.E( ) - n_repeat - n_short
          << " normal" << endl;

     // ---- 3. Map reads onto the graph: minimap2 or internal seeds ----

     // Shared per-read bookkeeping (both mappers).
     int64_t n_reads = 0, total_bp = 0, n_ge = 0;
     int64_t n_hiq_reads = 0;
     // Per-link: gap sum, gap^2 sum, #reads, #high-quality reads.
     struct LStat { int64_t gap_sum = 0, gap_sq = 0, count = 0, hiq = 0; };
     std::map< std::pair<int,int>, LStat > link_map;
     std::map< int, std::map<int,int> > trans_votes;
     std::map< std::pair<int,int>, LStat > jump_gap;
     std::map< std::pair<int,int>, std::map<int,int> > pair_votes;
     // P4: per-junction read gap samples (for ONT-consensus
     // filling): junction -> (read name, [r1,r2) of the gap
     // content on that read).  The samples are flank-anchored:
     // each starts exactly at the left edge's read-alignment end,
     // so a plain column majority over samples is the consensus.
     std::map< std::pair<int,int>,
          vec< std::pair< String, std::pair<int,int> > > > gap_reads;
     struct Trans { int e1, e2; Bool adjacent; int gap; };
     auto vote_transitions = [&]( vec<Trans> const& tr )
     {    int n_links_read = 0;
          vec< std::pair<int,int> > tseq;
          for ( auto const& t : tr )
          {    int e1 = t.e1, e2 = t.e2;
               trans_votes[e1][e2]++;
               trans_votes[ inv[e2] ][ inv[e1] ]++;
               tseq.push_back( {e1, e2} );
               if ( t.adjacent ) continue;   // standard overlap
               auto jkeys = { std::make_pair(e1,e2),
                    std::make_pair( int(inv[e2]), int(inv[e1]) ) };
               for ( auto const& key : jkeys )
               {    auto& js = jump_gap[key];
                    js.gap_sum += t.gap;
                    js.gap_sq += int64_t(t.gap) * t.gap;
                    js.count += 1;
                    js.hiq += 1;    }
               if ( edge_class[e1] == EDGE_NORMAL
                    && edge_class[e2] == EDGE_NORMAL )
               {    auto& entry = link_map[ {e1,e2} ];
                    entry.gap_sum += t.gap;
                    entry.gap_sq += int64_t(t.gap) * t.gap;
                    entry.count += 1;
                    entry.hiq += 1;
                    n_links_read++;    }    }
          for ( int u = 0; u + 1 < tseq.isize( ); u++ )
          {    if ( tseq[u].second != tseq[u+1].first ) continue;
               int R = tseq[u].second;
               int a2 = tseq[u].first, b2 = tseq[u+1].second;
               if ( a2 == R || b2 == R ) continue;
               pair_votes[ {R, a2} ][b2]++;
               pair_votes[ { int(inv[R]), int(inv[b2]) } ]
                    [ int(inv[a2]) ]++;    }
          {    int u = 0;
               while ( u < tseq.isize( ) )
               {    if ( edges[ tseq[u].second ].isize( ) >= 2000 )
                    {    u++;    continue;    }
                    int a3 = tseq[u].first;
                    int v = u;
                    while ( v < tseq.isize( )
                         && edges[ tseq[v].second ].isize( ) < 2000 )
                         v++;
                    if ( v < tseq.isize( ) )
                    {    int b3 = tseq[v].second;
                         for ( int i = u; i < v; i++ )
                         {    int x = tseq[i].second;
                              pair_votes[ {x, a3} ][b3]++;
                              pair_votes[ { int(inv[x]),
                                   int(inv[b3]) } ]
                                   [ int(inv[a3]) ]++;    }    }
                    u = v;    }    }
          return n_links_read;    };
     String line, header, raw, seq;
     Ofstream( rstat, out_head + ".readstats.txt" );
     rstat << "# id len n_anchors n_lookups hit_rate span"
           << " chain_cov n_blocks n_links n_norm_edges\n";

     if ( mapper == "mm2" )
     {    // ---- 3a. minimap2 mapping: one-strand reference, PAF in ----

          // One strand = one orientation of every genomic path.
          // Build maximal lines over the FULL doubled graph,
          // resolving branches toward the longest successor (the
          // backbone), so lines run through variant bubbles and
          // include their short interior edges.  Then keep one
          // orientation of each line/mirror pair: the reference
          // carries each genomic path once, in one strand.
          vec< vec<int> > lines;
          vec<Bool> used( hb.E( ), False );
          auto build_line = [&]( int e0 )
          {    vec<int> ln;
               int cur = e0;
               while ( cur >= 0 && !used[cur] )
               {    used[cur] = True;
                    ln.push_back(cur);
                    int w = to_right[cur];
                    if ( hb.FromSize(w) == 0 ) break;
                    int nxt = -1, nlen = -1;
                    for ( int j = 0; j < hb.FromSize(w); j++ )
                    {    int e2 = hb.IFrom( w, j );
                         if ( used[e2] ) continue;
                         int l2 = edges[e2].isize( );
                         if ( l2 > nlen
                              || ( l2 == nlen && e2 < nxt ) )
                         {    nlen = l2;    nxt = e2;    }    }
                    if ( nxt < 0 ) break;
                    cur = nxt;    }
               if ( !ln.empty( ) ) lines.push_back(ln);    };
          for ( int e0 = 0; e0 < hb.E( ); e0++ )
               if ( !used[e0] ) build_line(e0);
          {    int nall = lines.isize( );
               vec<Bool> keep( lines.isize( ), True );
               for ( int i = 0; i < lines.isize( ); i++ )
               {    vec<int> const& L = lines[i];
                    int k = L.isize( );
                    for ( int t = 0; t < k; t++ )
                    {    int a = L[t], b = inv[ L[k-1-t] ];
                         if ( a != b )
                         {    keep[i] = ( a < b );
                              break;    }    }    }
               vec< vec<int> > kept;
               for ( int i = 0; i < lines.isize( ); i++ )
                    if ( keep[i] ) kept.push_back( lines[i] );
               lines = kept;
               cout << Date( ) << ": mm2: one strand: " << nall
                    << " lines -> " << lines.isize( )
                    << " after inv-dedup" << endl;    }

          // Offset tables: edge j of line i occupies line-offsets
          // [lcum, lcum+len); the mirror line (inv-reversed edges)
          // has its own table for '-' mappings.
          vec< vec<int64_t> > lcum( lines.isize( ) );
          vec< vec<int> > mline( lines.isize( ) );
          vec< vec<int64_t> > mcum( lines.isize( ) );
          vec<int64_t> lseq_len( lines.isize( ), 0 );
          for ( int i = 0; i < lines.isize( ); i++ )
          {    vec<int> const& le = lines[i];
               int64_t c = 0;
               for ( int j = 0; j < le.isize( ); j++ )
               {    lcum[i].push_back(c);
                    c += edges[ le[j] ].isize( ) - OVER;    }
               int k = le.isize( );
               lseq_len[i] = lcum[i][k-1]
                    + edges[ le[k-1] ].isize( );
               for ( int t = 0; t < k; t++ )
               {    mline[i].push_back( inv[ le[k-1-t] ] );
                    int64_t lo = lcum[i][k-1-t]
                         + edges[ le[k-1-t] ].isize( );
                    mcum[i].push_back( lseq_len[i] - lo );    }    }
          {    int64_t rtot = 0;
               vec<int64_t> rlens;
               for ( int i = 0; i < lines.isize( ); i++ )
               {    rlens.push_back( lseq_len[i] );
                    rtot += lseq_len[i];    }
               Sort(rlens);
               int64_t rn50 = 0;
               {    int64_t acc = 0;
                    for ( int i = rlens.isize( ) - 1; i >= 0; i-- )
                    {    acc += rlens[i];
                         if ( acc >= rtot / 2 )
                         {    rn50 = rlens[i];
                              break;    }    }    }
               cout << Date( ) << ": mm2: " << lines.isize( )
                    << " one-strand lines (N50 "
                    << ToStringAddCommas(rn50) << " bp, total "
                    << ToStringAddCommas(rtot) << " bp)" << endl;    }

          String ref = out_head + ".mm2ref.fa";
          {    Ofstream( rout, ref );
               for ( int i = 0; i < lines.isize( ); i++ )
               {    vec<int> const& le = lines[i];
                    rout << ">p_" << i << "\n";
                    int64_t written = 0;
                    for ( int j = 0; j < le.isize( ); j++ )
                    {    bvec const& bv = edges[ le[j] ];
                         int from = ( j == 0 ? 0 : OVER );
                         for ( int b = from; b < bv.isize( ); b++ )
                         {    rout << Base::val2Char( bv[b] );
                              if ( ++written % 70 == 0 )
                                   rout << "\n";    }    }
                    if ( written % 70 != 0 ) rout << "\n";    }    }

          if ( !IsRegularFile( mm2_bin ) )
               FatalErr( "minimap2 binary not found: " << mm2_bin
                    << "\nbuild it with:  make -C minimap2/minimap2"
                    << " aarch64=1   (x86: drop aarch64=1)\n"
                    << "or use MAPPER=seed for the internal aligner." );
          String mcmd = mm2_bin
               + " -x map-ont -N 0 --paf-no-hit -c -t 4 "
               + ref + " " + reads_file;
          cout << Date( ) << ": mm2: " << mcmd << endl;
          fast_pipe_ifstream in(mcmd);

          // PAF columns: qname qlen qstart qend strand tname tlen
          // tstart tend mlen blen mapq [tags].
          struct MSeg { int e; int o1, o2; int r1, r2; };
          int dbg_printed = 0;
          int64_t n_tr_same = 0, n_tr_rep = 0, n_tr_deep = 0,
               n_tr_ovl = 0, n_tr_oob = 0, n_tr_slop = 0,
               n_tr_ok = 0;
          auto flush_read = [&]( String const& qname, int qlen,
               vec<MSeg>& segs )
          {    if ( qname.empty( ) ) return;
               int64_t id = n_reads++;
               total_bp += qlen;
               if ( qlen < min_read || segs.isize( ) == 0 ) return;
               std::sort( segs.begin( ), segs.end( ),
                    []( MSeg const& a, MSeg const& b )
                    {    return a.r1 < b.r1;    } );
               // Coarse-path mode: keep only long per-edge segments
               // (MIN_SEG=1000 = the classic map-to-contigs rule:
               // only >1kb alignments vote; the small variant and
               // hub edges ride along inside long line mappings but
               // do not become waypoints).
               if ( min_seg > 0 )
               {    vec<MSeg> kept;
                    for ( auto const& s : segs )
                         if ( s.o2 - s.o1 >= min_seg )
                              kept.push_back(s);
                    segs = kept;    }
               if ( segs.isize( ) == 0 ) return;
               n_ge++;
               // Mapping quality: aligned fraction of the read.
               int mapped_bp = 0;
               for ( auto const& s : segs )
                    mapped_bp += s.r2 - s.r1;
               double chain_cov =
                    double(mapped_bp) / double( Max( 1, qlen ) );
               Bool hiq_read = ( chain_cov >= min_chain_cov );
               int n_norm_edges = 0;
               {    std::set<int> ne;
                    for ( auto const& s : segs )
                         if ( edge_class[ s.e ] == EDGE_NORMAL )
                              ne.insert( s.e );
                    n_norm_edges = ne.size( );    }
               int n_links_read = 0;
               Bool dbg = ( VERBOSE && !debug_read.empty( )
                    && qname.Contains( debug_read ) );
               if ( dbg )
               {    cout << "  [dbg] read " << qname << " len="
                         << qlen << " segs=" << segs.size( )
                         << " chain_cov=" << fixed
                         << setprecision(3) << chain_cov
                         << " hiq=" << ( hiq_read ? "yes" : "no" )
                         << resetiosflags(ios::fixed) << endl;
                    int j2 = 0;
                    for ( auto const& s : segs )
                         cout << "  [dbg]   seg[" << j2++
                              << "] edge " << s.e << " (len "
                              << edges[ s.e ].isize( ) << ") off "
                              << s.o1 << "-" << s.o2 << "  read "
                              << s.r1 << "-" << s.r2 << endl;    }
               if ( hiq_read )
               {    n_hiq_reads++;
                    if ( VERBOSE && segs.isize( ) >= 3
                         && dbg_printed < 5 )
                    {    dbg_printed++;
                         cout << "  mm2 read " << qname << " ("
                              << segs.size( ) << " segs):";
                         for ( auto const& s : segs )
                              cout << " e" << s.e << "[" << s.o1
                                   << "-" << s.o2 << "@"
                                   << s.r1 << "-" << s.r2 << "]";
                         cout << endl;    }
                    vec<Trans> tr;
                    int c_same = 0, c_rep = 0, c_deep = 0,
                         c_ovl = 0, c_oob = 0, c_slop = 0,
                         c_ok = 0;
                    for ( int i = 0; i + 1 < segs.isize( ); i++ )
                    {    MSeg const& a = segs[i];
                         MSeg const& b = segs[i+1];
                         if ( a.e == b.e ) {    c_same++;    continue;    }
                         // Repeats are never trusted as waypoints.
                         if ( edge_class[a.e] == EDGE_REPEAT
                              || edge_class[b.e] == EDGE_REPEAT )
                         {    c_rep++;    continue;    }
                         int dist = b.r1 - a.r2;
                         int e1 = a.e, e2 = b.e;
                         Bool adjacent =
                              ( to_right[e1] == to_left[e2] );
                         if ( dbg )
                              cout << "  [dbg]   pair e" << e1
                                   << "@r" << a.r1 << "-" << a.r2
                                   << " -> e" << e2 << "@r"
                                   << b.r1 << "-" << b.r2
                                   << "  dist=" << dist
                                   << ( adjacent ? "  ADJACENT"
                                        : "  NOT-ADJACENT" )
                                   << endl;
                         // Adjacent edges share K-1 bases, so their
                         // mappings overlap on the read by ~OVER; a
                         // deeper overlap (or any overlap between
                         // non-adjacent edges) is not a thread.
                         if ( dist < -OVER - 50 )
                         {    if ( dbg )
                                   cout << "  [dbg]     REJECT"
                                        << " (deep overlap)"
                                        << endl;
                              c_deep++;    continue;    }
                         if ( !adjacent && dist <= 0 )
                         {    if ( dbg )
                                   cout << "  [dbg]     REJECT"
                                        << " (overlap on"
                                        << " non-adjacent edges)"
                                        << endl;
                              c_ovl++;    continue;    }
                         if ( adjacent )
                         {    // Junction algebra sanity (linear,
                              // indel-free approximation; the seg
                              // boundaries themselves are
                              // cigar-exact).  The K-1 overlap must
                              // roughly match the read distance.
                              int len_a = edges[e1].isize( );
                              int expected = a.o2 + dist - len_a
                                   + OVER;
                              if ( abs( b.o1 - expected ) > Slop(dist) )
                              {    if ( dbg )
                                        cout << "  [dbg]     REJECT"
                                             << " (junction offset"
                                             << " off by "
                                             << abs( b.o1 - expected )
                                             << " > "
                                             << Slop(dist) << ")"
                                             << endl;
                                   c_slop++;    continue;    }    }
                         if ( dbg )
                              cout << "  [dbg]     ACCEPT ("
                                   << ( adjacent ? "adjacency vote"
                                        : "jump vote, gap=" )
                                   << ( adjacent ? "" :
                                        ToString(dist) )
                                   << ")" << endl;
                         c_ok++;
                         tr.push_back( { e1, e2, adjacent, dist } );
                         // P4: remember the read span of the gap
                         // content for ONT-consensus filling.
                         if ( !adjacent )
                              gap_reads[ {e1,e2} ].push_back(
                                   { qname, { a.r2, b.r1 } } );    }
                    n_tr_same += c_same;    n_tr_rep += c_rep;
                    n_tr_deep += c_deep;    n_tr_ovl += c_ovl;
                    n_tr_oob += c_oob;    n_tr_slop += c_slop;
                    n_tr_ok += c_ok;
                    n_links_read = vote_transitions(tr);
                    if ( dbg )
                    {    for ( auto const& t : tr )
                         {    cout << "  [dbg]   stored e" << t.e1
                                   << " -> e" << t.e2
                                   << ( t.adjacent
                                        ? " (adjacency)" : " (jump)" )
                                   << endl;
                              cout << "  [dbg]     trans_votes["
                                   << t.e1 << "][" << t.e2 << "] = "
                                   << trans_votes[ t.e1 ][ t.e2 ]
                                   << "   mirror trans_votes["
                                   << inv[t.e2] << "][" << inv[t.e1]
                                   << "] = " << trans_votes[
                                        inv[t.e2] ][ inv[t.e1] ]
                                   << endl;
                              if ( !t.adjacent )
                              {    auto jit = jump_gap.find(
                                        { t.e1, t.e2 } );
                                   if ( jit != jump_gap.end( ) )
                                        cout << "  [dbg]     "
                                             << "jump_gap[e" << t.e1
                                             << ",e" << t.e2
                                             << "] count="
                                             << jit->second.count
                                             << " mean_gap="
                                             << ( jit->second.count
                                                  ? double(
                                                  jit->second.gap_sum)
                                                  / jit->second.count
                                                  : 0.0 ) << endl;
                                   if ( edge_class[t.e1]
                                             == EDGE_NORMAL
                                        && edge_class[t.e2]
                                             == EDGE_NORMAL )
                                        cout << "  [dbg]     "
                                             << "link_map[e" << t.e1
                                             << ",e" << t.e2
                                             << "] count="
                                             << link_map[
                                                  { t.e1, t.e2 } ]
                                                  .count << endl;    }    }
                         // Entry-exit pairing votes from this read's
                         // transition triples (a -> R -> b).
                         for ( int u = 0; u + 1 < tr.isize( ); u++ )
                              if ( tr[u].e2 == tr[u+1].e1 )
                              {    int R = tr[u].e2;
                                   int a2 = tr[u].e1,
                                        b2 = tr[u+1].e2;
                                   cout << "  [dbg]   pair_votes["
                                        << "e" << R << ", entry e"
                                        << a2 << "] -> exit e" << b2
                                        << " = " << pair_votes[
                                             { R, a2 } ][ b2 ]
                                        << "   mirror [e"
                                        << inv[R] << ", entry e"
                                        << inv[b2] << "] -> exit e"
                                        << inv[a2] << " = "
                                        << pair_votes[
                                             { int(inv[R]),
                                               int(inv[b2]) } ]
                                             [ int(inv[a2]) ]
                                        << endl;    }    }    }
               int span = segs.back( ).r2 - segs.front( ).r1;
               rstat << id << " " << qlen << " " << segs.size( )
                    << " 0 " << fixed << setprecision(3)
                    << chain_cov << " " << span << " " << chain_cov
                    << " " << segs.size( ) << " " << n_links_read
                    << " " << n_norm_edges << "\n"
                    << resetiosflags(ios::fixed);    };

          String pline, cur_name;
          int cur_qlen = 0;
          vec<MSeg> cur_segs;
          while ( true )
          {    getline( in, pline );
               if ( in.fail( ) ) break;
               if ( pline.empty( ) ) continue;
               // Tab-split the PAF line.
               vec<String> f;
               {    size_t p = 0;
                    while ( true )
                    {    size_t q = pline.find( '\t', p );
                         if ( q == String::npos )
                         {    f.push_back( pline.substr( p ) );
                              break;    }
                         f.push_back( pline.substr( p, q - p ) );
                         p = q + 1;    }    }
               if ( f.size( ) < 12 ) continue;
               String qname = f[0];
               if ( qname != cur_name )
               {    flush_read( cur_name, cur_qlen, cur_segs );
                    cur_name = qname;
                    cur_segs.clear( );    }
               cur_qlen = f[1].Int( );
               if ( f[5] == "*" ) continue;   // unmapped read
               if ( VERBOSE && !debug_read.empty( )
                    && f[0].Contains( debug_read ) )
                    cout << "  [dbg] paf: " << pline << endl;
               int qs = f[2].Int( ), qe = f[3].Int( );
               String strand = f[4];
               if ( f[5].size( ) < 3 ) continue;
               int p = f[5].substr( 2 ).Int( );   // "p_<id>"
               int ts = f[7].Int( ), te = f[8].Int( );
               int mlen = f[9].Int( );    // residue matches
               int blen = f[10].Int( );   // alignment block length
               // (f[11] is MAPQ)
               if ( p < 0 || p >= lines.isize( ) ) continue;
               // Mapping-quality gates.
               if ( qe - qs < min_map_len ) continue;
               int nm = -1;
               for ( int u = 12; u < f.isize( ); u++ )
                    if ( f[u].size( ) >= 5
                         && f[u].substr( 0, 5 ) == "NM:i:" )
                    {    nm = f[u].substr( 5 ).Int( );
                         break;    }
               double ident = ( nm >= 0 && blen > 0 )
                    ? 1.0 - double(nm) / blen
                    : ( blen > 0 ? double(mlen) / blen : 0.0 );
               if ( ident < min_map_id ) continue;
               // Cigar: exact read<->target offsets for the interior
               // line-edge boundaries (linear interpolation would
               // drift with ONT indels).
               String cigs;
               for ( int u = 12; u < f.isize( ); u++ )
                    if ( f[u].size( ) >= 5
                         && f[u].substr( 0, 5 ) == "cg:Z:" )
                    {    cigs = f[u].substr( 5 );
                         break;    }
               // M-op table: (target start, query-z start, length).
               vec< std::pair< std::pair<int64_t,int64_t>,
                    int64_t> > mops;
               if ( cigs.size( ) >= 3 )
               {    int64_t T = ts, Z = 0;
                    size_t ci = 0;
                    while ( ci < cigs.size( ) )
                    {    int64_t n = 0;
                         Bool got = False;
                         while ( ci < cigs.size( )
                              && cigs[ci] >= '0'
                              && cigs[ci] <= '9' )
                         {    n = n * 10 + ( cigs[ci] - '0' );
                              ci++;
                              got = True;    }
                         if ( !got || ci >= cigs.size( ) ) break;
                         char cop = cigs[ci++];
                         if ( cop == 'M' || cop == '='
                              || cop == 'X' )
                         {    mops.push_back(
                                   { { T, Z }, n } );
                              T += n;
                              Z += n;    }
                         else if ( cop == 'I' ) Z += n;
                         else if ( cop == 'D' ) T += n;    }    }
               // Read position of a line offset (in the cigar's
               // query-z space: '+' z counts from qs, '-' z counts
               // from the reverse-complemented read segment).
               auto read_at = [&]( int64_t Tb ) -> int
               {    if ( mops.empty( ) )
                    {    if ( strand == "+" )
                              return qs + int(Tb - ts);
                         double sc = double(qe - qs)
                              / double( Max( int64_t(1),
                                   int64_t(te) - ts ) );
                         return qs + int( sc * ( te - 1 - Tb ) );    }
                    int lo = 0, hi = mops.isize( ) - 1, best = -1;
                    while ( lo <= hi )
                    {    int mid = ( lo + hi ) / 2;
                         if ( mops[mid].first.first <= Tb )
                         {    best = mid;
                              lo = mid + 1;    }
                         else hi = mid - 1;    }
                    if ( best < 0 ) return qs;
                    int64_t z = mops[best].first.second
                         + Min( Tb - mops[best].first.first,
                              mops[best].second );
                    return ( strand == "+" )
                         ? qs + int(z) : qe - 1 - int(z);    };
               // Split the line mapping into per-edge segments:
               // '+' walks the line forward, '-' walks the mirror
               // line (inv-reversed edges, flipped offsets).
               vec<int> const& le = lines[p];
               vec<int64_t> const& lc = lcum[p];
               int64_t S = lseq_len[p];
               if ( strand == "+" )
               {    for ( int j = 0; j < le.isize( ); j++ )
                    {    int64_t st = Max( int64_t(ts), lc[j] );
                         int64_t en = Min( int64_t(te),
                              lc[j] + edges[ le[j] ].isize( ) );
                         if ( en <= st ) continue;
                         int r1 = read_at(st);
                         int r2 = read_at(en);
                         if ( r2 <= r1 ) continue;
                         cur_segs.push_back( { le[j],
                              int( st - lc[j] ), int( en - lc[j] ),
                              r1, r2 } );    }    }
               else
               {    vec<int> const& me = mline[p];
                    vec<int64_t> const& mc = mcum[p];
                    int64_t o1t = S - te, o2t = S - ts;
                    for ( int j = 0; j < me.isize( ); j++ )
                    {    int64_t st = Max( o1t, mc[j] );
                         int64_t en = Min( o2t,
                              mc[j] + edges[ me[j] ].isize( ) );
                         if ( en <= st ) continue;
                         // Mirror offset b: inclusive start at
                         // line offset S-1-st; exclusive end at
                         // line offset S-en (read pos qe-z(S-en)).
                         int r1 = read_at( S - 1 - st );
                         int r2 = read_at( S - en ) + 1;
                         if ( r2 <= r1 ) continue;
                         cur_segs.push_back( { me[j],
                              int( st - mc[j] ), int( en - mc[j] ),
                              r1, r2 } );    }    }    }
          flush_read( cur_name, cur_qlen, cur_segs );
          cout << Date( ) << ": mm2 transitions: ok=" << n_tr_ok
               << " same=" << n_tr_same << " repeat=" << n_tr_rep
               << " deep-overlap=" << n_tr_deep
               << " nonadj-overlap=" << n_tr_ovl
               << " junction-oob=" << n_tr_oob
               << " junction-slop=" << n_tr_slop << endl;
     }
     else
     {
     // ---- 2. Build k-mer index ----

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
     cout << Date( ) << ": index built, "
          << ToStringAddCommas(index.size( )) << " entries" << endl;

     // ---- 3. Stream reads; seed; block-chain; collect raw links ----

     String cmd = "cat " + reads_file;
     if ( reads_file.Contains( ".gz", -1 ) ) cmd = "zcat " + reads_file;
     fast_pipe_ifstream in(cmd);

     auto sameOrAdjacent = [&]( Anchor const& a, Anchor const& b )
     {    int dist = int(b.r) - int(a.r);
          if ( dist <= 0 ) return false;
          if ( a.e == b.e )
          {    int diff = int(b.o) - int(a.o);
               return diff > 0 && abs( diff - dist ) <= Slop(dist);    }
          if ( to_right[ a.e ] != to_left[ b.e ] ) return false;
          int len_a = edges[ a.e ].isize( );
          int expected = int(a.o) + dist - len_a + OVER;
          if ( expected < 0 ) return false;
          if ( expected + int(K) > edges[ b.e ].isize( ) ) return false;
          return abs( int(b.o) - expected ) <= Slop(dist);    };

     auto process = [&]( String const& header, String const& raw )
     {    int64_t id = n_reads++;
          total_bp += raw.size( );
          if ( (int64_t) raw.size( ) < min_read ) return;
          n_ge++;

          seq.resize(0);
          seq.reserve(raw.size( ));
          for ( int i = 0; i < raw.isize( ); i++ )
               seq.push_back( ValidBase(raw[i]) ? raw[i] : 'A' );
          bvec rd(seq);

          // Seed with exact k-mer matches; k-mers occurring more than
          // max_mult times in the graph are ambiguous (repeat copies /
          // diploid variant pairs) and are not used as anchors.
          vec<Anchor> anchors;
          int next_anchor_pos = 0;
          int n_lookups = 0, n_hits = 0;
          uint64_t const mask = ( uint64_t(1) << (2*K) ) - 1;
          uint64_t code = 0;
          for ( size_t j = 0; j < rd.size( ); j++ )
          {    code = ( ( code << 2 ) | rd[j] ) & mask;
               if ( j+1 < K ) continue;
               int pos = int(j) - int(K) + 1;
               if ( pos < next_anchor_pos ) continue;
               n_lookups++;
               auto lo = std::lower_bound( index.begin(), index.end(),
                    KEntry{ code, 0, 0 },
                    []( KEntry const& a, KEntry const& b )
                    {    return a.code < b.code;    } );
               auto hi = std::upper_bound( index.begin(), index.end(),
                    KEntry{ code, 0, 0 },
                    []( KEntry const& a, KEntry const& b )
                    {    return a.code < b.code;    } );
               if ( lo != hi ) n_hits++;
               if ( lo == hi ) continue;
               if ( hi - lo > max_mult ) continue;
               for ( auto it = lo; it != hi; ++it )
                    anchors.push_back( { uint32_t(pos), it->edge,
                         it->pos } );
               next_anchor_pos = pos + anchor_step;    }

          // Mapping quality of this read.
          //  hit_rate  = fraction of sampled k-mers that matched the
          //              graph (identity proxy; ONT errors leave ~70%
          //              of 19-mers intact on well-mapped sequence).
          //  chain_cov = anchored span / read length (coverage proxy).
          double hit_rate = n_lookups ? double(n_hits) / n_lookups : 0.0;
          double chain_cov = 0.0;
          if ( anchors.size( ) >= 2 )
          {    int span = int(anchors.back().r) + int(K)
                    - int(anchors.front().r);
               chain_cov = double(span) / double(raw.size( ));    }
          Bool hiq_read = ( chain_cov >= min_chain_cov );

          if ( anchors.size( ) < 2 ) return;

          // Greedy block-chain.
          vec< std::pair<int,int> > blocks;
          int block_start = 0;
          for ( int i = 1; i < anchors.isize( ); i++ )
          {    if ( !sameOrAdjacent( anchors[i-1], anchors[i] ) )
               {    blocks.push_back( { block_start, i-1 } );
                    block_start = i;    }    }
          blocks.push_back( { block_start, anchors.isize( ) - 1 } );

          // Extract scaffold links, filtering repeats and short edges.
          // A link is only trusted if both flanking blocks are well
          // supported (>= min_anchor_block anchors each): a break
          // between single-anchor blocks is usually a spurious seed,
          // not a real graph gap.
          int n_links_read = 0, n_norm_edges = 0;
          {    std::set<int> ne;
               for ( auto const& a : anchors )
                    if ( edge_class[ a.e ] == EDGE_NORMAL )
                         ne.insert( int(a.e) );
               n_norm_edges = ne.size( );    }
          if ( hiq_read )
          {    // Vote transitions from consecutive ANCHOR pairs.
               // Chained pairs that cross an edge boundary vote the
               // graph adjacency (this is the thread through the
               // graph: variant bubbles, short edges included).
               // Chain breaks vote a jump (real gap) when both
               // flanking blocks are well supported.
               vec<int> blk( anchors.isize( ), 0 );
               for ( int t = 0; t < blocks.isize( ); t++ )
                    for ( int i = blocks[t].first;
                         i <= blocks[t].second; i++ )
                         blk[i] = t;
               vec<Trans> tr;
               if ( min_seg > 0 )
               {    // Coarse mode, the same rule as the mm2 MIN_SEG
                    // filter: only long edges are waypoints.  Small
                    // variant and hub edges cannot become waypoints
                    // or be over-crossed; their anchors are bridged
                    // over, so the vote is waypoint-to-waypoint with
                    // the read distance as the gap.  Both flanking
                    // waypoint occurrences must be supported by at
                    // least min_anchor_block anchors -- the analog
                    // of the mm2 mapping-length gate.
                    vec<int> wo_e, wo_first, wo_last, wo_cnt;
                    for ( int i = 0; i < anchors.isize( ); i++ )
                    {    Anchor const& a = anchors[i];
                         if ( edges[ a.e ].isize( ) < min_seg )
                              continue;
                         if ( !wo_e.empty( )
                              && wo_e.back( ) == int(a.e) )
                         {    wo_last.back( ) = int(a.r);
                              wo_cnt.back( )++;    }
                         else
                         {    wo_e.push_back( int(a.e) );
                              wo_first.push_back( int(a.r) );
                              wo_last.push_back( int(a.r) );
                              wo_cnt.push_back(1);    }    }
                    for ( int i = 0; i + 1 < wo_e.isize( ); i++ )
                    {    if ( wo_cnt[i] < min_anchor_block
                              || wo_cnt[i+1] < min_anchor_block )
                              continue;
                         int e1 = wo_e[i], e2 = wo_e[i+1];
                         if ( edge_class[e1] == EDGE_REPEAT
                              || edge_class[e2] == EDGE_REPEAT )
                              continue;
                         Bool adjacent =
                              ( to_right[e1] == to_left[e2] );
                         tr.push_back( { e1, e2, adjacent,
                              wo_first[i+1] - wo_last[i] } );    }    }
               else
               for ( int i = 1; i < anchors.isize( ); i++ )
               {    Anchor const& a = anchors[i-1];
                    Anchor const& b = anchors[i];
                    if ( a.e == b.e ) continue;
                    // Repeats are never trusted as waypoints.
                    if ( edge_class[ a.e ] == EDGE_REPEAT ) continue;
                    if ( edge_class[ b.e ] == EDGE_REPEAT ) continue;
                    int e1 = int(a.e), e2 = int(b.e);
                    Bool adjacent = ( to_right[e1] == to_left[e2] );
                    Bool chained = sameOrAdjacent( a, b );
                    if ( chained && adjacent )
                         tr.push_back( { e1, e2, True, 0 } );
                    else
                    {    // Break: require support on both flanks --
                         // single-anchor blocks are usually spurious
                         // seeds, not real transitions.
                         int na = blocks[ blk[i-1] ].second
                              - blocks[ blk[i-1] ].first + 1;
                         int nb = blocks[ blk[i] ].second
                              - blocks[ blk[i] ].first + 1;
                         if ( na < min_anchor_block
                              || nb < min_anchor_block ) continue;
                         tr.push_back( { e1, e2, adjacent,
                              int(b.r) - int(a.r) } );    }    }
               n_links_read += vote_transitions(tr);    }
          if ( hiq_read ) n_hiq_reads++;

          rstat << id << " " << raw.size( ) << " " << anchors.size( )
                << " " << n_lookups << " " << fixed << setprecision(3)
                << hit_rate << " "
                << int(anchors.back().r) + int(K) - int(anchors.front().r)
                << " " << chain_cov << " " << blocks.size( ) << " "
                << n_links_read << " " << n_norm_edges << "\n"
                << resetiosflags(ios::fixed);    };

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
               getline( in, line );
               getline( in, line );
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
     }

     cout << Date( ) << ": " << n_reads << " reads (" << n_ge
          << " >= " << min_read << ", " << n_hiq_reads
          << " hi-q), " << link_map.size( )
          << " normal-jump links, " << trans_votes.size( )
          << " edges with thread votes" << endl;
     rstat.close( );

     // ---- 4. Link quality diagnostics ----

     // Normal-to-normal jumps (true scaffold gaps) with per-link gap
     // consistency; written for eyeballing and for P4 (gap filling).

     {    Ofstream( lqout, out_head + ".linkstats.txt" );
          lqout << "# left right n_reads n_hiq gap_mean gap_sd\n";
          for ( auto const& kv : link_map )
          {    int left = kv.first.first, right = kv.first.second;
               LStat const& st = kv.second;
               if ( st.count == 0 ) continue;
               double mean = double(st.gap_sum) / st.count;
               double var = double(st.gap_sq) / st.count - mean*mean;
               double sd = var > 0 ? std::sqrt(var) : 0.0;
               lqout << left << " " << right << " " << st.count
                    << " " << st.hiq << " " << fixed
                    << setprecision(1) << mean << " " << sd
                    << "\n" << resetiosflags(ios::fixed);    }    }

     // ---- 5. Consensus successor map (read threading) ----

     // Every consecutive pair of well-supported blocks votes the
     // transition e1 -> e2 (in both mirror worlds, so reads of both
     // orientations agree).  The successor of an edge is its
     // most-voted transition; branches (variant bubbles, repeat
     // exits) are resolved by read consensus, which is diploid-aware:
     // both haplotypes vote their own paths, and the majority path
     // forms the primary scaffold while minority paths surface as
     // separate short paths.

     std::map<int,int> succ, succ_sup;
     for ( auto const& kv : trans_votes )
     {    int best = -1, bestv = 0;
          for ( auto const& s : kv.second )
               if ( s.second > bestv ) { bestv = s.second; best = s.first; }
          if ( best >= 0 && bestv >= min_support )
          {    succ[kv.first] = best;
               succ_sup[kv.first] = bestv;    }    }
     std::set<int> has_pred;
     for ( auto const& kv : succ ) has_pred.insert(kv.second);
     int n_adj_succ = 0, n_jump_succ = 0;
     for ( auto const& kv : succ )
          if ( to_right[kv.first] == to_left[kv.second] ) n_adj_succ++;
          else n_jump_succ++;
     cout << Date( ) << ": threading: " << trans_votes.size( )
          << " voted edges, " << succ.size( ) << " successors ("
          << n_adj_succ << " adjacent, " << n_jump_succ
          << " jump), min votes " << min_support << endl;

     // ---- 6. Walk scaffold paths ----

     // Start at canonical normal edges without a predecessor and
     // follow the consensus successors through the whole graph
     // (short edges, variant bubbles included).  A walk that returns
     // to its start closes a circle.

     vec<int> scaffold_nodes;
     for ( int e = 0; e < hb.E( ); e++ )
     {    if ( e > inv[e] ) continue;
          if ( edge_class[e] != EDGE_NORMAL ) continue;
          scaffold_nodes.push_back(e);    }
     cout << Date( ) << ": " << scaffold_nodes.size( )
          << " canonical scaffold nodes" << endl;

     // Order starts by descending edge length: the first (longest)
     // walk claims the main genomic path before any other.
     std::sort( scaffold_nodes.begin( ), scaffold_nodes.end( ),
          [&]( int a, int b )
          {    return edges[a].isize( ) > edges[b].isize( );    } );

     // Vote-based predecessor map (for backward extension).
     std::map<int,int> pred, pred_sup;
     for ( auto const& kv : succ )
     {    int e2 = kv.second;
          int sup = succ_sup[kv.first];
          auto it = pred.find(e2);
          if ( it == pred.end( ) || sup > pred_sup[e2] )
          {    pred[e2] = kv.first;
               pred_sup[e2] = sup;    }    }

     // Bubble skip: at a branch with no read evidence, if all
     // branches quickly reconverge to one vertex (a variant bubble),
     // continue through the most-voted branch.
     auto bubble_skip = [&]( int cur, int& nxt )
     {    int w = to_right[cur];
          if ( hb.FromSize(w) < 2 ) return False;
          vec< vec<int> > branches;
          vec<int> endv;
          for ( int j = 0; j < hb.FromSize(w); j++ )
          {    vec<int> p;
               int e = hb.IFrom( w, j );
               int64_t size = 0;
               Bool ok = True;
               while ( True )
               {    p.push_back(e);
                    size += edges[e].isize( );
                    if ( p.isize( ) > 10 || size > 20000 )
                    {    ok = False;    break;    }
                    int v = to_right[e];
                    if ( hb.FromSize(v) == 0 )
                    {    ok = False;    break;    }
                    if ( hb.FromSize(v) >= 2 ) break;
                    e = hb.IFrom( v, 0 );    }
               if ( !ok ) return False;
               branches.push_back(p);
               endv.push_back( to_right[ p.back( ) ] );    }
          for ( int i = 1; i < endv.isize( ); i++ )
               if ( endv[i] != endv[0] ) return False;
          int best = 0; int64_t bestv = -1;
          for ( int i = 0; i < branches.isize( ); i++ )
          {    int64_t votes = 0;
               for ( int e : branches[i] )
               {    auto it = trans_votes.find(e);
                    if ( it != trans_votes.end( ) )
                         for ( auto const& s : it->second )
                              votes += s.second;    }
               if ( votes > bestv ) { bestv = votes; best = i; }    }
          nxt = branches[best][0];
          return True;    };

     // Branch guess: at a fork with no read evidence at all (an
     // ONT coverage gap), score each branch by the read support
     // reachable downstream and take the best.  The sequence stays
     // complete and non-duplicated; only the local order can be
     // wrong, which is unavoidable without evidence.
     auto guess_branch = [&]( int cur, int& nxt )
     {    int w = to_right[cur];
          if ( hb.FromSize(w) < 2 ) return False;
          double best_score = 0; int best = -1;
          for ( int j = 0; j < hb.FromSize(w); j++ )
          {    int b = hb.IFrom( w, j );
               double score = 0;
               {    auto tv = trans_votes.find(b);
                    if ( tv != trans_votes.end( ) )
                    {    int64_t votes = 0;
                         for ( auto const& s : tv->second )
                              votes += s.second;
                         if ( double(votes) > score )
                              score = double(votes);    }    }
               std::set<int> seen_v;
               std::deque< std::pair<int,int> > q;   // (vertex, hops)
               seen_v.insert( to_right[b] );
               q.push_back( { to_right[b], 1 } );
               while ( !q.empty( ) )
               {    auto [v, hops] = q.front( );
                    q.pop_front( );
                    if ( hops > 30 ) continue;
                    for ( int u = 0; u < hb.FromSize(v); u++ )
                    {    int e2 = hb.IFrom( v, u );
                         auto tv = trans_votes.find(e2);
                         if ( tv != trans_votes.end( ) )
                         {    int64_t votes = 0;
                              for ( auto const& s : tv->second )
                                   votes += s.second;
                              double sc = double(votes) / ( 1 + hops );
                              if ( sc > score ) score = sc;    }
                         int v2 = to_right[e2];
                         if ( !seen_v.count(v2) && hops < 30 )
                         {    seen_v.insert(v2);
                              q.push_back( { v2, hops+1 } );    }    }    }
               if ( score > best_score )
               {    best_score = score;    best = b;    }    }
          if ( best < 0 ) return False;
          nxt = best;
          return True;    };

     std::set<int> visited;
     vec< vec<int> > paths;
     vec<Bool> circular;
     int64_t n_fallback = 0, n_bfallback = 0, n_repeat_cross = 0,
          n_bubble = 0, n_guess = 0;
     // Global exit accounting with multiplicity: each genomic
     // crossing of a repeat consumes one slot of one of its voted
     // exits; an exit may be reused up to the number of distinct
     // voted predecessors of the edge (the copy-count estimate for
     // collapsed repeats, from read evidence rather than the walk's
     // own behavior).
     std::map<int, std::map<int,int>> exit_use;
     std::map<int,int> entry_count;
     for ( auto const& kv : trans_votes )
          for ( auto const& s : kv.second )
               if ( s.second >= 1 )
                    entry_count[ s.first ]++;

     // Hub diagnostics: multi-exit edges with their pairing evidence.
     if ( VERBOSE )
     {    cout << Date( ) << ": repeat hubs:" << endl;
          for ( auto const& kv : trans_votes )
          {    if ( kv.second.size( ) < 3 ) continue;
               int e = kv.first;
               cout << "  hub " << e << " (len "
                    << edges[e].isize( ) << ", entries="
                    << entry_count[e] << "):";
               for ( auto const& s : kv.second )
                    cout << " " << s.first << ":" << s.second;
               for ( auto const& pk : pair_votes )
                    if ( pk.first.first == e )
                         for ( auto const& s : pk.second )
                              cout << " | p(" << pk.first.second
                                   << ")->" << s.first << ":"
                                   << s.second;
               cout << endl;    }    }
     // Alternate walking and extension: each round starts new walks
     // at still-unvisited scaffold nodes (threading intervals the
     // previous round could not reach) and then merges what they
     // found.  Stops when a round adds neither paths nor joins.
     int64_t round_joins = 0;
     for ( int round = 0; round < 8; round++ )
     {    int paths_before = paths.isize( );
          int64_t joins_before = round_joins;
          for ( int pass = 0; pass < ( round == 0 ? 2 : 1 ); pass++ )
     {    for ( int start : scaffold_nodes )
          {    if ( visited.count(start) ) continue;
               if ( pass == 0 && has_pred.count(start) ) continue;
               vec<int> path;
               Bool circ = False;
               std::map<int,int> appends;
               int cur = start;
               int prev_e = -1;
               int region_entry = -1;   // entry edge of current short run
               Bool region_pass = False;   // inside a legit repeat crossing
               int64_t guard = 0;
               while ( cur >= 0
                    && guard++ < int64_t(4) * hb.E( ) )
               {    // In-path revisit?  (a repeat copy already in
                    // this path.)
                    int in_path = -1;
                    for ( int t = 0; t < path.isize( ); t++ )
                         if ( path[t] == cur )
                         {    in_path = t;    break;    }
                    if ( visited.count(cur) && in_path == 0 )
                    {    // Returned to our own start: full circle.
                         circ = True;
                         break;    }
                    Bool was_visited = ( visited.count(cur) > 0 );
                    int cap = entry_count.count(cur)
                         ? Max( 1, entry_count[cur] ) : 1;
                    if ( was_visited )
                    {    // A repeat copy (this path's or another's).
                         Bool long_edge =
                              edges[cur].isize( ) >= 2000;
                         Bool cross = False;
                         if ( !long_edge && prev_e >= 0 )
                         {    // Spanning-read evidence for this entry:
                              // cross this copy of the repeat.
                              if ( region_entry >= 0
                                   && pair_votes.count(
                                        { cur, region_entry } ) )
                                   cross = True;
                              // Continuation of a legit crossing
                              // through the repeat's shared short
                              // interior.
                              else if ( region_pass )
                              {    auto it = trans_votes.find(prev_e);
                                   if ( it != trans_votes.end( )
                                        && it->second.count(cur) )
                                        cross = True;    }    }
                         if ( !cross )
                         {    // Slot rule: exit reuse is bounded by
                              // the distinct voted predecessors.
                              Bool slots = False;
                              auto it = trans_votes.find(cur);
                              if ( it != trans_votes.end( ) )
                                   for ( auto const& s
                                        : it->second )
                                   {    if ( s.second < min_support )
                                             continue;
                                        if ( exit_use[cur][s.first]
                                             < cap )
                                        {    slots = True;
                                             break;    }    }
                              if ( !slots )
                              {    if ( in_path < 0 ) break;
                                   // Came around to our own path:
                                   // the walk is a closed cycle.
                                   // Rotate so the closure point is
                                   // the first edge.
                                   circ = True;
                                   if ( in_path > 0 )
                                   {    vec<int> rot;
                                        for ( int t = in_path;
                                             t < path.isize( ); t++ )
                                             rot.push_back(
                                                  path[t] );
                                        for ( int t = 0;
                                             t < in_path; t++ )
                                             rot.push_back( path[t] );
                                        path = rot;    }
                                   break;    }    }
                         if ( appends[cur] >= 8 )
                              break;    // degenerate loop: stop
                         n_repeat_cross++;    }
                    else
                    {    visited.insert(cur);
                         if ( inv[cur] != cur )
                              visited.insert(inv[cur]);    }
                    appends[cur]++;
                    path.push_back(cur);
                    region_pass = was_visited;
                    // Track the entry of the current short run (for
                    // multi-edge repeat pairing).
                    if ( edges[cur].isize( ) < 2000 )
                    {    if ( region_entry < 0 ) region_entry = prev_e;    }
                    else region_entry = -1;
                    // Successor choice: exits paired with our entry
                    // by spanning reads take priority (each genomic
                    // copy takes its own continuation); otherwise
                    // the best global vote.
                    int nxt = -1;
                    {    int64_t best = 0;
                         auto pit = pair_votes.find(
                              { cur, prev_e } );
                         auto pit2 = pair_votes.find(
                              { cur, region_entry } );
                         auto it = trans_votes.find(cur);
                         if ( it != trans_votes.end( ) )
                              for ( auto const& s : it->second )
                              {    if ( s.second < min_support )
                                        continue;
                                   if ( exit_use[cur][s.first]
                                        >= 2 * cap ) continue;
                                   int64_t sc = s.second;
                                   if ( pit != pair_votes.end( ) )
                                   {    auto pv = pit->second.find(
                                             s.first );
                                        if ( pv != pit->second.end( ) )
                                             sc += int64_t(1000000)
                                                  * pv->second;    }
                                   if ( pit2 != pair_votes.end( ) )
                                   {    auto pv = pit2->second.find(
                                             s.first );
                                        if ( pv != pit2->second.end( ) )
                                             sc += int64_t(1000000)
                                                  * pv->second;    }
                                   if ( sc > best )
                                   {    best = sc;
                                        nxt = s.first;    }    }
                         if ( nxt >= 0 )
                              exit_use[cur][nxt]++;    }
                    if ( nxt < 0 )
                    {    // No read evidence here: follow the graph
                         // when it is unambiguous (non-branching
                         // vertex), or skip through a variant bubble
                         // whose branches reconverge.
                         int w = to_right[cur];
                         if ( hb.FromSize(w) == 1 )
                         {    nxt = hb.IFrom( w, 0 );
                              n_fallback++;    }
                         else if ( hb.FromSize(w) >= 2
                              && bubble_skip( cur, nxt ) )
                              n_bubble++;
                         if ( nxt < 0 && hb.FromSize(w) >= 2
                              && guess_branch( cur, nxt ) )
                              n_guess++;
                         if ( nxt < 0 ) break;    }
                    prev_e = cur;
                    cur = nxt;    }
               // Backward extension: prepend predecessors (voted, or
               // the unique graph predecessor at a non-branching
               // vertex).  This captures the segment upstream of the
               // start edge in the same path.
               if ( !path.empty( ) && !circ )
               {    vec<int> pre;
                    while ( true )
                    {    int f = path.front( );
                         int p = -1;
                         auto it2 = pred.find(f);
                         if ( it2 != pred.end( ) ) p = it2->second;
                         else
                         {    int v = to_left[f];
                              if ( hb.ToSize(v) == 1 )
                              {    p = hb.ITo( v, 0 );
                                   n_bfallback++;    }    }
                         if ( p < 0 ) break;
                         if ( visited.count(p) )
                         {    if ( p == path.back( ) )
                                   circ = True;
                              break;    }
                         visited.insert(p);
                         if ( inv[p] != p ) visited.insert(inv[p]);
                         pre.push_back(p);    }
                    if ( !pre.empty( ) )
                    {    Reverse(pre);
                         pre.append(path);
                         path = pre;    }
                    // Circle check: does the terminus connect back
                    // to the front?
                    if ( !circ && path.size( ) > 1 )
                    {    int b = path.back( );
                         auto it3 = succ.find(b);
                         if ( ( it3 != succ.end( )
                                   && it3->second == path.front( ) )
                              || to_right[b]
                                   == to_left[ path.front( ) ] )
                              circ = True;    }    }
               if ( !path.empty( ) )
               {    paths.push_back(path);
                    circular.push_back(circ);    }    }    }

     // ---- Circle extension ----

     // Grow the longest paths by absorbing donor suffixes: the
     // receiver's terminal edge votes for an edge inside a donor, and
     // the donor's suffix from that edge is appended (the donor keeps
     // its prefix).  The receiver only appends, donors only shrink,
     // so the process always terminates; a receiver whose terminus
     // votes for its own front closes a circular scaffold.

     {    vec<int> order;
          for ( int i = 0; i < paths.isize( ); i++ )
               order.push_back(i);
          std::sort( order.begin( ), order.end( ),
               [&]( int a, int b )
               {    int64_t la = 0, lb = 0;
                    for ( int e : paths[a] ) la += edges[e].isize( );
                    for ( int e : paths[b] ) lb += edges[e].isize( );
                    return la > lb;    } );
          std::set<int> dead, closed;
          int64_t n_merge = 0;
          for ( int ri : order )
          {    if ( dead.count(ri) || closed.count(ri)
                    || paths[ri].empty( ) ) continue;
               int64_t guard = 0;
               while ( guard++ < 4000 && !closed.count(ri)
                    && !paths[ri].empty( ) )
               {    int64_t before = n_merge;
                    int t = paths[ri].back( );
                    auto it = trans_votes.find(t);
                    if ( it != trans_votes.end( ) )
                    {    vec< std::pair<int,int> > cand;
                         for ( auto const& s : it->second )
                              if ( s.second >= min_support )
                                   cand.push_back(
                                        { s.first, s.second } );
                         std::sort( cand.begin( ), cand.end( ),
                              []( std::pair<int,int> const& a,
                                   std::pair<int,int> const& b )
                              {    return a.second > b.second;    } );
                         for ( auto const& c : cand )
                         {    int x = c.first;
                              // Self-landing (front or interior):
                              // the path is a closed cycle; rotate
                              // so the landing edge is first.
                              int selfp = -1;
                              for ( int t2 = 0;
                                   t2 < paths[ri].isize( ); t2++ )
                                   if ( paths[ri][t2] == x )
                                   {    selfp = t2;    break;    }
                              if ( selfp >= 0 )
                              {    closed.insert(ri);
                                   circular[ri] = True;
                                   if ( selfp > 0 )
                                   {    vec<int> rot;
                                        for ( int t2 = selfp;
                                             t2 < paths[ri].isize( );
                                             t2++ )
                                             rot.push_back(
                                                  paths[ri][t2] );
                                        for ( int t2 = 0;
                                             t2 < selfp; t2++ )
                                             rot.push_back(
                                                  paths[ri][t2] );
                                        paths[ri] = rot;    }
                                   if ( VERBOSE )
                                        cout << "  ext: receiver "
                                             << ri << " closes at "
                                             << x << " (" << c.second
                                             << " votes)" << endl;
                                   break;    }
                              Bool done = False;
                              for ( int j = 0;
                                   j < paths.isize( ); j++ )
                              {    if ( j == ri || dead.count(j)
                                        || closed.count(j)
                                        || paths[j].empty( ) )
                                        continue;
                                   int p = -1;
                                   for ( int u = 0;
                                        u < paths[j].isize( ); u++ )
                                        if ( paths[j][u] == x )
                                        {    p = u;    break;    }
                                   if ( p < 0 )
                                   {    if ( VERBOSE && c.second >= 5 )
                                             cout << "  ext: receiver "
                                                  << ri << " wants "
                                                  << x << " ("
                                                  << c.second
                                                  << " votes): "
                                                  << "not in live "
                                                  << "paths"
                                                  << endl;
                                        continue;    }
                                   // The absorbed suffix must not
                                   // overlap the receiver besides
                                   // the joint edge itself.
                                   std::set<int> redges(
                                        paths[ri].begin( ),
                                        paths[ri].end( ) );
                                   Bool overlap = False;
                                   for ( int u = p;
                                        u < paths[j].isize( ); u++ )
                                   {    int e = paths[j][u];
                                        if ( redges.count(e)
                                             && e != x )
                                        {    overlap = True;
                                             break;    }    }
                                   if ( overlap )
                                   {    if ( VERBOSE && c.second >= 5 )
                                             cout << "  ext: receiver "
                                                  << ri << " wants "
                                                  << x << " ("
                                                  << c.second
                                                  << " votes): "
                                                  << "overlap in "
                                                  << "path " << j
                                                  << endl;
                                        continue;    }
                                   for ( int u = p;
                                        u < paths[j].isize( ); u++ )
                                        paths[ri].push_back(
                                             paths[j][u] );
                                   if ( p == 0 ) dead.insert(j);
                                   else paths[j].resize(p);
                                   n_merge++;
                                   done = True;
                                   break;    }
                              if ( done || closed.count(ri) )
                                   break;    }    }
                    if ( closed.count(ri) ) break;
                    if ( n_merge == before ) break;   // no progress
               }    }
          // Drop absorbed-away paths.
          vec< vec<int> > live;
          vec<Bool> live_circ;
          for ( int i = 0; i < paths.isize( ); i++ )
               if ( !dead.count(i) )
               {    live.push_back( paths[i] );
                    live_circ.push_back( circular[i] );    }
          paths = live;
          circular = live_circ;
          cout << Date( ) << ": circle extension: " << n_merge
               << " joins" << endl;
          round_joins += n_merge;    }
          if ( paths.isize( ) == paths_before
               && round_joins == joins_before )
               break;    }

     cout << Date( ) << ": built " << paths.size( )
          << " scaffold paths (" << n_fallback
          << " graph-successor steps, " << n_bubble
          << " bubble skips, " << n_guess << " branch guesses, "
          << n_repeat_cross << " repeat crossings)" << endl;


     std::map< std::pair<int,int>, vec<String> > gap_samples;
     // ---- 6b. P4 phase 2: ONT-consensus gap sequences ----

     // For every junction with >= 2 flank-anchored read samples,
     // extract the reads' gap subsequences and take a column
     // majority.  The consensus both validates graph-threaded
     // fills (phase 1 accepts a candidate only if it is consistent
     // with the reads) and fills the junctions the graph cannot.
     std::map< std::pair<int,int>, String > gap_cons;
     {    // Which reads carry which junction spans?  A junction's
          // spanning reads come in BOTH orientations: reads that
          // threaded (a -> b) directly, and rc reads that threaded
          // the mirror junction (inv[b] -> inv[a]) -- the mirror
          // samples are reverse-complemented and pooled, so every
          // gap contributes all of its spanning reads.
          struct GapWant { int a, b, r1, r2; Bool rc; };
          std::map< String, vec<GapWant> > want;
          for ( auto const& kv : gap_reads )
          {    int a = kv.first.first, b = kv.first.second;
               int nd = kv.second.size( ), nm = 0;
               auto mit = gap_reads.find( { inv[b], inv[a] } );
               if ( mit != gap_reads.end( ) )
                    nm = mit->second.size( );
               if ( nd + nm < 2 ) continue;
               for ( auto const& s : kv.second )
               {    GapWant w = { a, b, s.second.first,
                         s.second.second, False };
                    want[ s.first ].push_back(w);    }
               if ( mit != gap_reads.end( ) )
                    for ( auto const& s : mit->second )
                    {    GapWant w = { a, b, s.second.first,
                              s.second.second, True };
                         want[ s.first ].push_back(w);    }    }
          if ( !want.empty( ) )
          {    // Stream the reads once and pull the spans.
               String cmd2 = "cat " + reads_file;
               if ( reads_file.Contains( ".gz", -1 ) )
                    cmd2 = "zcat " + reads_file;
               fast_pipe_ifstream in2(cmd2);
               auto rc_in_place = []( String& s )
               {    int L2 = s.isize( );
                    for ( int z = 0; z < L2 / 2; z++ )
                    {    char c1 = s[z], c2 = s[ L2 - 1 - z ];
                         s[z] = ( c2=='A' ? 'T' : c2=='T' ? 'A'
                              : c2=='C' ? 'G' : c2=='G' ? 'C' : 'N' );
                         s[ L2 - 1 - z ] =
                              ( c1=='A' ? 'T' : c1=='T' ? 'A'
                              : c1=='C' ? 'G' : c1=='G' ? 'C'
                              : 'N' );    }
                    if ( L2 % 2 == 1 )
                         s[ L2 / 2 ] =
                              ( s[ L2 / 2 ]=='A' ? 'T'
                              : s[ L2 / 2 ]=='T' ? 'A'
                              : s[ L2 / 2 ]=='C' ? 'G'
                              : s[ L2 / 2 ]=='G' ? 'C' : 'N' );    };
               String line2, hdr, seq2;
               while ( true )
               {    getline( in2, line2 );
                    if ( in2.fail( ) ) break;
                    if ( line2.empty( ) ) continue;
                    if ( line2[0] == '>' )
                    {    if ( !hdr.empty( ) && want.count(hdr) )
                              for ( auto const& w : want[hdr] )
                              {    String sub = seq2.substr(
                                        w.r1, w.r2 - w.r1 );
                                   if ( w.rc ) rc_in_place(sub);
                                   gap_samples[ { w.a, w.b } ]
                                        .push_back(sub);    }
                         hdr = line2.substr( 1 );
                         size_t sp = hdr.find( ' ' );
                         if ( sp != String::npos )
                              hdr.resize(sp);
                         seq2.resize(0);    }
                    else if ( line2[0] != '@' )
                         seq2 += line2;    }
               if ( !hdr.empty( ) && want.count(hdr) )
                    for ( auto const& w : want[hdr] )
                    {    String sub = seq2.substr(
                              w.r1, w.r2 - w.r1 );
                         if ( w.rc ) rc_in_place(sub);
                         gap_samples[ { w.a, w.b } ]
                              .push_back(sub);    }    }    }
     // Consensus per junction via the internal MultipleAligner (the
     // same engine AssembleGaps uses for PE gap consensus): the
     // longest sample seeds the alignment, every other sample is
     // aligned to the evolving consensus with a phred-scored ONT
     // error model (subst/del/ins), so the samples' indel drift is
     // corrected rather than trusted as column offsets -- a plain
     // column majority scored only ~59% identity on the E. coli
     // fills because each read's flank-alignment end drifts by its
     // own accumulated indels.
     for ( auto const& kv : gap_samples )
     {    vec<String> const& smp = kv.second;
          if ( smp.size( ) < 2 ) continue;
          int best = 0;
          for ( int i = 1; i < smp.isize( ); i++ )
               if ( smp[i].isize( ) > smp[best].isize( ) )
                    best = i;
          if ( smp[best].isize( ) < 100 ) continue;
          // sanitize: ACGT only
          bvec cons0;
          for ( int i = 0; i < smp[best].isize( ); i++ )
          {    char c = smp[best][i];
               if ( c == 'A' ) cons0.push_back(0);
               else if ( c == 'C' ) cons0.push_back(1);
               else if ( c == 'G' ) cons0.push_back(2);
               else if ( c == 'T' ) cons0.push_back(3);    }
          if ( cons0.size( ) < 100 ) continue;
          MultipleAligner ma(
               Scorer( 0.005, 0.008, 0.003 ), gvec( cons0 ) );
          for ( int i = 0; i < smp.isize( ); i++ )
          {    if ( i == best ) continue;
               bvec rd;
               for ( int j = 0; j < smp[i].isize( ); j++ )
               {    char c = smp[i][j];
                    if ( c == 'A' ) rd.push_back(0);
                    else if ( c == 'C' ) rd.push_back(1);
                    else if ( c == 'G' ) rd.push_back(2);
                    else if ( c == 'T' ) rd.push_back(3);    }
               if ( rd.size( ) < 100 ) continue;
               ma.addRead( rd.begin( ), rd.end( ) );    }
          bvec cons = ma.getConsensus( );
          String out;
          out.reserve( cons.size( ) );
          for ( size_t z = 0; z < cons.size( ); z++ )
               out.push_back( Base::val2Char( cons[z] ) );
          if ( out.isize( ) >= 100 )
               gap_cons[ kv.first ] = out;    }
     cout << Date( ) << ": gap consensus: "
          << gap_cons.size( ) << " junctions with read-backed gap"
          << " sequences" << endl;
     // ---- 6c. P4 phase 1: gap threading ----

     // The coarse waypoints make the walk jump over the small edges
     // between backbone neighbors; the graph holds most of that
     // sequence.  For each jump with a read gap estimate, search
     // the graph for a connecting interior path of the right size
     // and splice it in.  Jumps with no consistent path stay as
     // N-gaps (true gaps, reported for ONT-consensus filling).
     int64_t n_gaps_filled = 0, n_gaps_kept = 0, n_bp_filled = 0,
          n_gaps_ambig = 0;
     {    Ofstream( gout, out_head + ".gaps.txt" );
          gout << "# left_edge right_edge gap_mean gap_sd n_reads"
               << " (kept N-gaps: ONT-consensus targets)\n";
          // Duplication guard: 19-mer codes of the deep flanks
          // (outside the K-1 vertex-overlap zones, which are shared
          // by construction).  A candidate interior sharing a
          // sampled 19-mer with a flank duplicates placed
          // sequence -- e.g. an inverted-repeat side copy -- and is
          // rejected even when its size matches the read gap.
          auto kmer_codes = [&]( bvec const& bv, int from, int to,
               std::set<uint64_t>& codes )
          {    int const KK = 19;
               if ( to - from < KK ) return;
               uint64_t mask = ( uint64_t(1) << (2*KK) ) - 1;
               uint64_t c = 0;
               for ( int j = from; j < to; j++ )
               {    c = ( ( c << 2 ) | bv[j] ) & mask;
                    if ( j - from >= KK - 1 ) codes.insert(c);    }    };
          for ( int pi = 0; pi < paths.isize( ); pi++ )
          {    vec<int>& path = paths[pi];
               std::set<int> pset( path.begin( ), path.end( ) );
               vec<int> out;
               for ( int i = 0; i < path.isize( ); i++ )
               {    if ( i > 0 )
                    {    int a = out.back( ), b = path[i];
                         if ( to_right[a] != to_left[b] )
                         {    // Jump: read gap estimate.
                              int gap = 100;
                              auto jit = jump_gap.find( { a, b } );
                              if ( jit != jump_gap.end( )
                                   && jit->second.count > 0 )
                                   gap = Max( 1, int( double(
                                        jit->second.gap_sum)
                                        / jit->second.count ) );
                              // Read gap estimates are precise
                              // (sd <= 12 bp on E. coli), so the
                              // acceptance slop only needs to cover
                              // ONT indel drift over the gap span.
                              int slop = 100 + gap / 10;
                              int64_t budget =
                                   gap + slop + OVER;
                              // Depth-first enumeration of interior
                              // graph paths a -> ... -> b.
                              vec< vec<int> > cands;
                              int64_t explored = 0;
                              std::function<void(int, vec<int>,
                                   int64_t)> dfs =
                                   [&]( int v, vec<int> ip,
                                        int64_t sz )
                              {    if ( cands.isize( ) >= 64 )
                                        return;
                                   if ( ++explored > 100000 )
                                        return;
                                   if ( v == to_left[b] )
                                   {    if ( ip.isize( ) > 0
                                             && abs( int( sz - OVER
                                                  - gap ) )
                                                  <= slop )
                                             cands.push_back(ip);
                                        return;    }
                                   if ( hb.FromSize(v) == 0 )
                                        return;
                                   for ( int j = 0;
                                        j < hb.FromSize(v); j++ )
                                   {    int e2 = hb.IFrom( v, j );
                                        if ( pset.count(e2) )
                                             continue;
                                        // Also skip mirror copies:
                                        // the doubled graph holds
                                        // both orientations, and an
                                        // interior containing the
                                        // involution of a placed
                                        // edge would duplicate its
                                        // sequence (inverted).
                                        if ( pset.count( inv[e2] ) )
                                             continue;
                                        if ( edge_class[e2]
                                                  == EDGE_REPEAT )
                                             continue;
                                        int64_t sz2 = sz
                                             + edges[e2].isize( )
                                             - OVER;
                                        if ( sz2 > budget )
                                             continue;
                                        if ( ip.isize( ) >= 24 )
                                             continue;
                                        vec<int> ip2 = ip;
                                        ip2.push_back(e2);
                                        dfs( to_right[e2], ip2,
                                             sz2 );    }    };
                              dfs( to_right[a], vec<int>( ), 0 );
                              if ( VERBOSE )
                              {    cout << "  fill e" << a
                                        << " -> e" << b
                                        << " gap=" << gap
                                        << " cands="
                                        << cands.isize( );
                                   for ( int c = 0;
                                        c < cands.isize( ); c++ )
                                   {    cout << " [";
                                        for ( int e2 : cands[c] )
                                             cout << " " << e2
                                                  << "("
                                                  << edges[e2]
                                                       .isize( )
                                                  << ")";
                                        cout << " ]";    }
                                   cout << endl;    }
                              if ( cands.isize( ) > 0 )
                              {    // Duplication guard first.
                                   std::set<uint64_t> dup;
                                   {    bvec const& va = edges[a];
                                        int la = va.isize( );
                                        kmer_codes( va,
                                             Max( 0, la - 1000 ),
                                             Max( 0, la - 250 ),
                                             dup );
                                        bvec const& vb = edges[b];
                                        int lb = vb.isize( );
                                        kmer_codes( vb,
                                             Min( 250, lb ),
                                             Min( 1000, lb ),
                                             dup );    }
                                   vec< vec<int> > pass;
                                   for ( int c = 0;
                                        c < cands.isize( ); c++ )
                                   {    Bool dup_hit = False;
                                        for ( int e2 : cands[c] )
                                        {    bvec const& v =
                                                  edges[e2];
                                             int le = v.isize( );
                                             int from =
                                                  Min( 250, le );
                                             int to =
                                                  Max( from,
                                                       le - 250 );
                                             int step = Max( 1,
                                                  ( to - from )
                                                  / 10 );
                                             uint64_t mask =
                                                  ( uint64_t(1)
                                                  << 38 ) - 1;
                                             for ( int p = from;
                                                  p + 19 <= to;
                                                  p += step )
                                             {    uint64_t c2 = 0;
                                                  for ( int j = p;
                                                       j < p + 19;
                                                       j++ )
                                                       c2 =
                                                       ( ( c2 << 2 )
                                                       | v[j] )
                                                       & mask;
                                                  if ( dup.count(
                                                            c2 ) )
                                                  {    dup_hit =
                                                            True;
                                                       break;    }    }
                                             if ( dup_hit )
                                                  break;    }
                                        if ( !dup_hit )
                                             pass.push_back(
                                                  cands[c] );    }
                                   // Reads-consistency: the junction's
                                   // ONT consensus is ground truth; a
                                   // candidate route whose sequence
                                   // contradicts it (e.g. the wrong
                                   // arm's copy of an inverted
                                   // repeat) is rejected even when
                                   // it is the only candidate -- the
                                   // gap is then filled from the
                                   // reads themselves.
                                   {    auto cit = gap_cons.find(
                                             { a, b } );
                                        if ( VERBOSE && cit
                                             != gap_cons.end( ) )
                                             cout << "  cons e" << a
                                                  << "->e" << b
                                                  << " len="
                                                  << cit->second
                                                       .isize( )
                                                  << endl;
                                        if ( cit != gap_cons.end( )
                                             && pass.isize( ) > 0 )
                                        {    String const& cons =
                                                  cit->second;
                                             vec<uint64_t> probes;
                                             uint64_t mask =
                                                  ( uint64_t(1)
                                                  << 38 ) - 1;
                                             for ( int p = 0;
                                                  p + 19
                                                       <= cons
                                                       .isize( );
                                                  p += Max( 20,
                                                       cons.isize( )
                                                       / 20 ) )
                                             {    Bool ok2 = True;
                                                  uint64_t c2 = 0;
                                                  for ( int j = p;
                                                       j < p + 19;
                                                       j++ )
                                                  {    char ch =
                                                            cons[j];
                                                       if ( ch
                                                            == 'A' )
                                                            c2 =
                                                            ( c2 << 2 );
                                                       else if ( ch
                                                            == 'C' )
                                                            c2 =
                                                            ( c2 << 2 )
                                                            | 1;
                                                       else if ( ch
                                                            == 'G' )
                                                            c2 =
                                                            ( c2 << 2 )
                                                            | 2;
                                                       else if ( ch
                                                            == 'T' )
                                                            c2 =
                                                            ( c2 << 2 )
                                                            | 3;
                                                       else
                                                       {    ok2 =
                                                                 False;
                                                            break;    }    }
                                                  if ( ok2 )
                                                       probes
                                                       .push_back(
                                                            c2 );    }
                                             if ( probes.isize( )
                                                       >= 5 )
                                             {    vec< vec<int> >
                                                        ok;
                                                  for ( int c = 0;
                                                       c <
                                                       pass.isize( );
                                                       c++ )
                                                  {    std::set<
                                                            uint64_t >
                                                            cset;
                                                       for ( int e2 :
                                                            pass[c] )
                                                            kmer_codes(
                                                                 edges[
                                                                 e2 ],
                                                                 0,
                                                                 edges[
                                                                 e2 ]
                                                                 .isize( ),
                                                                 cset );
                                                       int hits = 0;
                                                       for ( auto c2 :
                                                            probes )
                                                            if ( cset
                                                                 .count(
                                                                 c2 ) )
                                                                 hits++;
                                                       if ( hits * 2
                                                            >= (int)
                                                                 probes
                                                                 .isize( ) )
                                                            ok
                                                            .push_back(
                                                                 pass[
                                                                 c ] );    }
                                                  pass = ok;    }    }    }
                                   if ( pass.isize( ) == 0 )
                                   {    n_gaps_kept++;
                                        double sd = 0.0;
                                        if ( jit != jump_gap.end( )
                                             && jit->second.count
                                                  > 0 )
                                        {    double mean =
                                                  double(
                                                  jit->second.gap_sum)
                                                  / jit->second.count;
                                             double var =
                                                  double(
                                                  jit->second.gap_sq)
                                                  / jit->second.count
                                                  - mean * mean;
                                             sd = var > 0
                                                  ? std::sqrt(var)
                                                  : 0.0;
                                             gout << a << " " << b
                                                  << " " << mean
                                                  << " " << sd
                                                  << " "
                                                  << jit->second
                                                       .count
                                                  << "\n";    }
                                        else
                                             gout << a << " " << b
                                                  << " " << gap
                                                  << " 0 0\n";    }
                                   else
                                   {    cands = pass;
                                   // Pick the interior whose chain
                                   // size best matches the read gap
                                   // (ties: fewest edges).
                                   int best = 0, bestd = 0;
                                   for ( int c = 0;
                                        c < cands.isize( ); c++ )
                                   {    int64_t sz = 0;
                                        for ( int e2 : cands[c] )
                                             sz += edges[e2].isize( )
                                                  - OVER;
                                        int d = abs( int(
                                             sz - OVER - gap ) );
                                        if ( c == 0 || d < bestd
                                             || ( d == bestd
                                                  && cands[c]
                                                       .isize( )
                                                       < cands[best]
                                                       .isize( ) ) )
                                        {    bestd = d;
                                             best = c;    }    }
                                   // Mirror-route ambiguity: if the
                                   // candidate set contains two
                                   // interiors that are involution
                                   // images of each other (a fwd/rc
                                   // route pair through an inverted
                                   // repeat), the graph cannot tell
                                   // which one is placed -- keep
                                   // the N-gap for read-based
                                   // filling.  Variant-bubble
                                   // haplotypes are not mirrors and
                                   // are filled normally.
                                   Bool mirror_pair = False;
                                   for ( int c1 = 0;
                                        c1 < cands.isize( ); c1++ )
                                   {    for ( int c2 = c1 + 1;
                                             c2 < cands.isize( );
                                             c2++ )
                                        {    int k =
                                                  cands[c1].isize( );
                                             if ( k
                                                  != cands[c2]
                                                       .isize( ) )
                                                  continue;
                                             Bool m = True;
                                             for ( int t = 0;
                                                  t < k; t++ )
                                                  if ( cands[c2]
                                                       [ k-1-t ]
                                                       != inv[
                                                       cands[c1]
                                                       [ t ] ] )
                                                  {    m = False;
                                                       break;    }
                                             if ( m )
                                             {    mirror_pair =
                                                       True;
                                                  break;    }    }
                                        if ( mirror_pair ) break;    }
                                   if ( mirror_pair )
                                   {    // Ambiguous fwd/rc route:
                                        // keep the N-gap for
                                        // read-based filling.
                                        n_gaps_ambig++;
                                        if ( jit != jump_gap.end( )
                                             && jit->second.count
                                                  > 0 )
                                        {    double mean =
                                                  double(
                                                  jit->second.gap_sum)
                                                  / jit->second.count;
                                             double var =
                                                  double(
                                                  jit->second.gap_sq)
                                                  / jit->second.count
                                                  - mean * mean;
                                             double sd = var > 0
                                                  ? std::sqrt(var)
                                                  : 0.0;
                                             gout << a << " " << b
                                                  << " " << mean
                                                  << " " << sd
                                                  << " "
                                                  << jit->second
                                                       .count
                                                  << " ambiguous"
                                                  << "\n";    }
                                        else
                                             gout << a << " " << b
                                                  << " " << gap
                                                  << " 0 0"
                                                  << " ambiguous"
                                                  << "\n";    }
                                   else
                                   {
                                   for ( int e2 : cands[best] )
                                   {    out.push_back(e2);
                                        n_bp_filled +=
                                             edges[e2].isize( )
                                             - OVER;    }
                                   n_gaps_filled++;    }    }    }
                              else
                              {    n_gaps_kept++;
                                   double sd = 0.0;
                                   if ( jit != jump_gap.end( )
                                        && jit->second.count > 0 )
                                   {    double mean = double(
                                             jit->second.gap_sum)
                                             / jit->second.count;
                                        double var =
                                             double(
                                             jit->second.gap_sq)
                                             / jit->second.count
                                             - mean * mean;
                                        sd = var > 0
                                             ? std::sqrt(var)
                                             : 0.0;
                                        gout << a << " " << b
                                             << " " << mean
                                             << " " << sd << " "
                                             << jit->second.count
                                             << "\n";    }
                                   else
                                        gout << a << " " << b
                                             << " " << gap
                                             << " 0 0\n";    }    }    }
                    out.push_back( path[i] );    }
               path = out;    }    }
     cout << Date( ) << ": gap threading: filled "
          << n_gaps_filled << " gaps ("
          << ToStringAddCommas(n_bp_filled)
          << " bp from the graph), " << n_gaps_kept
          << " kept as N-gaps (" << n_gaps_ambig
          << " ambiguous routes)" << endl;
     // Path termini diagnostics: why did each walk stop?
     if ( VERBOSE )
     {    Ofstream( pout, out_head + ".paths.txt" );
          for ( int i = 0; i < paths.isize( ); i++ )
          {    vec<int> const& path = paths[i];
               if ( path.empty( ) ) continue;
               pout << "path " << i << " (" << path.size( )
                    << " edges, circular="
                    << ( circular[i] ? "yes" : "no" ) << "):";
               for ( int e : path )
                    pout << " " << e << "("
                         << edges[e].isize( ) << ")";
               pout << "\n";
               int t = path.back( );
               int w = to_right[t];
               cout << "  path " << i << " (" << path.size( )
                    << " edges) ends at edge " << t << " (len "
                    << edges[t].isize( ) << ", class "
                    << ( edge_class[t] == EDGE_NORMAL ? "N"
                         : edge_class[t] == EDGE_SHORT ? "S" : "R" )
                    << "), out-degree " << hb.FromSize(w)
                    << ", succ votes:";
               auto it = trans_votes.find(t);
               if ( it != trans_votes.end( ) )
                    for ( auto const& s : it->second )
                         cout << " " << s.first << ":" << s.second;
               int f = path.front( );
               int v = to_left[f];
               cout << " | starts at edge " << f << " (len "
                    << edges[f].isize( ) << "), in-degree "
                    << hb.ToSize(v) << ", pred votes:";
               for ( int j = 0; j < hb.ToSize(v); j++ )
               {    int pe = hb.ITo( v, j );
                    auto it2 = trans_votes.find(pe);
                    int votes = 0;
                    if ( it2 != trans_votes.end( ) )
                    {    auto vs = it2->second.find(f);
                         if ( vs != it2->second.end( ) )
                              votes = vs->second;    }
                    cout << " " << pe << ":" << votes;    }
               cout << endl;    }    }

     // ---- 7. Output ----

     {    Ofstream( fout, out_head + ".scaffold.fasta" );
          Ofstream( fills, out_head + ".ontfills.fasta" );
          vec<int> scaffold_lens;
          int n_circ = 0, n_jumps = 0, n_gaps_ont = 0;
          int64_t n_bp_ont = 0;
          int64_t n_gap_bp = 0;
          for ( int i = 0; i < paths.isize( ); i++ )
          {    vec<int> const& path = paths[i];
               if ( path.empty( ) ) continue;
               // Assemble the sequence: graph-adjacent edges overlap
               // by K-1 bases at their shared vertex; jumps get the
               // read-estimated gap as N's (P4 will fill these).
               String seq_out;
               for ( int j = 0; j < path.isize( ); j++ )
               {    int e = path[j];
                    bvec const& bv = edges[e];
                    size_t from = 0;
                    if ( j > 0 && to_right[ path[j-1] ] == to_left[e] )
                         from = ( bv.size( ) < size_t(OVER) )
                              ? bv.size( ) : size_t(OVER);
                    for ( size_t b = from; b < bv.size( ); b++ )
                         seq_out.push_back(
                              Base::val2Char( bv[b] ) );
                    if ( j < path.isize( ) - 1
                         && to_right[e] != to_left[ path[j+1] ] )
                    {    n_jumps++;
                         // P4 phase 2: an ONT-consensus gap sequence
                         // takes precedence over N's.
                         auto cit = gap_cons.find(
                              { e, path[j+1] } );
                         if ( cit != gap_cons.end( ) )
                         {    n_gaps_ont++;
                              n_bp_ont += cit->second.isize( );
                              seq_out.append( cit->second );
                              // Record the fill for consensus
                              // quality evaluation.
                              fills << ">ontfill_e" << e
                                    << "_e" << path[j+1]
                                    << " len="
                                    << cit->second.isize( )
                                    << "\n";
                              for ( int z = 0;
                                   z < cit->second.isize( );
                                   z += 70 )
                                   fills << cit->second
                                        .substr( z, 70 )
                                        << "\n";    }
                         else
                         {    if ( VERBOSE )
                                   cout << "  N-gap e" << e
                                        << " -> e" << path[j+1]
                                        << ": cons="
                                        << ( gap_cons.count(
                                                  { e, path[j+1] } )
                                                  ? "yes" : "no" )
                                        << " samples="
                                        << ( gap_reads.count(
                                                  { e, path[j+1] } )
                                                  ? (int) gap_reads[
                                                       { e,
                                                       path[j+1] } ]
                                                       .size( )
                                                  : 0 )
                                        << endl;    int gap = 100;
                              auto it = jump_gap.find(
                                   { e, path[j+1] } );
                              if ( it != jump_gap.end( )
                                   && it->second.count > 0 )
                                   gap = Max( 1, int( double(
                                        it->second.gap_sum)
                                        / it->second.count ) );
                              gap = Min( gap, 20000 );
                              n_gap_bp += gap;
                              seq_out.append(
                                   String( gap, 'N' ) );    }    }    }
               fout << ">scaffold_" << i << " edges=" << path.size( )
                    << " len=" << seq_out.size( )
                    << ( circular[i] ? " circular=True" : "" )
                    << "\n";
               for ( int b = 0; b < seq_out.isize( ); b += 70 )
                    fout << seq_out.substr( b, 70 ) << "\n";
               scaffold_lens.push_back( seq_out.isize( ) );
               if ( circular[i] ) n_circ++;    }

           Sort( scaffold_lens );
           int64_t scaffold_total = 0;
           for ( int l : scaffold_lens ) scaffold_total += l;
           int64_t n50 = 0;
           {    int64_t acc = 0;
                for ( int i = scaffold_lens.isize( ) - 1; i >= 0; i-- )
                {    acc += scaffold_lens[i];
                     if ( acc >= scaffold_total / 2 )
                     {    n50 = scaffold_lens[i];
                          break;    }    }    }
           cout << "\n===================================================="
                << "=======================\n";
           cout << "P3 scaffold stats (read threading)\n";
           cout << "  repeats detected: " << n_repeat
                << " (threshold " << repeat_threshold << ")\n";
           cout << "  short edges (< " << min_line << " bp): " << n_short
                << " (threaded through, not scaffold nodes)\n";
           cout << "  hi-q reads: " << n_hiq_reads << " / " << n_ge
                << "\n";
           cout << "  normal-jump links: " << link_map.size( ) << "\n";
           cout << "  thread successors: " << succ.size( ) << " ("
                << n_adj_succ << " adjacent, " << n_jump_succ
                << " jump)\n";
           cout << "  scaffold paths: " << paths.size( ) << " ("
                << n_circ << " circular)\n";
           cout << "  gap jumps in output: " << n_jumps << " ("
                << n_gaps_ont << " ONT-filled, "
                << ToStringAddCommas(n_bp_ont) << " bp; "
                << ToStringAddCommas(n_gap_bp) << " N bp)\n";
           cout << "  scaffold total: "
                << ToStringAddCommas(scaffold_total) << " bp\n";
           cout << "  scaffold N50: " << ToStringAddCommas(n50) << " bp\n";
           cout << "  PE-only N50: 133,298 bp\n";
           cout << "  improvement: "
                << fixed << setprecision(2)
                << double(n50) / 133298.0 << "x\n";
           cout << "=========================================================="
                << "=============\n"
                << resetiosflags(ios::fixed);    }
}
