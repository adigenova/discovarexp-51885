///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// See also: NhoodInfoServer

// This uses files in the assembly directory plus files in the parent directory:
// subsam.names
// genome.names
// genome.names_alt
// genome.ambint
// genome.ambint_alt

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "FastIfstream.h"
#include "TokenizeString.h"

#include "paths/long/large/tools/NhoodInfoCore.h"


void NhoodInfoCore( int argc, char *argv[] ) 
{
     // Note that some argument defaults are mirrored by 
     // nhood_info_state::Initialize.

    const char *DOC =
	"NhoodInfo is a visualization tool for DISCOVAR de novo assemblies. "
	"Regions of interest can be selected using reference co-ordinates, "
	"short alignable sequences, line ids, or edge ids. "
	"The annoated localized graph is generated as a GraphViz dot file, or "
	"optionally also as a PDF or PNG file.";

     BeginCommandArguments;
     CommandEnvVarPrefix("NHOODINFO");
     CommandDoc( DOC );
     CommandArgument_String_OrDefault_Doc(DIR_IN, ".", "location of the DISCOVAR "
          "final assembly directory containing a.hbx file");
     CommandArgument_String_Abbr_Doc(OUT, O, 
	  "output assembly graph files as OUT.{dot,fasta} ");
     CommandArgument_String_Abbr_OrDefault_Doc(OUT_PREFIX, OP, "",
          "If set, output files are OUT_PREFIX/OUT.{dot,fasta};\n"
          "Can also be set via environment variable NHOODINFO_OUT_PREFIX");
     CommandArgument_String_Abbr_OrDefault_Doc(SEEDS, S, "",
          "see online documentation");
     CommandArgument_Int_OrDefault_Doc(RANDOM_SEED, -1,
          "set seed for random number generator");
     CommandArgument_String_OrDefault_Doc(SEEDS_MINUS, "", 
          "exclude these seeds");
     CommandArgument_Int_Abbr_OrDefault_Doc(DEPTH, D, 5, 
          "expand to this depth");
     CommandArgument_Bool_OrDefault_Doc(SHOW_INV, False, 
          "show inv edge ids");
     CommandArgument_Bool_Abbr_OrDefault_Doc(GREEN, G, False,
          "color seeds green");
     CommandArgument_Bool_OrDefault_Doc(BNG, False, "show BNG restriction sites");
     CommandArgument_Bool_OrDefault_Doc(ALTREF, False,
          "show alignments to alternate reference if available");
     CommandArgument_Bool_OrDefault_Doc(SHOW_CN, True,
          "show predicted copy number for edges");
     CommandArgument_Bool_OrDefault_Doc(SHOW_ALIGN, True, "show edge alignments");
     CommandArgument_Bool_OrDefault_Doc(PURPLE1, False,
          "for multi-sample assemblies, color edges purple and bold if only the "
          "first sample is present for the given edge");
     CommandArgument_Bool_Abbr_OrDefault_Doc(EXT, E, False,
          "extend along lines");
     CommandArgument_Bool_OrDefault_Doc(COUNT, False, 
          "show number of reads supporting each edge");
     CommandArgument_Bool_Abbr_OrDefault_Doc(NEATO, N, False, 
          "use neato layout for DOT file (tends to form nice circles)");
     CommandArgument_Bool_Abbr_OrDefault_Doc(INTERACTIVE, I, False, 
          "enter interactive mode");
     CommandArgument_Bool_Abbr_OrDefault_Doc(LABEL_CONTIGS, LC, False, 
          "label contigs");
     CommandArgument_Double_Abbr_OrDefault_Doc(FONTSIZE, FS, 20,
          "edge label font size");
     CommandArgument_Double_OrDefault_Doc(SCALE, 1.3,
          "scale for edge width, arrowhead size and vertex diameter");
     CommandArgument_Bool_OrDefault_Doc(PNG, False, 
          "generate PNG file using dot. Requires dot 2.28 pr higher." );
     CommandArgument_Bool_OrDefault_Doc(PDF, False,
          "generate PDF file using dot. Requires dot 2.28 or higher.");
     CommandArgument_Bool_OrDefault_Doc(REL, False,
          "use relative numbering for graph edges");
     CommandArgument_Bool_OrDefault_Doc(LINENO, False, "show lines numbers");
     CommandArgument_Bool_OrDefault_Doc(SEQ_LOOKUP, False, 
          "build assembly lookup table to allow fast lookup by sequence - "
          "slow initialization;");
     CommandArgument_Bool_OrDefault_Doc(COV2, False, 
          "multiply last coverage by two; a temporary hack");
     EndCommandArguments;


     // Check dot executable, etc.

     if ( PDF || PNG ) TestDot( );

     if ( SEEDS == "" && !INTERACTIVE ) {
	 cout << "You need to specify SEEDS if you're not using "
	      << "INTERACTIVE mode." << endl;
          Scram(1);
     }

     if ( OUT_PREFIX != "" ) 
	 OUT = OUT_PREFIX + "/" + OUT;


     NhoodInfoEngine engine;

     char* cwd = get_current_dir_name( );
     cout << "running from " << cwd << endl;
     free(cwd);
     if (INTERACTIVE) 
	 cout << "(loading)" << endl;
     else 
	 cout << Date( ) << ": loading" << endl;

     engine.Initialize(DIR_IN, SEQ_LOOKUP, EXT, COV2);

     if (INTERACTIVE)
	 cout << "Interactive mode, type H for help, Ctrl-C to exit." << endl;

     // Set state.

     engine.SetState(SEEDS, DEPTH, EXT, COUNT, GREEN, BNG, PURPLE1, NEATO, REL, 
          LINENO, SHOW_INV, SHOW_CN, SHOW_ALIGN, FONTSIZE, SCALE, ALTREF, COV2);

     engine.RunAsClient(INTERACTIVE, DEPTH, OUT, PNG, PDF, RANDOM_SEED, 
          SEEDS_MINUS, LABEL_CONTIGS);

     // Done.

     cout << Date( ) << ": done" << endl;    
}

void NhoodInfoEngine::Initialize(const String& DIR_IN, const bool SEQ_LOOKUP, 
     const bool EXT, const bool COV2) 
{
     String head;
     if ( isReadable( DIR_IN + "/a.hbx" ) ) head = DIR_IN + "/a";
     else
     {    cout << "I can't find your assembly.  Please check the DIR_IN argument."
  	       << endl;
         Scram(1);    }
     BinaryReader::readFile( head + ".hbx", &hb );

     // Load inversion.

     BinaryReader::readFile( DIR_IN + "/a.inv", &inv );

     // Build lookup table.


     if ( SEQ_LOOKUP == true )
	 {    cout << Date( ) << ": building sequence lookup table" << endl;
	     vecbasevector edges( hb.E( ) );
	     for ( int e = 0; e < hb.E( ); e++ )
		 if ( e <= inv[e] ) edges[e] = hb.EdgeObject(e);
	     MakeKmerLookup0( edges, kmers_plus );
	     cout << Date( ) << ": done" << endl;    }

     // Load lines if available.

     if ( IsRegularFile( DIR_IN + "/a.lines" ) )
	 {    BinaryReader::readFile( head + ".lines", &lines );
	     GetTol( hb, lines, tol );
	     GetLineLengths( hb, lines, llens );    }
     else if (EXT)
	 {    cout << "EXT: can't find lines." << endl;
	     Scram(1);    }

     // Load genome record names if available.

     String gnf, gnf_alt, line;
     if ( IsRegularFile( DIR_IN + "/../genome.names" ) ) 
          gnf = DIR_IN + "/../genome.names";
     else if ( IsRegularFile( DIR_IN + "/genome.names" ) ) 
          gnf = DIR_IN + "/genome.names";
     if ( gnf != "" )
     {    fast_ifstream in(gnf);
          while(1)
	  {    getline( in, line );
	       if ( in.fail( ) ) break;
	       genome_names.push_back(line);    }    }
     if ( IsRegularFile( DIR_IN + "/../genome.names_alt" ) ) 
          gnf_alt = DIR_IN + "/../genome.names_alt";
     else if ( IsRegularFile( DIR_IN + "/genome.names_alt" ) ) 
          gnf_alt = DIR_IN + "/genome.names_alt";
     if ( gnf_alt != "" )
     {    fast_ifstream in(gnf_alt);
          while(1)
	  {    getline( in, line );
	       if ( in.fail( ) ) break;
	       genome_names_alt.push_back(line);    }    }

     // Load ambint if available.

     String xnf, xnf_alt;
     if ( IsRegularFile( DIR_IN + "/../genome.ambint" ) ) 
          xnf = DIR_IN + "/../genome.ambint";
     else if ( IsRegularFile( DIR_IN + "/genome.ambint" ) ) 
          xnf = DIR_IN + "/genome.ambint";
     if ( xnf != "" ) BinaryReader::readFile( xnf, &ambint );
     if ( IsRegularFile( DIR_IN + "/../genome.ambint_alt" ) ) 
          xnf_alt = DIR_IN + "/../genome.ambint_alt";
     else if ( IsRegularFile( DIR_IN + "/genome.ambint_alt" ) ) 
          xnf_alt = DIR_IN + "/genome.ambint_alt";
     if ( xnf_alt != "" ) BinaryReader::readFile( xnf, &ambint_alt );

     // Load alignments and coverage.

     if ( IsRegularFile( DIR_IN + "/a.aligns" ) )
          BinaryReader::readFile( DIR_IN + "/a.aligns", &hits );
     if ( IsRegularFile( DIR_IN + "/a.aligns_alt" ) )
          BinaryReader::readFile( DIR_IN + "/a.aligns_alt", &hits_alt );
     if ( IsRegularFile( DIR_IN + "/a.cov" ) )
	 BinaryReader::readFile( DIR_IN + "/a.cov", &cov );
     else if ( IsRegularFile( DIR_IN + "/a.covs" ) )
     {    BinaryReader::readFile( DIR_IN + "/a.covs", &covs );
          if (COV2)
          {    int ns = covs.size( );
               for ( int e = 0; e < covs[0].isize( ); e++ )
                    covs[ns-1][e].Set( covs[ns-1][e].Cov( ) * 2 );    }    }

     // Load edge counts.  Note that the file a.countsb is not required.  It is
     // in a.final but not a.fin.

     if ( IsRegularFile( DIR_IN + "/a.countsb" ) )
          BinaryReader::readFile( DIR_IN + "/a.countsb", &count );
     if ( IsRegularFile( DIR_IN + "/../subsam.names" ) )
          BinaryReader::readFile( DIR_IN + "/../subsam.names", &subsam_names );
     else if ( IsRegularFile( DIR_IN + "/subsam.names" ) )
          BinaryReader::readFile( DIR_IN + "/subsam.names", &subsam_names );
     else 
     {    cout << "I can't find the file subsam.names.  Perhaps you have "
               << "a very old assembly." << endl;
          cout << "Giving up." << endl;
          Scram(1);    }    }

void NhoodInfoEngine::SetState(const String& SEEDS,
			  const int DEPTH,
			  const bool EXT,
			  const bool COUNT,
			  const bool GREEN,
                          const bool BNG,
                          const bool PURPLE1,
			  const bool NEATO,
			  const bool REL,
                          const bool LINENO,
			  const bool SHOW_INV,
                          const bool SHOW_CN,
                          const bool SHOW_ALIGN,
			  const double FONTSIZE,
			  const double SCALE,
                          const bool ALTREF,
                          const bool COV2) {
    state.SEEDS = SEEDS;
    state.DEPTH = DEPTH;
    state.EXT = EXT;
    state.COUNT = COUNT;
    state.GREEN = GREEN;
    state.BNG = BNG;
    state.SHOW_CN = SHOW_CN;
    state.SHOW_ALIGN = SHOW_ALIGN;
    state.PURPLE1 = PURPLE1;
    state.NEATO = NEATO;
    state.REL = REL;
    state.LINENO = LINENO;
    state.SHOW_INV = SHOW_INV;
    state.FONTSIZE = FONTSIZE;
    state.SCALE = SCALE;
    state.ALTREF = ALTREF;
    state.COV2 = COV2;
}

void ReportEdges( const int nedges, ostream& out, const Bool ambflag )
{    out << "found " << nedges << " edges" << endl;
     if ( nedges == 0 && !ambflag )
     {    out << "Huh, no edges.  If you're looking for an interval on "
               << "the reference sequence,\n"
               << "you might want to try using a larger interval." << endl;    }    }

void NhoodInfoEngine::RunAsClient(const Bool INTERACTIVE, const int DEPTH, 
			     const String& OUT, const bool PNG, const bool PDF, 
			     const int RANDOM_SEED,
			     const String& SEEDS_MINUS, const Bool LABEL_CONTIGS ){

     while(1)
     {    if (INTERACTIVE)
	  {    cout << "? ";
	       flush(cout);
	       String line;
	       getline( cin, line );
	       if ( cin.fail() )  break;
	       if ( !state.SetState( line, cout, False, hb, *this ) ) continue;    }

	  // Find seeds and build neighborhood.

	  vec<int> seeds;
	  const int max_edges = 5000;
	  ostringstream tout;
          Bool ambflag = False;
	  if ( !DefineSeeds( hb, inv, kmers_plus, lines, tol, 
               state.ALTREF ? genome_names_alt : genome_names,
               ambint, ambflag, state.ALTREF ? hits_alt : hits,
               state, RANDOM_SEED, SEEDS_MINUS, seeds, max_edges, tout ) )
	  {    cout << tout.str( );
               if ( !INTERACTIVE ) break;
	       continue;    }
	  cout << tout.str( );
	  vec<Bool> invisible;
	  BuildNhood( hb, seeds, state.DEPTH, invisible, max_edges );
	  int nedges = 0;
	  for ( int e = 0; e < hb.E( ); e++ )
	       if ( !invisible[e] ) nedges++;
	  if ( nedges > max_edges )
	  {    cout << "found > " << max_edges << " edges" << endl;
	       cout << "That's too many edges to view." << endl;
               if ( !INTERACTIVE ) break;
	       continue;    }
          ReportEdges( nedges, cout, ambflag );

	  // Make output.

	  vec<int> used;
	  for ( int e = 0; e < hb.E( ); e++ )
	       if ( !invisible[e] ) used.push_back(e);
          state.SetUsed(used);
	  MakeFasta( hb, used, OUT + ".fasta" );
	  MakeDot( hb, inv, lines, tol, llens, cov, covs, 
               state.ALTREF ? hits_alt : hits,
               state.ALTREF ? genome_names_alt : genome_names, used, state, seeds, 
               invisible, subsam_names, count, LABEL_CONTIGS, OUT + ".dot" );
	  if (PNG) SystemSucceed( "dot -Tpng -o " + OUT + ".png " + OUT + ".dot" );
	  if (PDF) SystemSucceed( "dot -Tpdf -o " + OUT + ".pdf " + OUT + ".dot" );
	  if ( !INTERACTIVE ) break;    }

}


void NhoodInfoEngine::RunAsServer(const String& SERVER_DIR, 
     const Bool LABEL_CONTIGS )
{
	const int nthreads = omp_get_max_threads( );
	ForceAssert( IsDirectory(SERVER_DIR) );
	for ( int t = 0; t < nthreads; t++ ) 
	    {    String dir = SERVER_DIR + "/" + ToString(t+1);
		if ( !IsDirectory(dir) ) Mkdir777(dir);    }
	Ofstream(PID, SERVER_DIR +"/pid");
	PID << getpid() << endl;
	PID.close();
#pragma omp parallel for
	for ( int t = 0; t < nthreads; t++ )
	    {    String dir = SERVER_DIR + "/" + ToString(t+1);
	while(1)
	    {    vec<String> all = AllFiles(dir);

	// Delete old files and sort by age.

	vec<Bool> to_delete( all.size( ), False );
	vec<Bool> age;
	for ( int i = 0; i < all.isize( ); i++ )
	    {    String fn = dir + "/" + all[i];
	double agei = AgeInMinutes(fn);
	if ( agei > 5 )
	    {    Remove(fn);
	to_delete[i] = True;    }
	else age.push_back(agei);    }
	EraseIf( all, to_delete );
	ReverseSortSync( age, all );

	// Find first unsatisfied query.

	String req, O;
	for ( int i = 0; i < all.isize( ); i++ )
	    {    if ( !all[i].Contains( ".req", -1 ) ) continue;
	req = FirstLineOfFile( dir + "/" + all[i] );
	Remove( dir + "/" + all[i] );
	O = dir + "/" + all[i].RevBefore( ".req" );
	break;    }
	if ( req == "" )
	    {    usleep(100000);
	continue;    }

	// Process query.

	class nhood_info_state state;
	state.Initialize( );
	state.SCALE = 2;
#pragma omp critical
	{    cout << Date( ) << ": #" << t+1 << " processing "
		  << req << endl;    }
	Ofstream( tout, O + ".txt2");
	if ( !state.SetState( req, tout, True, hb, *this ) ) 
	    {
#pragma omp critical
		{    cout << Date( ) << ": #" << t+1 << " done"
			  << endl;    }
		Rename(O+".txt2",O+".txt");
		continue;
	    }
	    vec<int> seeds;
	    const int max_edges = 2000;
            Bool ambflag = False;
	    if ( !DefineSeeds( hb, inv, kmers_plus, lines, tol, 
                 state.ALTREF ? genome_names_alt : genome_names,
                 ambint, ambflag, state.ALTREF ? hits_alt : hits,
                 state, -1, "", seeds, max_edges, tout ) )
		{    cout << "Could not parse seeds." << endl;
#pragma omp critical
		{    cout << Date( ) << ": #" << t+1 << " done"
			  << endl;    }
		Rename(O+".txt2",O+".txt");
		continue;
	    }
	    vec<Bool> invisible;
	    BuildNhood( hb, seeds, state.DEPTH, invisible, max_edges );
	    int nedges = 0;
	    for ( int e = 0; e < hb.E( ); e++ )
		if ( !invisible[e] ) nedges++;
	    if ( nedges > max_edges )
		{    tout << "found > " << max_edges << " edges" << endl;
		    tout << "That's too many edges to view in demo." << endl;
#pragma omp critical
		    {    cout << Date( ) << ": #" << t+1 << " done"
			      << endl;    }
		    Rename(O+".txt2",O+".txt");
		    continue;    }
                ReportEdges( nedges, tout, ambflag );
		vec<int> used;
		for ( int e = 0; e < hb.E( ); e++ )
		    if ( !invisible[e] ) used.push_back(e);
                state.SetUsed(used);
		MakeDot( hb, inv, lines, tol, llens, cov, covs, 
                     state.ALTREF ? hits_alt : hits,
                     state.ALTREF ? genome_names_alt : genome_names,
                     used, state, seeds, invisible, subsam_names,
 	             count, LABEL_CONTIGS, O + ".dot" );
		int status = System( 
		    "dot -Tsvg -o " + O + ".svg.1 " + O + ".dot");
		if ( status != 0 )
                    {    tout << "There was a problem generating the image." << endl
                              << "We have never seen this happen before." << endl
                              << "Maybe try again.  We will investigate." << endl;
		    cout << Date( ) << ": #" << t+1 
			 << " image generation failed" << endl;
		    Rename(O+".txt2",O+".txt");
		    continue;    }
		Remove( O + ".dot" );
		Rename( O + ".svg.1", O + ".svg" );
#pragma omp critical
		{    cout << Date( ) << ": #" << t+1 << " done"
			  << endl;    }
		Rename(O+".txt2",O+".txt");
		continue;
		
	    }    // while (1)
	
	    }     // for nthreads
	
}
