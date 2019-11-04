// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef PRINTMEMSTATS
#define PRINTMEMSTATS
//
//  Usage:
//    printmemstats MaxMS               //  Track maximum memory usage statistics
//    MaxMS.PrintMemStats()             //  Print out Memory Stats if they are higher than last time
//    MaxMS.PrintMemStats( ostream &o ) //  Print out Memory Stats to ostream &o
//

#include <sys/resource.h>
#include <iostream>
#if ! defined( __GNUC__ ) && __GNUC__ > 2
#include <ostream>
#endif

bool PrintMemSt( std::ostream &o, bool force ) ; 
#ifdef __GNUC__
#define ForcePrintMem \
     PrintMemSt( cout, true );  \
     cout << " at " <<__PRETTY_FUNCTION__ \
               << ", line " << __LINE__ \
               << " of file " << __FILE__ << " (FORCED) \n";
#define PrintMem \
     if ( PrintMemSt( cout, false )  ) \
     {    cout << " at " <<__PRETTY_FUNCTION__ \
               << ", line " << __LINE__ \
               << " of file " << __FILE__ << "\n"; }
#else
#define ForcePrintMem \
     PrintMemSt( cout, true );  \
     cout << " at line " << __LINE__ \
               << " of file " << __FILE__ << " (FORCED) \n";
#define PrintMem \
     if ( PrintMemSt( cout, false )  ) \
     {    cout << " at line " << __LINE__ \
               << " of file " << __FILE__ << "\n"; }
#endif

extern long maxrss_;
extern long ixrss_;
extern long idrss_;
extern long isrss_;
extern long minflt_;

//  inline void PrintMem( ostream &o = cout ) { if ( PrintMemSt( o ) ) o << " at line "<<
//  								     __LINE__ << " in file " <<
//  								     __FILE__ << endl ; }
#endif

