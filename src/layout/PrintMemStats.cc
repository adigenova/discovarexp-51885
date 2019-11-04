// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#include "layout/PrintMemStats.h"
#include <sys/resource.h>
#include "system/Assert.h"
#include <iomanip>

long maxrss_ = 8000 ;
long ixrss_ = 20 ;
long idrss_ = 20 ;
long isrss_ = 20 ;
long minflt_ = 20 ;

const double min_ratio = 1.1 ; 

bool PrintMemSt( std::ostream &o, bool force ) {
  struct rusage self_usage ;
  (void) getrusage(  RUSAGE_SELF, &self_usage ) ; 
  bool print_this = force ; 

  if ( (double) self_usage.ru_maxrss > min_ratio * (double) maxrss_ ) print_this = true ; 
#ifdef KEN
  if ( (double) self_usage.ru_minflt > min_ratio * (double) minflt_ ) print_this = true ; 
#endif
  if ( print_this ) {
    o << " max memory usage = " << self_usage.ru_maxrss << "K" ;
    if ( self_usage.ru_maxrss > min_ratio * (double) maxrss_ ) 
      o << " UP " << std::setw(4) << std::setprecision(2) <<
	(int) ( 100.0 * ( (double) self_usage.ru_maxrss - (double) maxrss_ ) / (double) maxrss_ ) << "% " ;
#ifdef KEN
    o  << " page reclaims = " << self_usage.ru_minflt ; 
    if ( self_usage.ru_minflt > min_ratio * (double) minflt_ ) 
      o << " UP " <<
	(int) ( 100.0 * ( (double) self_usage.ru_minflt - (double) minflt_ ) / (double) minflt_ ) << "% " ;
#endif
    maxrss_ = self_usage.ru_maxrss ;
    minflt_ = self_usage.ru_minflt ;
  }
  return print_this ; 
}

