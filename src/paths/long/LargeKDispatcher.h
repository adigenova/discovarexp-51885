///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * LargeKDispatcher.h
 *
 *  Created on: May 14, 2013
 *      Author: tsharpe
 */

#ifndef PATHS_LONG_LARGEKDISPATCHER_H_
#define PATHS_LONG_LARGEKDISPATCHER_H_

#include "system/System.h"

class BigK
{
    static int constexpr gK[] =
    {20, 24, 40, 48, 60, 72, 80, 84, 88, 100, 128, 144, 172, 200, 320, 368, 400, 
     500, 544, 640, 720, 1000, 1200, 1600, 2000, 4000, 10000};
    // Case 1:  DANGER, WILL ROBINSON.  If you add a value to the array, you
    //          must also add another case statement at the end of the switch.
    //          If you forget, you won't get a dispatch on the last K value,
    //          you'll get a fatal err at runtime.  This would be bad.
    //          Keep the array size and the switch size in sync!!!
    // Case 2:  If you remove a value from the array, you'll have to remove
    //          the final case statement.  But if you forget, the compiler will
    //          warn you because you'll be referring to a non-existent array
    //          element.  So you'll fix it.
    // Case 3:  If you modify a value in the array, everything is fine.
    //          Nothing else changes.

public:

    // iterator pair over allowable values
    static int const* begin() { return gK; }
    static int const* end() { return gK+sizeof(gK)/sizeof(gK[0]); }

    // requirements:  a functor templated on K.  it may have any arg list,
    // and any return type, including void.
    // I.e.,
    // template <int K> struct FunctorT
    // { return_t operator()( arg1_t, arg2_t, etc ); };
    //
    // what it does:  instantiate and call the functor for the value of K (which
    // must be one of the allowable values) that you have passed as an arg.
    //
    // how you invoke it:
    // return_t ret = BigK::dispatch<FunctorT>(k, arg1, arg2, etc);
    //
    template <template <int K> class C, typename... Args>
    static auto dispatch( int k, Args&&... args ) -> decltype( C<40>()(args...) )
    {
        switch (k)
        {
        case gK[0]: return C<gK[0]>()(std::forward<Args>(args)...);
        case gK[1]: return C<gK[1]>()(std::forward<Args>(args)...);
        case gK[2]: return C<gK[2]>()(std::forward<Args>(args)...);
        case gK[3]: return C<gK[3]>()(std::forward<Args>(args)...);
        case gK[4]: return C<gK[4]>()(std::forward<Args>(args)...);
        case gK[5]: return C<gK[5]>()(std::forward<Args>(args)...);
        case gK[6]: return C<gK[6]>()(std::forward<Args>(args)...);
        case gK[7]: return C<gK[7]>()(std::forward<Args>(args)...);
        case gK[8]: return C<gK[8]>()(std::forward<Args>(args)...);
        case gK[9]: return C<gK[9]>()(std::forward<Args>(args)...);
        case gK[10]: return C<gK[10]>()(std::forward<Args>(args)...);
        case gK[11]: return C<gK[11]>()(std::forward<Args>(args)...);
        case gK[12]: return C<gK[12]>()(std::forward<Args>(args)...);
        case gK[13]: return C<gK[13]>()(std::forward<Args>(args)...);
        case gK[14]: return C<gK[14]>()(std::forward<Args>(args)...);
        case gK[15]: return C<gK[15]>()(std::forward<Args>(args)...);
        case gK[16]: return C<gK[16]>()(std::forward<Args>(args)...);
        case gK[17]: return C<gK[17]>()(std::forward<Args>(args)...);
        case gK[18]: return C<gK[18]>()(std::forward<Args>(args)...);
        case gK[19]: return C<gK[19]>()(std::forward<Args>(args)...);
        case gK[20]: return C<gK[20]>()(std::forward<Args>(args)...);
        case gK[21]: return C<gK[21]>()(std::forward<Args>(args)...);
        case gK[22]: return C<gK[22]>()(std::forward<Args>(args)...);
        case gK[23]: return C<gK[23]>()(std::forward<Args>(args)...);
        case gK[24]: return C<gK[24]>()(std::forward<Args>(args)...);
        case gK[25]: return C<gK[25]>()(std::forward<Args>(args)...);
        case gK[26]: return C<gK[26]>()(std::forward<Args>(args)...);
        default: FatalErr("Illegal value " << k << " for K in BigK dispatcher.");
        }
    }
};

#endif /* PATHS_LONG_LARGEKDISPATCHER_H_ */
