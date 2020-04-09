#pragma once
#if !defined(_LATTICE_MACROS)
#define _LATTICE_MACROS

#define LASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::terminate(); \
        } \
    } while (false)



#endif ///_LATTICE_MACROS