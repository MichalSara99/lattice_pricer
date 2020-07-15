#if !defined(_LATTICE_DISABLE_WARNINGS)
#define _LATTICE_DISABLE_WARNINGS
#if defined(_MSC_VER)
	#define DISABLE_WARNING_PUSH           __pragma(warning( push ))
	#define DISABLE_WARNING_POP            __pragma(warning( pop )) 
	#define DISABLE_WARNING(warningNumber) __pragma(warning( disable : warningNumber ))

	#define DISABLE_WARNING_UNREFERENCED_FORMAL_PARAMETER    DISABLE_WARNING(4100)
	#define DISABLE_WARNING_UNREFERENCED_FUNCTION            DISABLE_WARNING(4505)
#endif 
#endif ///_LATTICE_DISABLE_WARNINGS