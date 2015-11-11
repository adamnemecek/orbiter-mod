// hyperoval.C
// 
// Anton Betten
//
// previous version Dec 6, 2004
// revised June 19, 2006
// revised Aug 17, 2008
//
// Searches for arcs and hyperovals
//
//

#include "orbiter.h"
#include "discreta.h"
#include "hyperoval.h"

// global data:

INT t0; // the system time when the program started


int main(int argc, const char **argv)
{
	t0 = os_ticks();
	
	
	{
	arc_generator Gen;
	
	Gen.init(argc, argv);
	
	Gen.main(argc, argv);
	
	
	}
	//the_end(t0);
	the_end_quietly(t0);
}



