// sub.C
//
// Jayant Apte
// November 16, 2015

#include "orbiter.h"
#include "discreta.h"
#include "sub.h"

// global data:

INT t0; // the system time when the program started
const BYTE *version = "sub.C version 5/18/2009";

void print_usage()
{
	cout << "usage: sub.out [options] -n <n> -klist <k_1 k_2...k_i> -q <q> -d <d>" << endl;
	cout << "where options can be:" << endl;

	generator gen;

	gen.usage();

}


int main(int argc, const char **argv)
{
	cout << version << endl;
	t0 = os_ticks();


	{
	sub_generator sg;

	sg.init(argc, argv);

	sg.main();

	//sg.gen->print_tree();

	}
	the_end(t0);
	//the_end_quietly(t0);
}

