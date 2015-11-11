// rls_extend.C
// 
// Anton Betten
// 1/2/13
//
// 
//
//

#include "orbiter.h"
#include "discreta.h"

#include "regular_ls.h"


// global data:

INT t0; // the system time when the program started



int main(int argc, const char **argv)
{
	INT i;
	INT verbose_level = 0;
	INT f_file = FALSE;
	const BYTE *fname;
	INT f_N = FALSE;
	INT N = 0;
	INT f_K = FALSE;
	INT K = 0;
	INT f_R = FALSE;
	INT R = 0;
	INT f_lambda_reached = FALSE;
	INT f_depth = FALSE;
	INT depth = 0;
	INT f_lexorder_test = FALSE;
	INT f_memory_debug = FALSE;
	INT f_single_case = FALSE;
	INT single_case = 0;

	t0 = os_ticks();

	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-v") == 0) {
			verbose_level = atoi(argv[++i]);
			cout << "-v " << verbose_level << endl;
			}
		else if (strcmp(argv[i], "-file") == 0) {
			f_file = TRUE;
			fname = argv[++i];
			cout << "-file " << fname << endl;
			}
		else if (strcmp(argv[i], "-N") == 0) {
			f_N = TRUE;
			N = atoi(argv[++i]);
			cout << "-N " << N << endl;
			}
		else if (strcmp(argv[i], "-K") == 0) {
			f_K = TRUE;
			K = atoi(argv[++i]);
			cout << "-K " << K << endl;
			}
		else if (strcmp(argv[i], "-R") == 0) {
			f_R = TRUE;
			R = atoi(argv[++i]);
			cout << "-R " << R << endl;
			}
		else if (strcmp(argv[i], "-lambda_reached") == 0) {
			f_lambda_reached = TRUE;
			cout << "-lambda_reached" << endl;
			}
		else if (strcmp(argv[i], "-depth") == 0) {
			f_depth = TRUE;
			depth = atoi(argv[++i]);
			cout << "-depth " << depth << endl;
			}
		else if (strcmp(argv[i], "-lex") == 0) {
			f_lexorder_test = TRUE;
			cout << "-lex" << endl;
			}
		else if (strcmp(argv[i], "-memory_debug") == 0) {
			f_memory_debug = TRUE;
			cout << "-memory_debug" << endl;
			}
		else if (strcmp(argv[i], "-single_case") == 0) {
			f_single_case = TRUE;
			single_case = atoi(argv[++i]);
			cout << "-single_case " << single_case << endl;
			}
		}

	if (!f_file) {
		cout << "Please use option -file <fname>" << endl;
		exit(1);
		}
	if (!f_K) {
		cout << "please specify K using -K <K>" << endl;
		exit(1);
		}
	if (!f_N) {
		cout << "please specify N using -N <N>" << endl;
		exit(1);
		}
	if (!f_R) {
		cout << "please specify R using -R <R>" << endl;
		exit(1);
		}
	if (!f_depth) {
		cout << "please specify depth using -depth <depth>" << endl;
		exit(1);
		}

	if (f_memory_debug) {
		start_memory_debug();
		}
	//INT f_v = (verbose_level >= 1);
	{
	regular_ls_generator Gen;

	Gen.init_basic(argc, argv, 0 /*verbose_level*/);
	
	Gen.init_group(0 /*verbose_level*/);

	Gen.init_action_on_k_subsets(K, 0 /*verbose_level*/);
	Gen.onr = R;

	Gen.extend(fname, 
		f_single_case, single_case, 
		N, K, R, f_lambda_reached, depth, f_lexorder_test, verbose_level);

	}
	if (f_memory_debug) {
		registry_dump_sorted();
		}

	the_end(t0);
	//the_end_quietly(t0);
}



