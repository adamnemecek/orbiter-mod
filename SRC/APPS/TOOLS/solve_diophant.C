// solve_diophant.C
// 
// Anton Betten
// May 17, 2015
//
// 

#include "orbiter.h"
#include "discreta.h"


// global data:

INT t0; // the system time when the program started


int main(int argc, char **argv)
{
	INT i;
	t0 = os_ticks();
	INT verbose_level = 0;
	INT f_file = FALSE;	
	const BYTE *fname = NULL;
	//INT f_maxdepth = FALSE;
	//INT maxdepth = 0;
	//INT print_interval = 1000;
	INT f_list_of_cases = FALSE;
	const BYTE *fname_list_of_cases = NULL;
	const BYTE *fname_template = NULL;
	INT f_prefix = FALSE;
	const BYTE *prefix = NULL;
	INT f_output_file = FALSE;
	const BYTE *output_file = NULL;
	INT f_print = FALSE;
	INT f_draw = FALSE;
	INT f_tree = FALSE;
	INT f_decision_nodes_only = FALSE;
	const BYTE *fname_tree = NULL;
	INT f_coordinates = FALSE;
	INT xmax_in = ONE_MILLION;
	INT ymax_in = ONE_MILLION;
	INT xmax_out = ONE_MILLION;
	INT ymax_out = ONE_MILLION;
	INT f_output_solution_raw = FALSE;
	INT f_test = FALSE;
	const BYTE *solution_file = NULL;

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
		else if (strcmp(argv[i], "-tree") == 0) {
			f_tree = TRUE;
			f_decision_nodes_only = FALSE;
			fname_tree = argv[++i];
			cout << "-tree " << fname_tree << endl;
			}
		else if (strcmp(argv[i], "-tree_decision_nodes_only") == 0) {
			f_tree = TRUE;
			f_decision_nodes_only = TRUE;
			fname_tree = argv[++i];
			cout << "-tree_decision_nodes_only " << fname_tree << endl;
			}
		else if (strcmp(argv[i], "-list_of_cases") == 0) {
			f_list_of_cases = TRUE;
			fname_list_of_cases = argv[++i];
			fname_template = argv[++i];
			cout << "-list_of_cases " << fname_list_of_cases << " " << fname_template << endl;
			}
		else if (strcmp(argv[i], "-prefix") == 0) {
			f_prefix = TRUE;
			prefix = argv[++i];
			cout << "-prefix " << prefix << endl;
			}
		else if (strcmp(argv[i], "-output_file") == 0) {
			f_output_file = TRUE;
			output_file = argv[++i];
			cout << "-output_file " << output_file << endl;
			}
		else if (strcmp(argv[i], "-print") == 0) {
			f_print = TRUE;
			cout << "-print " << endl;
			}
		else if (strcmp(argv[i], "-draw") == 0) {
			f_draw = TRUE;
			cout << "-draw " << endl;
			}
		else if (strcmp(argv[i], "-output_solution_raw") == 0) {
			f_output_solution_raw = TRUE;
			cout << "-output_solution_raw " << endl;
			}
		else if (strcmp(argv[i], "-coordinates") == 0) {
			f_coordinates = TRUE;
			xmax_in = atoi(argv[++i]);
			ymax_in = atoi(argv[++i]);
			xmax_out = atoi(argv[++i]);
			ymax_out = atoi(argv[++i]);
			cout << "-coordinates " << xmax_in << " " << ymax_in << " " << xmax_out << " " << ymax_out << endl;
			}
		else if (strcmp(argv[i], "-test") == 0) {
			f_test = TRUE;
			solution_file = argv[++i];
			cout << "-test " << solution_file << endl;
			}
		
		}


	if (f_file) {

		//INT search_steps, decision_steps, nb_sol, dt;

		diophant *Dio;

		Dio = new diophant;

		Dio->read_compact_format(fname, verbose_level);

		Dio->append_equation();

		INT m, j;
	
		m = Dio->m - 1;
		for (j = 0; j < Dio->n; j++) {
			Dio->Aij(m, j) = 1;
			}
		Dio->type[m] = t_EQ;
		Dio->RHS[m] = Dio->sum;

		if (f_print) {
			BYTE fname_base[1000];

			sprintf(fname_base, "%s", fname);
			replace_extension_with(fname_base, "_print.tex");
			{
				ofstream fp(fname_base);
				Dio->latex_it(fp);
			}
			cout << "Written file " << fname_base << " of size " << file_size(fname_base) << endl;
			}
		
		if (f_draw) {
			BYTE fname_base[1000];

			sprintf(fname_base, "%s", fname);
			replace_extension_with(fname_base, "_drawing");		
			Dio->draw_it(fname_base, xmax_in, ymax_in, xmax_out, ymax_out);
			}

		if (f_test) {
			cout << "testing solutions from file " << solution_file << endl;
			Dio->test_solution_file(solution_file, verbose_level);
			}

		else {
			Dio->solve_all_DLX_with_RHS(f_tree, fname_tree, verbose_level - 2);
			cout << "Found " << Dio->_resultanz << " solutions with " << Dio->nb_steps_betten << " backtrack steps" << endl;
			}

		
		delete Dio;

		}
	else {
		cout << "Please use options -file or -list_of_cases" << endl;
		exit(1);
		}

	cout << "solve_diophant.out is done" << endl;
	the_end(t0);
	//the_end_quietly(t0);

}

