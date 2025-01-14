// all_rainbow_cliques.C
// 
// Anton Betten
// October 28, 2012
//
// 
// previously called all_cliques.C
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
	INT f_maxdepth = FALSE;
	INT maxdepth = 0;
	INT print_interval = 1000;
	INT f_list_of_cases = FALSE;
	const BYTE *fname_list_of_cases = NULL;
	const BYTE *fname_template = NULL;
	INT f_prefix = FALSE;
	const BYTE *prefix = NULL;
	INT f_output_file = FALSE;
	const BYTE *output_file = NULL;
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
		}

	if (f_file) {

		INT search_steps, decision_steps, nb_sol, dt;

		colored_graph_all_cliques(fname, f_output_solution_raw, 
			f_draw, xmax_in, ymax_in, xmax_out, ymax_out, 
			f_output_file, output_file, 
			f_maxdepth, maxdepth, 
			f_tree, f_decision_nodes_only, fname_tree,  
			print_interval, 
			search_steps, decision_steps, nb_sol, dt, 
			verbose_level);
			// in GALOIS/colored_graph.C
		}
	else if (f_list_of_cases) {
		INT *list_of_cases;
		INT nb_cases;
		BYTE fname_sol[1000];
		BYTE fname_stats[1000];
		
		if (f_output_file) {
			sprintf(fname_sol, "%s", output_file);
			sprintf(fname_stats, "%s", output_file);
			replace_extension_with(fname_stats, "_stats.csv");
			}
		else {
			sprintf(fname_sol, "solutions_%s", fname_list_of_cases);
			sprintf(fname_stats, "statistics_%s", fname_list_of_cases);
			replace_extension_with(fname_stats, ".csv");
			}
		read_set_from_file(fname_list_of_cases, list_of_cases, nb_cases, verbose_level);
		cout << "nb_cases=" << nb_cases << endl;

		colored_graph_all_cliques_list_of_cases(list_of_cases, nb_cases, f_output_solution_raw, 
			f_draw, xmax_in, ymax_in, xmax_out, ymax_out, 
			fname_template, 
			fname_sol, fname_stats, 
			f_maxdepth, maxdepth, 
			f_prefix, prefix, 
			print_interval, verbose_level);
		
		FREE_INT(list_of_cases);
		cout << "all_rainbow_cliques.out written file " << fname_sol << " of size " << file_size(fname_sol) << endl;
		cout << "all_rainbow_cliques.out written file " << fname_stats << " of size " << file_size(fname_stats) << endl;
		}
	else {
		cout << "Please use options -file or -list_of_cases" << endl;
		exit(1);
		}

	cout << "all_rainbow_cliques.out is done" << endl;
	the_end(t0);
	//the_end_quietly(t0);

}



