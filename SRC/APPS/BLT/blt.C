// blt.C
// 
// Anton Betten
// 8/13/2006
//
// This is the main file for the classification of partial BLT-sets of a given size.
// It relies on the class blt_generator
// as well as the algorithm LIB/SNAKES_AND_LADDERS
// and the clique finder in LIB/GALOIS/clique_finder.C
//
//

#include "orbiter.h"
#include "discreta.h"

#include "blt.h"

#define MAX_FILES 1000

// global data:

INT t0; // the system time when the program started


int main(int argc, const char **argv)
{
	t0 = os_ticks();
	INT verbose_level = 0;
	INT f_q = FALSE;
	INT q = 0;
	INT f_poly = FALSE;
	const BYTE *poly = NULL;
	INT f_starter = FALSE;
	INT f_create_system = FALSE;
	INT f_create_graphs = FALSE;
	INT f_lexorder_test = FALSE;
	INT /*extend_level,*/ extend_orbit, extend_level_candidates;
	INT extend_orbit_r, extend_orbit_m;
	INT f_output_prefix = FALSE;
	const BYTE *output_prefix = "";
	INT i;
	INT f_Law71 = FALSE;
	//INT f_find_cliques = FALSE;
	INT f_build_db = FALSE;
	//INT level = 0;

	INT f_starter_depth = FALSE;
	INT starter_depth = 0;

#if 0
	INT f_read_solution_files = FALSE;
	const BYTE *fname[MAX_FILES];
	INT nb_files = 0;
#endif

	INT f_read_solutions = FALSE;
	INT f_read_solutions_after_split = FALSE;
	INT read_solutions_split_m = 0;
	
	INT f_read_statistics_after_split = FALSE;
	INT read_statistics_split_m = 0;

	INT f_compute_orbits = FALSE;
	INT f_isomorph_testing = FALSE;
	INT f_classification_graph = FALSE;
	INT f_event_file = FALSE; // -e <event file> option
	const BYTE *event_file_name;
	INT print_mod = 500;
	INT f_report = FALSE;
	INT f_subset_orbits = FALSE;
	INT f_eliminate_graphs_if_possible = FALSE;
	INT f_down_orbits = FALSE;
	INT f_draw_poset = FALSE;

	INT f_starter_directory_name = FALSE;
	const BYTE *starter_directory_name = "";
	INT f_has_starter_prefix = FALSE;
	const BYTE *starter_prefix = "";
	INT f_solution_directory_name = FALSE;
	const BYTE *solution_directory_name = "./";


	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-v") == 0) {
			verbose_level = atoi(argv[++i]);
			cout << "-v " << verbose_level << endl;
			}
		else if (strcmp(argv[i], "-q") == 0) {
			f_q = TRUE;
			q = atoi(argv[++i]);
			cout << "-q " << q << endl;
			}
		else if (strcmp(argv[i], "-starter") == 0) {
			f_starter = TRUE;
			cout << "-starter " << endl;
			}
		else if (strcmp(argv[i], "-starter_depth") == 0) {
			f_starter_depth = TRUE;
			starter_depth = atoi(argv[++i]);
			cout << "-starter_depth " << starter_depth << endl;
			}
		else if (strcmp(argv[i], "-poly") == 0) {
			f_poly = TRUE;
			poly = argv[++i];
			cout << "-poly " << poly << endl;
			}
		else if (strcmp(argv[i], "-create_system") == 0) {
			f_create_system = TRUE;
			//extend_level = atoi(argv[++i]);
			extend_orbit = atoi(argv[++i]);
			extend_level_candidates = atoi(argv[++i]);
			cout << "-create_system " << " " << extend_orbit << " " << extend_level_candidates << endl;
			}
		else if (strcmp(argv[i], "-create_graphs") == 0) {
			f_create_graphs = TRUE;
			//extend_level = atoi(argv[++i]);
			extend_orbit_r = atoi(argv[++i]);
			extend_orbit_m = atoi(argv[++i]);
			extend_level_candidates = atoi(argv[++i]);
			cout << "-create_graphs " << " " << extend_orbit_r << " " << extend_orbit_m << " " << extend_level_candidates << endl;
			}
		else if (strcmp(argv[i], "-lex") == 0) {
			f_lexorder_test = TRUE;
			cout << "-lex" << endl;
			}
		else if (strcmp(argv[i], "-output_prefix") == 0) {
			f_output_prefix = TRUE;
			output_prefix = argv[++i];
			cout << "-output_prefix " << output_prefix << endl;
			}
		else if (strcmp(argv[i], "-Law71") == 0) {
			f_Law71 = TRUE;
			cout << "-Law71" << endl;
			}

#if 0
		else if (strcmp(argv[i], "-find_cliques") == 0) {
			f_find_cliques = TRUE;
			extend_level = atoi(argv[++i]);
			extend_orbit = atoi(argv[++i]);
			cout << "-find_cliques " << extend_level << " " << extend_orbit << endl;
			}
#endif

		else if (strcmp(argv[i], "-build_db") == 0) {
			f_build_db = TRUE;
			//level = atoi(argv[++i]);
			cout << "-build_db " << endl;
			}
#if 0
		else if (strcmp(argv[i], "-read_solution_files") == 0) {
			f_read_solution_files = TRUE;
			level = atoi(argv[++i]);
			i++;
			while (i < argc) {
				fname[nb_files] = argv[i];
				if (strcmp(fname[nb_files], "-1") == 0) {
					break;
					}
				nb_files++;
				i++;
				}
			cout << "-read_solution_files ";
			INT j;
			for (j = 0; j < nb_files; j++) {
				cout << fname[j] << " ";
				}
			cout << endl;
			}
#endif
		else if (strcmp(argv[i], "-read_solutions") == 0) {
			f_read_solutions = TRUE;
			cout << "-read_solutions " << endl;
			}
		else if (strcmp(argv[i], "-read_solutions_after_split") == 0) {
			f_read_solutions_after_split = TRUE;
			read_solutions_split_m = atoi(argv[++i]);
			cout << "-read_solutions_after_split " << read_solutions_split_m << endl;
			}
		else if (strcmp(argv[i], "-read_statistics_after_split") == 0) {
			f_read_statistics_after_split = TRUE;
			read_statistics_split_m = atoi(argv[++i]);
			cout << "-read_statistics_after_split " << read_statistics_split_m << endl;
			}

		else if (strcmp(argv[i], "-compute_orbits") == 0) {
			f_compute_orbits = TRUE;
			//level = atoi(argv[++i]);
			cout << "-compute_orbits " << endl;
			}
		else if (strcmp(argv[i], "-isomorph_testing") == 0) {
			f_isomorph_testing = TRUE;
			//level = atoi(argv[++i]);
			cout << "-isomorph_testing " << endl;
			}
		else if (strcmp(argv[i], "-classification_graph") == 0) {
			f_classification_graph = TRUE;
			//level = atoi(argv[++i]);
			cout << "-make_classification_graph " << endl;
			}
		else if (strcmp(argv[i], "-e") == 0) {
			i++;
			f_event_file = TRUE;
			event_file_name = argv[i];
			cout << "-e " << event_file_name << endl;
			}
		else if (strcmp(argv[i], "-print_interval") == 0) {
			print_mod = atoi(argv[++i]);
			cout << "-print_interval " << print_mod << endl;
			}
		else if (strcmp(argv[i], "-report") == 0) {
			f_report = TRUE;
			//level = atoi(argv[++i]);
			cout << "-report " << endl;
			}
		else if (strcmp(argv[i], "-subset_orbits") == 0) {
			f_subset_orbits = TRUE;
			//level = atoi(argv[++i]);
			cout << "-subset_orbits " << endl;
			}
		else if (strcmp(argv[i], "-eliminate_early") == 0) {
			f_eliminate_graphs_if_possible = TRUE;
			cout << "-eliminate_early " << endl;
			}
		else if (strcmp(argv[i], "-down_orbits") == 0) {
			f_down_orbits = TRUE;
			//level = atoi(argv[++i]);
			cout << "-down_orbits " << endl;
			}
		else if (strcmp(argv[i], "-draw_poset") == 0) {
			f_draw_poset = TRUE;
			cout << "-draw_poset " << endl;
			}
		else if (strcmp(argv[i], "-starter_directory_name") == 0) {
			f_starter_directory_name = TRUE;
			starter_directory_name = argv[++i];
			cout << "-starter_directory_name " << starter_directory_name << endl;
			}
		else if (strcmp(argv[i], "-solution_directory_name") == 0) {
			f_solution_directory_name = TRUE;
			solution_directory_name = argv[++i];
			cout << "-solution_directory_name " << solution_directory_name << endl;
			}
		else if (strcmp(argv[i], "-starter_prefix") == 0) {
			f_has_starter_prefix = TRUE;
			starter_prefix = argv[++i];
			cout << "-starter_prefix " << starter_prefix << endl;
			}

		}

	INT f_v = (verbose_level >= 1);

	if (!f_q) {
		cout << "Please use option -q <q>" << endl;
		exit(1);
		}

	{
	blt_generator Gen;
	INT schreier_depth = 5;
	INT f_debug = FALSE;
	INT f_implicit_fusion = FALSE;
		
	Gen.init_basic(q, starter_directory_name, starter_prefix, 
		argc, argv, verbose_level);
	
	
	Gen.init_group(verbose_level);
	
	Gen.init2(verbose_level);
	
	INT f_use_invariant_subset_if_available = TRUE;

	if (Gen.f_override_schreier_depth) {
		schreier_depth = Gen.override_schreier_depth;
		}
	
	if (f_v) {
		cout << "init finished, calling main, schreier_depth = " << schreier_depth << endl;
		}


	if (f_starter) {

		INT depth;
		INT f_embedded = TRUE;

		depth = Gen.gen->main(t0, schreier_depth, 
			f_use_invariant_subset_if_available, 
			f_implicit_fusion, 
			f_debug, 
			Gen.gen->verbose_level);
		cout << "Gen.gen->main returns depth=" << depth << endl;
		//Gen.gen->print_data_structure_tex(depth, Gen.gen->verbose_level);
		if (f_draw_poset) {
			Gen.gen->draw_poset(Gen.gen->fname_base, depth, 0 /* data1 */, f_embedded, Gen.gen->verbose_level);
			}

#if 0
		if (q == 23) {
			INT set_in[3] = {1, 2728, 7448 };
			INT set_out[3];
			INT set_size = 3;
			INT level = 3;
			INT *Elt;
			INT f_implicit_fusion = TRUE;
			INT case_nb;

			Elt = NEW_INT(Gen.gen->A->elt_size_in_INT);
			case_nb = Gen.gen->trace_set(set_in, set_size, level, 
				set_out, Elt, 
				f_implicit_fusion, verbose_level);

			cout << "case_nb=" << case_nb << endl;
			cout << "canonical set=";
			INT_vec_print(cout, set_out, set_size);
			cout << endl;
			INT v5[5];
			for (i = 0; i < set_size; i++) {
				Gen.O->unrank_point(v5, 1, set_out[i], 0 /* verbose_level */);
				cout << i << " : " << set_out[i] << " : ";
				INT_vec_print(cout, v5, 5);
				cout << endl;
				}
			cout << "transporter:" << endl;
			Gen.gen->A->element_print(Elt, cout);
			}
#endif
		}
	else if (f_create_system) {
		BYTE fname[1000];
		//INT orbit;
		INT nb_orbits;
		INT width;


		if (!f_starter_depth) {
			cout << "Please use option -starter_depth <starter_depth>" << endl;
			exit(1);
			}

		sprintf(fname, "%s%s_lvl_%ld", starter_directory_name, starter_prefix, starter_depth);
		nb_orbits = count_number_of_orbits_in_file(fname, 0);
		cout << "There are " << nb_orbits << " starters" << endl;
		if (nb_orbits < 0) {
			cout << "Something is wrong" << endl;
			exit(1);
			}

		width = (INT)(log10(nb_orbits)) + 1;
		Gen.create_system(starter_depth, 
			extend_orbit, 
			extend_level_candidates, 
			output_prefix, width, verbose_level);
		}
	else if (f_create_graphs) {

		if (!f_starter_depth) {
			cout << "Please use option -starter_depth <starter_depth>" << endl;
			exit(1);
			}
		Gen.create_graphs(starter_depth, 
			extend_orbit_r, extend_orbit_m, 
			extend_level_candidates, 
			output_prefix, 
			f_lexorder_test, f_eliminate_graphs_if_possible, 
			verbose_level);
		}
	else if (f_Law71) {
		Gen.Law_71(verbose_level);
		}
#if 0
	else if (f_find_cliques) {
		Gen.find_cliques(extend_level, extend_orbit, verbose_level);
		}
#endif
	else if (f_build_db) {

		cout << "build_db" << endl;
		
		system("mkdir ISO");

		BYTE prefix[1000];
		sprintf(prefix, "%s%s", starter_directory_name, starter_prefix);

		cout << "prefix=" << prefix << endl;


		if (!f_starter_depth) {
			cout << "please use option -starter_depth <starter_depth>" << endl;
			exit(1);
			}

		cout << "starter_depth=" << starter_depth << endl;

		isomorph_build_db(Gen.A, Gen.A, Gen.gen, 
			Gen.target_size, prefix, (BYTE *) "ISO/", 
			starter_depth, verbose_level);
		}
#if 0
	else if (f_read_solution_files) {

		if (!f_starter_depth) {
			cout << "Please use option -starter_depth <starter_depth>" << endl;
			exit(1);
			}

		isomorph_read_solution_files_from_clique_finder(Gen.A, Gen.A, Gen.gen, 
			Gen.target_size, Gen.gen->fname_base, (BYTE *) "ISO/", 
			starter_depth, 
			fname, nb_files, verbose_level);
		}
#endif
	else if (f_read_solutions) {


		BYTE prefix[1000];
		sprintf(prefix, "%s%s", starter_directory_name, starter_prefix);

		cout << "prefix=" << prefix << endl;

		if (!f_starter_depth) {
			cout << "please use option -starter_depth <starter_depth>" << endl;
			exit(1);
			}

		cout << "starter_depth=" << starter_depth << endl;


		const BYTE *fname[1];
		BYTE fname1[1000];
		INT nb_files = 1;
		
		sprintf(fname1, "%ssolutions_%ld_%ld_0_1.txt", solution_directory_name, q, starter_depth);
		fname[0] = fname1;


		isomorph_read_solution_files_from_clique_finder(Gen.A, Gen.A, Gen.gen, 
			Gen.target_size, prefix, (const BYTE *) "ISO/", starter_depth, 
			fname, nb_files, verbose_level);
#if 0
		isomorph_read_solution_files(Gen.A, Gen.A, Gen.gen, 
			Gen.target_size, prefix, (const BYTE *)"ISO/", starter_depth, 
			(const BYTE **) fname, nb_files, verbose_level);
#endif
		}


	else if (f_read_solutions_after_split) {


		BYTE prefix[1000];
		sprintf(prefix, "%s%s", starter_directory_name, starter_prefix);

		cout << "prefix=" << prefix << endl;

		if (!f_starter_depth) {
			cout << "please use option -starter_depth <starter_depth>" << endl;
			exit(1);
			}

		cout << "starter_depth=" << starter_depth << endl;



		BYTE **fname;
		BYTE fname1[1000];
		INT nb_files = 0;

		nb_files = read_solutions_split_m;
		fname = NEW_PBYTE(nb_files);
		for (i = 0; i < read_solutions_split_m; i++) {
			sprintf(fname1, "%ssolutions_%ld_%ld_%ld_%ld.txt", solution_directory_name, q, starter_depth, i, read_solutions_split_m);
			fname[i] = NEW_BYTE(strlen(fname1) + 1);
			strcpy(fname[i], fname1);
			}
		cout << "Reading the following " << nb_files << " files:" << endl;
		for (i = 0; i < nb_files; i++) {
			cout << i << " : " << fname[i] << endl;
			}


		
		isomorph_read_solution_files_from_clique_finder(Gen.A, Gen.A, Gen.gen, 
			Gen.target_size, prefix, (const BYTE *) "ISO/", starter_depth, 
			(const BYTE **) fname, nb_files, verbose_level);
#if 0
		isomorph_read_solution_files(Gen.A, Gen.A, Gen.gen, 
			Gen.target_size, prefix, (const BYTE *)"ISO/", starter_depth, 
			(const BYTE **) fname, nb_files, verbose_level);
#endif
		}


	else if (f_read_statistics_after_split) {
		BYTE **fname;
		BYTE fname1[1000];
		INT nb_files = 0;

		cout << "f_read_statistics_after_split" << endl;
		cout << "starter_directory_name=" << starter_directory_name << endl;
		cout << "starter_prefix=" << starter_prefix << endl;
		cout << "solution_directory_name=" << solution_directory_name << endl;
		


		if (!f_starter_depth) {
			cout << "please use option -starter_depth <starter_depth>" << endl;
			exit(1);
			}

		cout << "starter_depth=" << starter_depth << endl;

		nb_files = read_statistics_split_m;
		fname = NEW_PBYTE(nb_files);
		for (i = 0; i < read_statistics_split_m; i++) {
			sprintf(fname1, "%ssolutions_%ld_%ld_%ld_%ld_stats.txt", solution_directory_name, q, starter_depth, i, read_solutions_split_m);
			fname[i] = NEW_BYTE(strlen(fname1) + 1);
			strcpy(fname[i], fname1);
			}
		cout << "Reading the following " << nb_files << " files:" << endl;
		for (i = 0; i < nb_files; i++) {
			cout << i << " : " << fname[i] << endl;
			}

		BYTE prefix[1000];

		sprintf(prefix, "%s%s", starter_directory_name, starter_prefix);
		
		isomorph_read_statistic_files(Gen.A, Gen.A, Gen.gen, 
			Gen.target_size, prefix, (BYTE *)"ISO/", starter_depth, 
			(const BYTE **) fname, nb_files, verbose_level);

		for (i = 0; i < nb_files; i++) {
			FREE_BYTE(fname[i]);
			}
		FREE_PBYTE(fname);
		}

	else if (f_compute_orbits) {

		BYTE prefix[1000];
		sprintf(prefix, "%s%s", starter_directory_name, starter_prefix);

		cout << "prefix=" << prefix << endl;

		if (!f_starter_depth) {
			cout << "please use option -starter_depth <starter_depth>" << endl;
			exit(1);
			}

		cout << "starter_depth=" << starter_depth << endl;

		isomorph_compute_orbits(Gen.A, Gen.A, Gen.gen, 
			Gen.target_size, prefix, (BYTE *) "ISO/", 
			starter_depth, verbose_level);
		}
	else if (f_isomorph_testing) {

		BYTE prefix[1000];
		sprintf(prefix, "%s%s", starter_directory_name, starter_prefix);

		cout << "prefix=" << prefix << endl;

		if (!f_starter_depth) {
			cout << "please use option -starter_depth <starter_depth>" << endl;
			exit(1);
			}

		cout << "starter_depth=" << starter_depth << endl;

		isomorph_testing(Gen.A, Gen.A, Gen.gen, 
			Gen.target_size, prefix, (BYTE *) "ISO/", 
			starter_depth, 
			f_event_file, event_file_name, print_mod, 
			verbose_level);
		}
	else if (f_classification_graph) {

		BYTE prefix[1000];
		sprintf(prefix, "%s%s", starter_directory_name, starter_prefix);

		cout << "prefix=" << prefix << endl;

		if (!f_starter_depth) {
			cout << "please use option -starter_depth <starter_depth>" << endl;
			exit(1);
			}

		cout << "starter_depth=" << starter_depth << endl;

		isomorph_classification_graph(Gen.A, Gen.A, Gen.gen, 
			Gen.target_size, 
			prefix, (BYTE *)"ISO/", 
			starter_depth, 
			verbose_level);
		}
	else if (f_report) {

		BYTE prefix[1000];
		sprintf(prefix, "%s%s", starter_directory_name, starter_prefix);

		cout << "prefix=" << prefix << endl;

		if (!f_starter_depth) {
			cout << "please use option -starter_depth <starter_depth>" << endl;
			exit(1);
			}

		isomorph_worker(Gen.A, Gen.A, Gen.gen, 
			Gen.target_size, prefix, (BYTE *) "ISO/", 
			callback_report, &Gen, 
			starter_depth, verbose_level);
		}
	else if (f_subset_orbits) {

		BYTE prefix[1000];
		sprintf(prefix, "%s%s", starter_directory_name, starter_prefix);

		cout << "prefix=" << prefix << endl;

		if (!f_starter_depth) {
			cout << "please use option -starter_depth <starter_depth>" << endl;
			exit(1);
			}

		isomorph_worker(Gen.A, Gen.A, Gen.gen, 
			Gen.target_size, prefix, (BYTE *) "ISO/", 
			callback_subset_orbits, &Gen, 
			starter_depth, verbose_level);
		}
	else if (f_down_orbits) {

		BYTE prefix[1000];
		sprintf(prefix, "%s%s", starter_directory_name, starter_prefix);

		cout << "prefix=" << prefix << endl;

		if (!f_starter_depth) {
			cout << "please use option -starter_depth <starter_depth>" << endl;
			exit(1);
			}

		isomorph_compute_down_orbits(Gen.A, Gen.A, Gen.gen, 
			Gen.target_size, 
			prefix, (BYTE *)"ISO/", 
			&Gen, 
			starter_depth, verbose_level);
		}

	cout << "cleaning up Gen" << endl;
	}


	the_end(t0);
	//the_end_quietly(t0);
}





