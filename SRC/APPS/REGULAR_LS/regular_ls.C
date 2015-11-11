// regular_ls.C
// 
// Anton Betten
// 1/1/13
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
	INT f_starter = FALSE;
	INT f_starter_size = FALSE;
	INT starter_size = 0;
	INT f_starter_directory_name = FALSE;
	const BYTE *starter_directory_name = NULL;
	INT f_draw_poset = FALSE;
	INT f_embedded = FALSE;
	INT f_build_db = FALSE;
	INT f_memory_debug = FALSE;
	//INT f_prefix = FALSE;
	//const BYTE *prefix = NULL;
	INT f_lift = FALSE;
	const BYTE *prefix_lift = "./";
	INT f_lex = FALSE;
	INT f_split = FALSE;
	INT split_r = 0, split_m = 0;
	INT f_solve = FALSE;
	INT f_save = FALSE;
	INT f_read_instead = FALSE;
	INT f_draw_system = FALSE;
	const BYTE *fname_system = NULL;
	INT f_write_tree = FALSE;
	const BYTE *fname_tree = NULL;
	INT f_solution_directory_name = FALSE;
	const BYTE *solution_directory_name = "./";
	INT f_read_solutions = FALSE;
	INT f_read_solutions_after_split = FALSE;
	INT read_solutions_split_m = 0;
	INT f_compute_orbits = FALSE;
	INT f_isomorph_testing = FALSE;
	INT f_event_file = FALSE; // -e <event file> option
	const BYTE *event_file_name;
	INT print_mod = 500;
	

	t0 = os_ticks();

	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-v") == 0) {
			verbose_level = atoi(argv[++i]);
			cout << "-v " << verbose_level << endl;
			}
		else if (strcmp(argv[i], "-starter") == 0) {
			f_starter = TRUE;
			cout << "-starter" << endl;
			}
		else if (strcmp(argv[i], "-starter_size") == 0) {
			f_starter_size = TRUE;
			starter_size = atoi(argv[++i]);
			cout << "-starter_size " << starter_size << endl;
			}
		else if (strcmp(argv[i], "-starter_directory_name") == 0) {
			f_starter_directory_name = TRUE;
			starter_directory_name = argv[++i];
			cout << "-starter_directory_name " << starter_directory_name << endl;
			}
		else if (strcmp(argv[i], "-draw_poset") == 0) {
			f_draw_poset = TRUE;
			cout << "-draw_poset " << endl;
			}
		else if (strcmp(argv[i], "-embedded") == 0) {
			f_embedded = TRUE;
			cout << "-embedded " << endl;
			}
		else if (strcmp(argv[i], "-build_db") == 0) {
			f_build_db = TRUE;
			cout << "-build_db" << endl;
			}
		else if (strcmp(argv[i], "-memory_debug") == 0) {
			f_memory_debug = TRUE;
			cout << "-memory_debug" << endl;
			}
		else if (strcmp(argv[i], "-lift") == 0) {
			f_lift = TRUE;
			prefix_lift = argv[++i];
			cout << "-lift " << prefix_lift << endl;
			}
		else if (strcmp(argv[i], "-lex") == 0) {
			f_lex = TRUE;
			cout << "-lex " << endl;
			}
		else if (strcmp(argv[i], "-split") == 0) {
			f_split = TRUE;
			split_r = atoi(argv[++i]);
			split_m = atoi(argv[++i]);
			cout << "-split " << split_r << " " << split_m << endl;
			}
		else if (strcmp(argv[i], "-solve") == 0) {
			f_solve = TRUE;
			cout << "-solve " << endl;
			}
		else if (strcmp(argv[i], "-save") == 0) {
			f_save = TRUE;
			cout << "-save " << endl;
			}
		else if (strcmp(argv[i], "-read") == 0) {
			f_read_instead = TRUE;
			cout << "-read " << endl;
			}
		else if (strcmp(argv[i], "-draw_system") == 0) {
			f_draw_system = TRUE;
			fname_system = argv[++i];
			cout << "-draw_system " << fname_system << endl;
			}
		else if (strcmp(argv[i], "-write_tree") == 0) {
			f_write_tree = TRUE;
			fname_tree = argv[++i];
			cout << "-write_tree " << fname_tree << endl;
			}
		else if (strcmp(argv[i], "-solution_directory_name") == 0) {
			f_solution_directory_name = TRUE;
			solution_directory_name = argv[++i];
			cout << "-solution_directory_name " << solution_directory_name << endl;
			}
		else if (strcmp(argv[i], "-read_solutions") == 0) {
			f_read_solutions = TRUE;
			cout << "-read_solutions " << endl;
			}
		else if (strcmp(argv[i], "-read_solutions_after_split") == 0) {
			f_read_solutions_after_split = TRUE;
			read_solutions_split_m = atoi(argv[++i]);
			cout << "-read_solutions_after_split " << read_solutions_split_m << endl;
			}
		else if (strcmp(argv[i], "-compute_orbits") == 0) {
			f_compute_orbits = TRUE;
			cout << "-compute_orbits " << endl;
			}
		else if (strcmp(argv[i], "-isomorph_testing") == 0) {
			f_isomorph_testing = TRUE;
			cout << "-isomorph_testing " << endl;
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
		}

	INT f_v = (verbose_level >= 1);

	if (f_memory_debug) {
		start_memory_debug();
		}
	{
	regular_ls_generator Gen;

	Gen.init_basic(argc, argv, verbose_level);
	
	Gen.init_group(verbose_level);


	Gen.init_action_on_k_subsets(Gen.k, 0 /*verbose_level*/);

	Gen.init_generator(Gen.n, FALSE, NULL, 
		Gen.A->Strong_gens, verbose_level);


	if (f_v) {
		cout << "init finished" << endl;
		}


	if (f_starter) {

		INT f_write_candidate_file = TRUE;
		
		if (!starter_size) {
			cout << "please use option -starter_size <starter_size>" << endl;
			exit(1);
			}
		Gen.compute_starter(starter_size, 
			starter_directory_name, 
			f_lex, f_write_candidate_file, 
			f_draw_poset, f_embedded, verbose_level);

		}

	else if (f_build_db) {
		
		if (!starter_size) {
			cout << "please use option -starter_size <starter_size>" << endl;
			exit(1);
			}

		sprintf(Gen.gen->fname_base, "%s%s", starter_directory_name, Gen.base_fname);

		isomorph_build_db(
			Gen.A, Gen.A2, Gen.gen,
			Gen.n /*target_size*/, 
			Gen.gen->fname_base, (BYTE *)"ISO/", 
			starter_size, verbose_level);
		}

	else if (f_lift) {
		
		cout << "Gen.base_fname=" << Gen.base_fname << endl;
		
		if (!starter_size) {
			cout << "please use option -starter_size <starter_size>" << endl;
			exit(1);
			}
		if (f_starter_directory_name == FALSE) {
			cout << "please use option -starter_directory_name <starter_directory_name>" << endl;
			exit(1);
			}
		if (f_solution_directory_name == FALSE) {
			cout << "please use option -solution_directory_name <solution_directory_name>" << endl;
			exit(1);
			}

		BYTE cmd[1000];

		sprintf(cmd, "mkdir %s", solution_directory_name);
		system(cmd);
		sprintf(cmd, "mkdir %s", prefix_lift);
		system(cmd);
		
		
		compute_lifts_new(
			Gen.A, Gen.A2, (void *) &Gen,
			Gen.base_fname, 
			starter_directory_name /* input_prefix */, prefix_lift /* output_prefix */, solution_directory_name /* solution_prefix */, 
			starter_size, Gen.n /* target_size */, 
			f_lex, f_split, split_r, split_m, 
			f_solve, f_save, f_read_instead, 
			f_draw_system, fname_system, 
			f_write_tree, fname_tree,
			rls_generator_lifting_prepare_function_new, 
			rls_generator_early_test_function, 
			(void *) &Gen, 
			FALSE /* f_has_solution_test_function */, 
			NULL,  
			(void *) &Gen,
			FALSE /* f_has_late_cleanup_function */, 
			NULL, 
			verbose_level);


		}

	else if (f_read_solutions) {
		BYTE *fname[1];
		BYTE fname1[1000];
		INT nb_files = 1;
		
		if (!starter_size) {
			cout << "please use option -starter_size <starter_size>" << endl;
			exit(1);
			}
		if (f_starter_directory_name == FALSE) {
			cout << "please use option -starter_directory_name <starter_directory_name>" << endl;
			exit(1);
			}
		if (f_solution_directory_name == FALSE) {
			cout << "please use option -solution_directory_name <solution_directory_name>" << endl;
			exit(1);
			}

		sprintf(Gen.gen->fname_base, "%s%s", starter_directory_name, Gen.base_fname);

		sprintf(fname1, "%s%s_depth_%ld_solutions.txt", 
			solution_directory_name, Gen.base_fname, starter_size);
		fname[0] = fname1;

		isomorph_read_solution_files(Gen.A, Gen.A2, Gen.gen, 
			Gen.n, Gen.gen->fname_base, (BYTE *)"ISO/", starter_size, 
			(const BYTE **) fname, nb_files, verbose_level);

		}
	else if (f_read_solutions_after_split) {
		BYTE **fname;
		BYTE fname1[1000];
		INT nb_files = 0;

		if (!starter_size) {
			cout << "please use option -starter_size <starter_size>" << endl;
			exit(1);
			}
		if (f_starter_directory_name == FALSE) {
			cout << "please use option -starter_directory_name <starter_directory_name>" << endl;
			exit(1);
			}
		if (f_solution_directory_name == FALSE) {
			cout << "please use option -solution_directory_name <solution_directory_name>" << endl;
			exit(1);
			}

		sprintf(Gen.gen->fname_base, "%s%s", starter_directory_name, Gen.base_fname);

		nb_files = read_solutions_split_m;
		fname = NEW_PBYTE(nb_files);
		for (i = 0; i < read_solutions_split_m; i++) {
			sprintf(fname1, "%s%s_depth_%ld_split_%ld_%ld_solutions.txt", 
				solution_directory_name, Gen.base_fname, starter_size, i, read_solutions_split_m);
			fname[i] = NEW_BYTE(strlen(fname1) + 1);
			strcpy(fname[i], fname1);
			}
		cout << "Reading the following " << nb_files << " files:" << endl;
		for (i = 0; i < nb_files; i++) {
			cout << i << " : " << fname[i] << endl;
			}
		isomorph_read_solution_files(Gen.A, Gen.A2, Gen.gen, 
			Gen.n, Gen.gen->fname_base, (BYTE *)"ISO/", starter_size, 
			(const BYTE **) fname, nb_files, verbose_level);
		for (i = 0; i < nb_files; i++) {
			FREE_BYTE(fname[i]);
			}
		FREE_PBYTE(fname);
		}

	else if (f_compute_orbits) {

		if (!starter_size) {
			cout << "please use option -starter_size <starter_size>" << endl;
			exit(1);
			}
		if (f_starter_directory_name == FALSE) {
			cout << "please use option -starter_directory_name <starter_directory_name>" << endl;
			exit(1);
			}

		sprintf(Gen.gen->fname_base, "%s%s", starter_directory_name, Gen.base_fname);


		isomorph_compute_orbits(Gen.A, Gen.A2, Gen.gen, 
			Gen.n, Gen.gen->fname_base, (BYTE *)"ISO/", starter_size, verbose_level);
		}
	else if (f_isomorph_testing) {
			
		if (!starter_size) {
			cout << "please use option -starter_size <starter_size>" << endl;
			exit(1);
			}
		if (f_starter_directory_name == FALSE) {
			cout << "please use option -starter_directory_name <starter_directory_name>" << endl;
			exit(1);
			}

		sprintf(Gen.gen->fname_base, "%s%s", starter_directory_name, Gen.base_fname);

		isomorph_testing(Gen.A, Gen.A2, Gen.gen, 
			Gen.n, Gen.gen->fname_base, (BYTE *)"ISO/", starter_size, 
			f_event_file, event_file_name, print_mod, verbose_level);
		}


	}
	if (f_memory_debug) {
		registry_dump_sorted();
		}

	the_end(t0);
	//the_end_quietly(t0);
}



