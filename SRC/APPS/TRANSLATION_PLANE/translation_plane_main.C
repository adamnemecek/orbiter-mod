// translation_plane_main.C
// 
// Anton Betten
// July 9, 2013
//
//
// 
//
//

#include "orbiter.h"
#include "discreta.h"
#include "translation_plane.h"


// global data:

INT t0; // the system time when the program started

//void do_Fano_subplanes(translation_plane &T, INT order, INT verbose_level);

#define MAX_FILES 1000

int main(int argc, const char **argv)
{
	INT i, j;
	INT verbose_level = 0;
	INT f_poly = FALSE;
	const BYTE *poly = NULL;
	INT f_order = FALSE;
	INT order = 0;
	INT f_dim_over_kernel = FALSE;
	INT dim_over_kernel = 0;
	INT f_clique_level = FALSE;
	INT clique_level = -1;
	INT f_make_spread = FALSE;
	INT type_of_spread = 0;
	INT f_starter = FALSE;
	INT f_depth = FALSE;
	INT depth = 0;
	INT f_identify = FALSE;
	INT identify_data[1000];
	INT identify_data_sz = 0;
	INT f_lift = FALSE;
	INT lift_level = FALSE;
	const BYTE *lift_prefix = "";
	INT f_lex = FALSE;
	INT f_build_db = FALSE;
	INT level = FALSE;
	INT f_read_solution_files = FALSE;
	const BYTE *solution_fname[MAX_FILES];
	INT nb_files = 0;
	INT f_compute_orbits = FALSE;
	INT f_isomorph_testing = FALSE;
	INT f_classification_graph = FALSE;
	INT f_CO = FALSE;
	INT f_report = FALSE;
	INT f_make_quotients = FALSE;
	INT f_extend_simple = FALSE;
	INT extend_starter[1000];
	INT starter_size;
	INT f_plane_type_klein = FALSE;
	const BYTE *fname_plane_type_klein;
	INT f_print_spread = FALSE;
	const BYTE *fname_print_spread;
	INT f_HMO = FALSE;
	const BYTE *fname_HMO;
	INT f_down_orbits = FALSE;
	INT down_orbits_level = 0;
	INT f_split = FALSE;
	INT split_r = 0;
	INT split_m = 1;
	INT f_solve = FALSE;
	INT f_save = FALSE;
	INT f_event_file = FALSE; // -e <event file> option
	const BYTE *event_file_name;
	INT print_mod = 1000;
	INT f_Fano = FALSE;
	INT f_recoordinatize = FALSE;
	INT f_print_representatives = FALSE;
	INT representatives_size = 0;
	const BYTE *representatives_fname = NULL;
	INT f_test_identify = FALSE;
	INT identify_level = 0;
	INT identify_nb_times = 0;
	INT f_draw_poset = FALSE;
	INT f_embedded = FALSE;
	INT f_print_data_structure = FALSE;
	INT f_draw_system = FALSE;
	const BYTE *fname_system = NULL;
	INT f_write_tree = FALSE;
	const BYTE *fname_tree = NULL;

	t0 = os_ticks();
	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-v") == 0) {
			verbose_level = atoi(argv[++i]);
			cout << "-v " << verbose_level << endl;
			}
		else if (strcmp(argv[i], "-poly") == 0) {
			f_poly = TRUE;
			poly = argv[++i];
			cout << "-poly " << poly << endl;
			}
		else if (strcmp(argv[i], "-order") == 0) {
			f_order = TRUE;
			order = atoi(argv[++i]);
			cout << "-order " << order << endl;
			}
		else if (strcmp(argv[i], "-dim_over_kernel") == 0) {
			f_dim_over_kernel = TRUE;
			dim_over_kernel = atoi(argv[++i]);
			cout << "-dim_over_kernel " << dim_over_kernel << endl;
			}
		else if (strcmp(argv[i], "-clique_level") == 0) {
			f_clique_level = TRUE;
			clique_level = atoi(argv[++i]);
			cout << "-clique_level " << clique_level << endl;
			}
		else if (strcmp(argv[i], "-FTWKB") == 0) {
			f_make_spread = TRUE;
			type_of_spread = SPREAD_OF_TYPE_FTWKB;
			cout << "-FTWKB" << endl;
			}
		else if (strcmp(argv[i], "-Kantor") == 0) {
			f_make_spread = TRUE;
			type_of_spread = SPREAD_OF_TYPE_KANTOR;
			cout << "-Kantor" << endl;
			}
		else if (strcmp(argv[i], "-DicksonKantor") == 0) {
			f_make_spread = TRUE;
			type_of_spread = SPREAD_OF_TYPE_DICKSON_KANTOR;
			cout << "-DicksonKantor" << endl;
			}
		else if (strcmp(argv[i], "-Hudson") == 0) {
			f_make_spread = TRUE;
			type_of_spread = SPREAD_OF_TYPE_HUDSON;
			cout << "-Hudson" << endl;
			}
		else if (strcmp(argv[i], "-Kantor2") == 0) {
			f_make_spread = TRUE;
			type_of_spread = SPREAD_OF_TYPE_KANTOR2;
			cout << "-Kantor2" << endl;
			}
		else if (strcmp(argv[i], "-Ganley") == 0) {
			f_make_spread = TRUE;
			type_of_spread = SPREAD_OF_TYPE_GANLEY;
			cout << "-Ganley" << endl;
			}
		else if (strcmp(argv[i], "-Law_Penttila") == 0) {
			f_make_spread = TRUE;
			type_of_spread = SPREAD_OF_TYPE_LAW_PENTTILA;
			cout << "-Law_Penttila" << endl;
			}
		else if (strcmp(argv[i], "-starter") == 0) {
			f_starter = TRUE;
			cout << "-starter " << endl;
			}
		else if (strcmp(argv[i], "-depth") == 0) {
			f_depth = TRUE;
			depth = atoi(argv[++i]);
			cout << "-depth " << depth << endl;
			}
		else if (strcmp(argv[i], "-identify") == 0) {
			INT a;
			
			f_identify = TRUE;
			j = 0;
			while (TRUE) {
				a = atoi(argv[++i]);
				if (a == -1) {
					break;
					}
				identify_data[j++] = a;
				}
			identify_data_sz = j;
			cout << "-identify ";
			INT_vec_print(cout, identify_data, identify_data_sz);
			cout << endl;
			}
		else if (strcmp(argv[i], "-test_identify") == 0) {
			f_test_identify = TRUE;
			identify_level = atoi(argv[++i]);
			identify_nb_times = atoi(argv[++i]);
			cout << "-test_identify " << identify_level << " " << identify_nb_times << endl;
			}
		else if (strcmp(argv[i], "-lift") == 0) {
			f_lift = TRUE;
			lift_level = atoi(argv[++i]);
			lift_prefix = argv[++i]; 
			cout << "-lift " << lift_level << " " << lift_prefix << endl;
			}
		else if (strcmp(argv[i], "-lex") == 0) {
			f_lex = TRUE;
			cout << "-lex" << endl;
			}
		else if (strcmp(argv[i], "-build_db") == 0) {
			f_build_db = TRUE;
			level = atoi(argv[++i]);
			cout << "-build_db " << level << endl;
			}
		else if (strcmp(argv[i], "-read_solution_files") == 0) {
			f_read_solution_files = TRUE;
			level = atoi(argv[++i]);
			i++;
			nb_files = 0;
			while (i < argc) {
				solution_fname[nb_files] = argv[i];
				cout << "solution_fname[nb_files]=" << solution_fname[nb_files] << endl;
				if (strcmp(solution_fname[nb_files], "-1") == 0) {
					break;
					}
				nb_files++;
				i++;
				}
			cout << "-read_solution_files ";
			for (j = 0; j < nb_files; j++) {
				cout << solution_fname[j] << " ";
				}
			cout << endl;
			}
		else if (strcmp(argv[i], "-compute_orbits") == 0) {
			f_compute_orbits = TRUE;
			level = atoi(argv[++i]);
			cout << "-compute_orbits " << level << endl;
			}
		else if (strcmp(argv[i], "-isomorph_testing") == 0) {
			f_isomorph_testing = TRUE;
			level = atoi(argv[++i]);
			cout << "-isomorph_testing " << level << endl;
			}
		else if (strcmp(argv[i], "-classification_graph") == 0) {
			f_classification_graph = TRUE;
			level = atoi(argv[++i]);
			cout << "-classification_graph " << level << endl;
			}
		else if (strcmp(argv[i], "-CO") == 0) {
			f_CO = TRUE;
			level = atoi(argv[++i]);
			cout << "-CO " << level << endl;
			}
		else if (strcmp(argv[i], "-report") == 0) {
			f_report = TRUE;
			cout << "-report " << endl;
			}
		else if (strcmp(argv[i], "-make_quotients") == 0) {
			f_make_quotients = TRUE;
			cout << "-make_quotients " << endl;
			}
		else if (strcmp(argv[i], "-extend_simple") == 0) {
			f_extend_simple = TRUE;
			i++;
			starter_size = 0;
			while (i < argc) {
				extend_starter[starter_size] = atoi(argv[i]);
				if (extend_starter[starter_size] == -1) {
					break;
					}
				starter_size++;
				i++;
				}
			cout << "-extend_simple ";
			for (j = 0; j < starter_size; j++) {
				cout << extend_starter[j] << " ";
				}
			cout << endl;
			}
		else if (strcmp(argv[i], "-plane_type_klein") == 0) {
			f_plane_type_klein = TRUE;
			fname_plane_type_klein = argv[++i];
			cout << "-plane_type_klein " << fname_plane_type_klein << endl;
			}
		else if (strcmp(argv[i], "-print_spread") == 0) {
			f_print_spread = TRUE;
			fname_print_spread = argv[++i];
			cout << "-print_spread " << fname_print_spread << endl;
			}
		else if (strcmp(argv[i], "-HMO") == 0) {
			f_HMO = TRUE;
			fname_HMO = argv[++i];
			cout << "-HMO " << fname_HMO << endl;
			}
		else if (strcmp(argv[i], "-down_orbits") == 0) {
			f_down_orbits = TRUE;
			down_orbits_level = atoi(argv[++i]);
			cout << "-down_orbits " << down_orbits_level << endl;
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
		else if (strcmp(argv[i], "-Fano") == 0) {
			f_Fano = TRUE;
			cout << "-Fano " << endl;
			}
		else if (strcmp(argv[i], "-recoordinatize") == 0) {
			f_recoordinatize = TRUE;
			cout << "-recoordinatize " << endl;
			}
		else if (strcmp(argv[i], "-print_representatives") == 0) {
			f_print_representatives = TRUE;
			representatives_size = atoi(argv[++i]);
			representatives_fname = argv[++i];
			cout << "-print_representatives" << representatives_size << " " << representatives_fname << endl;
			}
		else if (strcmp(argv[i], "-draw_poset") == 0) {
			f_draw_poset = TRUE;
			cout << "-draw_poset " << endl;
			}
		else if (strcmp(argv[i], "-embedded") == 0) {
			f_embedded = TRUE;
			cout << "-embedded " << endl;
			}
		else if (strcmp(argv[i], "-print_data_structure") == 0) {
			f_print_data_structure = TRUE;
			cout << "-print_data_structure " << endl;
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
		}

	if (!f_order) {
		cout << "please use option -order <order>" << endl;
		exit(1);
		}

	INT p, e, e1, n, k, q;
	
	factor_prime_power(order, p, e);
	cout << "order = " << order << " = " << p << "^" << e << endl;

	if (f_dim_over_kernel) {
		if (e % dim_over_kernel) {
			cout << "dim_over_kernel does not divide e" << endl;
			exit(1);
			}
		e1 = e / dim_over_kernel;
		n = 2 * dim_over_kernel;
		k = dim_over_kernel;
		q = i_power_j(p, e1);
		cout << "order=" << order << " n=" << n << " k=" << k << " q=" << q << endl;
		}
	else {
		n = 2 * e;
		k = e;
		q = p;
		cout << "order=" << order << " n=" << n << " k=" << k << " q=" << q << endl;
		}

	INT f_v = (verbose_level >= 1);
	finite_field *F;
	translation_plane T;

	F = new finite_field;

	F->init_override_polynomial(q, poly, 0 /* verbose_level */);

	T.read_arguments(argc, argv);
	

	T.init(order, n, k, F, f_recoordinatize, 0 /*MINIMUM(verbose_level - 1, 2)*/);
	
	T.init2(0 /*verbose_level*/);

	if (clique_level >= 0) {
		translation_plane_init_clique(&T, T.gen, clique_level, verbose_level);
		}
	
	
	if (f_make_spread) {
		T.write_spread_to_file(type_of_spread, verbose_level);
		}
	else if (f_starter) {
		if (!f_depth) {
			cout << "Please use option -depth <depth>" << endl;
			exit(1);
			}
		T.compute(depth, verbose_level);

#if 0
		BYTE fname[1000];

		sprintf(fname, "%s_lvl_%ld", T.gen->fname_base, depth);
		//T.gen->A->read_file_and_print_representatives(fname, FALSE);
#endif
		cout << "depth = " << depth << endl;
		cout << "spread_size = " << T.spread_size << endl;
	

		if (f_draw_poset) {
			if (f_v) {
				cout << "before gen->draw_poset" << endl;
				}
			T.gen->draw_poset(T.gen->fname_base, depth, 0 /* data1 */, f_embedded, verbose_level);
			}


		if (f_print_data_structure) {
			if (f_v) {
				cout << "before gen->print_data_structure_tex" << endl;
				}
			T.gen->print_data_structure_tex(depth, 0 /*gen->verbose_level*/);
			}


#if 0
		if (f_identify) {
			T.identify(identify_data, identify_data_sz, verbose_level);
			}
#endif

		}
	else if (f_identify) {
		if (!f_depth) {
			cout << "Please use option -depth <depth>" << endl;
			exit(1);
			}
		cout << "classifying translation planes" << endl;
		T.compute(order + 1, 0 /* verbose_level */);
		cout << "classifying translation planes done" << endl;

		//T.gen->print_node(5);
		INT *transporter;
		INT orbit_at_level;
		
		transporter = NEW_INT(T.gen->A->elt_size_in_INT);
		
		T.gen->identify(identify_data, identify_data_sz, transporter, orbit_at_level, verbose_level);

		FREE_INT(transporter);
		}
#if 0
	else if (f_Fano) {
		do_Fano_subplanes(T, order, verbose_level);
		}
#endif
	else if (f_test_identify) {
		if (!f_depth) {
			cout << "Please use option -depth <depth>" << endl;
			exit(1);
			}
		cout << "classifying translation planes" << endl;
		T.compute(order + 1, 0 /* verbose_level */);
		cout << "classifying translation planes done" << endl;

		T.gen->test_identify(identify_level, identify_nb_times, verbose_level);
		}
	else if (f_lift) {
		//T.compute(verbose_level);

		compute_lifts(T.A, T.A2, (void *) &T, 
			T.gen->fname_base, 
			lift_prefix, lift_prefix, lift_prefix, 
			lift_level, T.spread_size, 
			f_lex, f_split, split_r, split_m, 
			f_solve, f_save, FALSE /*f_read_instead*/, 
			f_draw_system, fname_system, 
			f_write_tree, fname_tree, 
			translation_plane_lifting_prepare_function,
			translation_plane_lifting_cleanup_function,
			translation_plane_lifting_early_test_function, 
			(void *) &T, 
			FALSE,  NULL, NULL,
			FALSE, NULL, 
			verbose_level);
			// TOP_LEVEL/extra.C

		}
	else if (f_build_db) {
		system("mkdir ISO");
		isomorph_build_db(T.A, T.A2, T.gen, 
			order + 1, T.gen->fname_base, (BYTE *)"ISO/", level, verbose_level);
		}
	else if (f_read_solution_files) {
		isomorph_read_solution_files(T.A, T.A2, T.gen, 
			order + 1 /* target_size */, T.gen->fname_base, (BYTE *)"ISO/", level, 
			solution_fname, nb_files, verbose_level);
		}
	else if (f_compute_orbits) {
		isomorph_compute_orbits(T.A, T.A2, T.gen, 
			order + 1 /* target_size */, T.gen->fname_base, (BYTE *)"ISO/", level, verbose_level);
		}
	else if (f_isomorph_testing) {
		isomorph_testing(T.A, T.A2, T.gen, 
			order + 1 /* target_size */, T.gen->fname_base, (BYTE *)"ISO/", level, 
			f_event_file, event_file_name, print_mod, verbose_level);
		}
	else if (f_classification_graph) {
		isomorph_classification_graph(T.A, T.A2, T.gen, 
			order + 1 /* target_size */, 
			T.gen->fname_base, (BYTE *)"ISO/", 
			level, 
			verbose_level);
		}
	else if (f_CO) {
		T.czerwinski_oakden(level, verbose_level);
		}
	else if (f_report) {
		if (!f_depth) {
			cout << "Please use option -depth <depth>" << endl;
			exit(1);
			}
		isomorph_worker(T.A, T.A2, T.gen, 
			order + 1 /* target_size */, T.gen->fname_base, (BYTE *)"ISO/", 
			translation_plane_callback_report, &T, 
			depth, verbose_level);

		//T.print_classification(order, level, f_select, select_first, select_len, verbose_level);
		}
	else if (f_make_quotients) {
		if (!f_depth) {
			cout << "Please use option -depth <depth>" << endl;
			exit(1);
			}
		isomorph_worker(T.A, T.A2, T.gen, 
			order + 1 /* target_size */, T.gen->fname_base, (BYTE *)"ISO/", 
			translation_plane_callback_make_quotients, &T, 
			depth, verbose_level);

		}
	else if (f_extend_simple) {
		translation_plane_extend_simple(&T, starter_size, extend_starter, FALSE /*f_lex*/, 
			FALSE /* f_write_graph_file */, 
			FALSE /* f_draw_graph */, 
			FALSE /* f_write_tree */, 
			FALSE /* f_decision_nodes_only */, 
			verbose_level);
		}
	else if (f_plane_type_klein) {
		T.test_plane_intersection_type_of_klein_image(
			fname_plane_type_klein, verbose_level);
		}
	else if (f_print_spread) {
		T.read_and_print_spread(fname_print_spread, verbose_level);
		}
	else if (f_HMO) {
		T.HMO(fname_HMO, verbose_level);
		}
	else if (f_down_orbits) {
		isomorph_compute_down_orbits(T.A, T.A2, T.gen, 
			T.spread_size, 
			T.gen->fname_base, (BYTE *)"ISO/", 
			&T, 
			down_orbits_level, verbose_level);
		}
	if (f_print_representatives) {
		orbit_rep *R;
		INT *M;
		INT no, nb;
		BYTE fname[1000];
		
		R = new orbit_rep;
		M = NEW_INT(T.k * T.n);

		sprintf(fname, "%s_lvl_%ld", representatives_fname, representatives_size);

		nb = count_number_of_orbits_in_file(fname, verbose_level);

		cout << "there are " << nb << " orbit representatives in the file " << fname << endl;
		for (no = 0; no < nb; no++) {
			R->init_from_file(T.A /*A_base*/, (BYTE *) representatives_fname, 
				representatives_size, no, representatives_size - 1/*level_of_candidates_file*/, 
				translation_plane_lifting_early_test_function, 
				&T, 
				verbose_level - 1
				);
			// R has: INT *candidates; INT nb_candidates;
	
			for (i = 0; i < representatives_size; i++) {
				cout << R->rep[i] << " ";
				}
			cout << endl;
			for (i = 0; i < representatives_size; i++) {
				cout << R->rep[i] << " = " << endl;
				T.Grass->unrank_INT_here(M, R->rep[i], 0/*verbose_level - 4*/);
				INT_matrix_print(M, T.k, T.n);
				}
			}
		}


//end:
	//the_end(t0);
	the_end_quietly(t0);
}

#if 0
void do_Fano_subplanes(translation_plane &T, INT order, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, j;
	BYTE fname[1000];

	if (f_v) {
		cout << "do_Fano_subplanes" << endl;
		}
	cout << "T.target_size = " << T.target_size << endl;
	
	T.compute(T.target_size, verbose_level);

	sprintf(fname, "%s_lvl_%ld", T.gen->fname_base, T.target_size);
	//T.gen->A->read_file_and_print_representatives(fname, FALSE);

	


	INT *Reps;
	BYTE **Aut_ascii;
	INT nb_reps, size;
			
	cout << "Reading spreads from file " << fname << endl;
	T.gen->A->read_representatives_and_strong_generators(fname, Reps, Aut_ascii, nb_reps, size, 0 /* verbose_level */);
	cout << "Read spreads from file " << fname << endl;


	translation_plane_via_andre_model *TP;
	INT *Nb_subplanes;
	//INT **Subplanes;

	TP = new translation_plane_via_andre_model[nb_reps];


	Nb_subplanes = NEW_INT(nb_reps);
	//Subplanes = NEW_PINT(nb_reps);
	//INT depth = 4;
	INT depth = 7;
			
	for (i = 0; i < nb_reps; i++) {
		cout << endl;
		cout << endl;
		cout << endl;
		cout << endl;
		cout << endl;
		cout << endl;
		cout << "Making translation plane from representative " << i << ":" << endl;

		vector_ge *gens;
		INT *tl;
				
		T.gen->A->get_generators_from_ascii_coding(Aut_ascii[i], gens, tl, verbose_level - 2);



		TP[i].init(Reps + i * T.target_size, T.q, T.k, T.F, gens, tl, verbose_level);

		BYTE prefix[1000];

		sprintf(prefix, "TP_%ld_", i);

		//TP[i].classify_arcs(prefix, depth, verbose_level);
		TP[i].classify_subplanes(prefix, verbose_level);

		Nb_subplanes[i] = TP[i].arcs->number_of_orbits_at_depth(depth);


		FREE_OBJECT(gens);
		FREE_INT(tl);
		}

	for (i = 0; i < nb_reps; i++) {
		cout << "Translation plane " << i << ":" << endl;
		TP[i].arcs->print_orbit_numbers(depth);
		}


#if 0
	for (i = 0; i < nb_reps; i++) {
		cout << "Translation plane " << i << ":" << endl;
		BYTE fname[1000];
		INT *Quadrangles;
		INT nb_quadrangles;
		INT size;

		sprintf(fname, "%s_lvl_%ld", TP[i].arcs->fname_base, depth);
		TP[i].An1->read_representatives(fname, Quadrangles, nb_quadrangles, size, verbose_level);
		Nb_subplanes[i] = 0;
		Subplanes[i] = NEW_INT(nb_quadrangles * 7);
		cout << "Found " << nb_quadrangles << " quadrangles, testing for subplanes" << endl;
		for (j = 0; j < nb_quadrangles; j++) {
			if (TP[i].check_if_quadrangle_defines_a_subplane(Quadrangles + j * 4, Subplanes[i] + 7 * Nb_subplanes[i], verbose_level)) {
				Nb_subplanes[i]++;
				}
			}
		}
	for (i = 0; i < nb_reps; i++) {
		cout << "Translation plane " << i << " has " << Nb_subplanes[i] << " quadrangles that span subplanes" << endl;
		cout << "The subplanes are:" << endl;
		INT_matrix_print(Subplanes[i], Nb_subplanes[i], 7);
		}


	for (i = 0; i < nb_reps; i++) {
		cout << "Translation plane " << i << " has " << Nb_subplanes[i] << " quadrangles that span subplanes" << endl;
		}
#endif
	for (i = 0; i < nb_reps; i++) {
		cout << "Translation plane " << i << " has " << Nb_subplanes[i] << " orbits of Fano subplanes" << endl;
		TP[i].arcs->compute_and_print_automorphism_group_orders(depth, cout);


		if (i == 4 || i == 5) {
			INT fst, node, size;
			INT *set;

			set = NEW_INT(depth);
			fst = TP[i].arcs->first_oracle_node_at_level[depth];
			for (j = 0; j < Nb_subplanes[i]; j++) {
				longinteger_object stab_order;
				node = fst + j;
				TP[i].arcs->stabilizer_order(node, stab_order);
				if ((stab_order.as_INT() % 14) == 0) {
					cout << "orbit " << j << " has stabilizer of order " << stab_order << " which is divisible by 14" << endl;
					TP[i].arcs->get_set(node, set, size);
					if (size != depth) {
						cout << "size != depth" << endl;
						exit(1);
						}
					cout << "The subplane is : ";
					INT_vec_print(cout, set, depth);
					cout << endl;


					vector_ge *stab_gens;
					INT *stab_tl;
							
					stab_gens = new vector_ge;
					stab_tl = NEW_INT(TP[i].An1->base_len);

					TP[i].arcs->get_stabilizer(stab_gens, stab_tl, depth, j, 0 /*verbose_level */);
					sims *S;
					longinteger_object Sgo;

					S = create_sims_from_generators_with_target_group_order_factorized(TP[i].An1, 
						stab_gens, stab_tl, TP[i].An1->base_len, 0 /* verbose_level */);
					S->group_order(Sgo);
					cout << "created group of order " << Sgo << endl;

					delete S;
					delete stab_gens;
					FREE_INT(stab_tl);
					}
				}
			FREE_INT(set);
			}

		cout << endl;
		}

	delete [] TP;
	FREE_INT(Nb_subplanes);

#if 0
	for (i = 0; i < nb_reps; i++) {
		FREE_INT(Subplanes[i]);
		}
	FREE_PINT(Subplanes);
#endif

	FREE_INT(Reps);
	for (i = 0; i < nb_reps; i++) {
		FREE_BYTE(Aut_ascii[i]);
		}
	FREE_PBYTE(Aut_ascii);
		
	if (f_v) {
		cout << "do_Fano_subplanes done" << endl;
		}

}
#endif



