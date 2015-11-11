// arc_generator.C
// 
// Anton Betten
//
// previous version Dec 6, 2004
// revised June 19, 2006
// revised Aug 17, 2008
// moved here from hyperoval.C May 10, 2013
//
// Searches for arcs and hyperovals in Desarguesian projective planes
//
//

#include "orbiter.h"
#include "discreta.h"
#include "hyperoval.h"

void arc_generator::read_arguments(int argc, const char **argv)
{
	verbose_level = 0;
	INT i, j;
	INT f_q = FALSE;
	//INT f_depth = FALSE;
	INT f_target_size = FALSE;

	for (i = 1; i < argc; i++) {
		
		if (strcmp(argv[i], "-v") == 0) {
			verbose_level = atoi(argv[++i]);
			cout << "-v " << verbose_level << endl;
			}
		else if (strcmp(argv[i], "-graph") == 0) {
			f_write_graph = TRUE;
			cout << "-graph " << endl;
			}
		else if (strcmp(argv[i], "-maxdepth") == 0) {
			f_maxdepth = TRUE;
			maxdepth = atoi(argv[++i]);
			cout << "-maxdepth " << maxdepth << endl;
			}
		else if (strcmp(argv[i], "-target_size") == 0) {
			f_target_size = TRUE;
			target_size = atoi(argv[++i]);
			cout << "-target_size " << target_size << endl;
			}
		else if (strcmp(argv[i], "-q") == 0) {
			f_q = TRUE;
			q = atoi(argv[++i]);
			cout << "-q " << q << endl;
			}
		else if (strcmp(argv[i], "-compute_starter") == 0) {
			f_compute_starter = TRUE;
			cout << "-compute_starter " << endl;
			}
		else if (strcmp(argv[i], "-starter_size") == 0) {
			f_starter_size = TRUE;
			starter_size = atoi(argv[++i]);
			cout << "-starter_size " << starter_size << endl;
			}
		else if (strcmp(argv[i], "-build_db") == 0) {
			f_build_db = TRUE;
			cout << "-build_db " << endl;
			}
		else if (strcmp(argv[i], "-draw_poset") == 0) {
			f_draw_poset = TRUE;
			cout << "-draw_poset " << endl;
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
		else if (strcmp(argv[i], "-read_solutions") == 0) {
			f_read_solutions = TRUE;
			cout << "-read_solutions " << endl;
			}
		else if (strcmp(argv[i], "-compute_orbits") == 0) {
			f_compute_orbits = TRUE;
			cout << "-compute_orbits " << endl;
			}
		else if (strcmp(argv[i], "-classify") == 0) {
			f_classify = TRUE;
			cout << "-classify " << endl;
			}
		else if (strcmp(argv[i], "-decomposition_matrix") == 0) {
			f_decomposition_matrix = TRUE;
			cout << "-decomposition_matrix " << endl;
			}
		else if (strcmp(argv[i], "-report") == 0) {
			f_report = TRUE;
			cout << "-report " << endl;
			}
		else if (strcmp(argv[i], "-solution_file") == 0) {
			solution_file[nb_solution_files] = argv[++i];
			cout << "-solution_file " << solution_file[nb_solution_files] << endl;
			nb_solution_files++;
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
		else if (strcmp(argv[i], "-no_arc_testing") == 0) {
			f_no_arc_testing = TRUE;
			cout << "-no_arc_testing " << endl;
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
		else if (strcmp(argv[i], "-recognize") == 0) {
			f_recognize = TRUE;
			recognize_set_sz = atoi(argv[++i]);
			recognize_set = NEW_INT(recognize_set_sz);
			for (j = 0; j < recognize_set_sz; j++) {
				recognize_set[j] = atoi(argv[++i]);
				}
			cout << "-recognize ";
			INT_vec_print(cout, recognize_set, recognize_set_sz);
			cout << endl;
			}
		}
	if (!f_q) {
		cout << "Please specify the field size using the option -q <q>" << endl;
		exit(1);
		}
	if (!f_target_size) {
		cout << "Please specify the target size using the option -target_size <target_size>" << endl;
		exit(1);
		}
}


void arc_generator::main(int argc, const char **argv)
{
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "arc_generator::main" << endl;
		}
	if (f_compute_starter) {

		prepare_prefix_starter();

		compute_starter(argc, argv, 
			f_recognize, recognize_set, recognize_set_sz, 
			verbose_level);
		
		}
	else if (f_build_db) {

		build_db(argc, argv, verbose_level);
		
		}
	else if (f_lift) {

		if (!f_starter_directory_name) {
			cout << "use option -starter_directory_name <starter_directory_name>" << endl;
			exit(1);
			}
		if (!f_solution_directory_name) {
			cout << "use option -solution_directory_name <solution_directory_name>" << endl;
			exit(1);
			}

		cout << "base_fname=" << base_fname << endl;
		
		compute_lifts_new(
			A, A, (void *) this,
			base_fname, 
			prefix_starter /* input_prefix */, prefix_lift /* output_prefix */, solution_directory_name /* solution_prefix */, 
			starter_size, target_size, 
			f_lex, f_split, split_r, split_m, 
			f_solve, f_save, f_read_instead, 
			f_draw_system, fname_system, 
			f_write_tree, fname_tree,
			arc_generator_lifting_prepare_function_new, 
			arc_generator_early_test_function, 
			(void *) this, 
			TRUE /* f_has_solution_test_function */, 
			callback_arc_test,  
			(void *) this,
			FALSE /* f_has_late_cleanup_function */, 
			NULL, 
			verbose_level);


		}
	else if (f_read_solutions) {
		BYTE prefix_iso[1000];
		BYTE fname_base[1000];
		generator *gen;


		prepare_prefix_iso(prefix_iso);
		prepare_fname_base(fname_base);

		gen = prepare_generator(argc, argv, verbose_level);
	
	
		
		isomorph_read_solution_files(A, A, gen, 
			target_size, 
			fname_base, prefix_iso, 
			starter_size, 
			solution_file, nb_solution_files, 
			verbose_level);

		delete gen;
		}
	else if (f_compute_orbits) {
		BYTE prefix_iso[1000];
		BYTE fname_base[1000];
		generator *gen;

		prepare_prefix_iso(prefix_iso);
		prepare_fname_base(fname_base);

		gen = prepare_generator(argc, argv, verbose_level);
		
		isomorph_compute_orbits(A, A, gen, 
			target_size, 
			fname_base, prefix_iso, 
			starter_size, 
			verbose_level);

		delete gen;
		}
	else if (f_classify) {
		BYTE prefix_iso[1000];
		BYTE fname_base[1000];
		generator *gen;

		prepare_prefix_iso(prefix_iso);
		prepare_fname_base(fname_base);

		gen = prepare_generator(argc, argv, verbose_level);
		
		isomorph_testing(A, A, gen, 
			target_size, 
			fname_base, prefix_iso, 
			starter_size, 
			FALSE, NULL, 1000, 
			verbose_level);

		delete gen;
		}
	else if (f_decomposition_matrix) {
		BYTE prefix_iso[1000];
		BYTE fname_base[1000];
		generator *gen;

		prepare_prefix_iso(prefix_iso);
		prepare_fname_base(fname_base);

		gen = prepare_generator(argc, argv, verbose_level);
		
		isomorph_classification_graph(A, A, gen, 
			target_size, 
			fname_base, prefix_iso, 
			starter_size, 
			verbose_level);

		delete gen;
		}
	else if (f_report) {
		BYTE prefix_iso[1000];
		BYTE fname_base[1000];
		generator *gen;

		prepare_prefix_iso(prefix_iso);
		prepare_fname_base(fname_base);

		gen = prepare_generator(argc, argv, verbose_level);
		
		isomorph_worker(A, A, gen, 
			target_size, 
			fname_base, prefix_iso, 
			callback_arc_report, this, 
			starter_size, verbose_level);

		delete gen;
		}
	if (f_v) {
		cout << "arc_generator::main done" << endl;
		}
}


arc_generator::arc_generator()
{
	null();
}

arc_generator::~arc_generator()
{
	if (Grass) {
		delete Grass;
		}
	if (A) {
		delete A;
		}
	if (A_on_lines) {
		delete A_on_lines;
		}
	null();
}

void arc_generator::null()
{

	f_starter_directory_name = FALSE;
	f_solution_directory_name = FALSE;
	Grass = NULL;
	A = NULL;
	A_on_lines = NULL;
	P2 = NULL;
	verbose_level = 0;
	f_maxdepth = FALSE;
	f_write_graph = FALSE;

	f_compute_starter = FALSE;
	f_build_db = FALSE;
	f_draw_poset = FALSE;
	f_lift = FALSE;
	f_read_solutions = FALSE;
	f_compute_orbits = FALSE;
	f_classify = FALSE;
	f_decomposition_matrix = FALSE;
	f_report = FALSE;
	f_recognize = FALSE;
	recognize_set = NULL;
	recognize_set_sz = NULL;

	nb_solution_files = 0;

	f_draw_system = FALSE;
	f_write_tree = FALSE;

	f_lex = FALSE;
	f_split = FALSE;
	split_r = 0;
	split_m = 1;
	f_solve = FALSE;
	f_save = FALSE;
	f_read_instead = FALSE;
	
	f_no_arc_testing = FALSE;
}

void arc_generator::freeself()
{
	if (A) {
		delete A;
		}
	if (A_on_lines) {
		delete A_on_lines;
		}
	if (recognize_set) {
		FREE_INT(recognize_set);
		}
	if (P2) {
		delete P2;
		}
	null();
}

void arc_generator::prepare_prefix_iso(BYTE *prefix_iso)
{
	sprintf(prefix_iso, "ISO/");
}

void arc_generator::prepare_prefix_starter()
{
	if (f_starter_directory_name) {
		sprintf(prefix_starter, "%s", starter_directory_name);
		}
	else {
		sprintf(prefix_starter, "./");
		}
}

void arc_generator::prepare_fname_base(BYTE *fname_base)
{
	sprintf(fname_base, "%s%s", prefix_starter, base_fname);
}

generator *arc_generator::prepare_generator(int argc, const char **argv, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	BYTE fname_base[1000];
	generator *gen;	

	if (f_v) {
		cout << "arc_generator::prepare_generator" << endl;
		cout << "arc_generator::prepare_generator starter_size = " << starter_size << endl;
		}
	prepare_fname_base(fname_base);

	gen = new generator;


	gen->read_arguments(argc, argv, 0);


	gen->depth = starter_size;
	gen->initialize(A, A,  
		A->Strong_gens, 
		starter_size, 
		fname_base, verbose_level - 1);


	if (f_no_arc_testing) {
		gen->init_check_func(placebo_test_function, 
			(void *)this /* candidate_check_data */);
		}
	else {
		gen->init_check_func(::check_arc, 
			(void *)this /* candidate_check_data */);
		}
	if (f_v) {
		cout << "arc_generator::prepare_generator done" << endl;
		}
	return gen;
}

void arc_generator::compute_starter(int argc, const char **argv, 
	INT f_recognize, INT *set, INT sz, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT t0 = os_ticks();
	generator *gen;

	
	BYTE cmd[1000];



	cout << "arc_generator::compute_starter prefix_starter=" << prefix_starter << endl;

	if (!f_starter_size) {
		cout << "please use option -starter_size <starter_size>" << endl;
		exit(1);
		}
	
	sprintf(cmd, "mkdir %s", prefix_starter);
	system(cmd);

	gen = prepare_generator(argc, argv, verbose_level);


	cout << "arc_generator::compute_starter base_fname=" << base_fname << endl;

	gen->f_print_function = TRUE;
	gen->print_function = print_arc;
	gen->print_function_data = this;
	
#if 0
	if (gen->f_extend) {
		do_extend(verbose_level, gen->verbose_level_group_theory);
		time_check(cout, t0);
		cout << endl;
		exit(0);
		}
#endif

	INT schreier_depth = 1000;
	INT f_use_invariant_subset_if_available = TRUE;
	INT f_debug = FALSE;
	INT f_implicit_fusion = FALSE;
	INT depth;
	INT f_embedded = TRUE;

	if (f_v) {
		cout << "arc_generator::compute_starter before generator_main" << endl;
		}

	depth = gen->main(t0, 
		schreier_depth, 
		f_use_invariant_subset_if_available, 
		f_implicit_fusion, 
		f_debug, 
		verbose_level);
	if (f_v) {
		cout << "arc_generator::compute_starter gen->main returns depth=" << depth << endl;
		}

	if (f_v) {
		cout << "arc_generator::compute_starter after gen->main" << endl;
		}


#if 0
	if (f_v) {
		cout << "arc_generator::compute_starter before gen->print_data_structure_tex" << endl;
		}

	//gen->print_data_structure_tex(depth, 0 /*gen->verbose_level */);
#endif

	if (f_draw_poset) {
		if (f_v) {
			cout << "arc_generator::compute_starter before gen->draw_poset" << endl;
			}

		gen->draw_poset(gen->fname_base, depth, 0 /* data1 */, f_embedded, 0 /* gen->verbose_level */);
		}


	if (f_recognize) {
		INT *canonical_set;
		INT *Elt_transporter;
		INT f_implicit_fusion = TRUE;

		cout << "recognize" << endl;
		cout << "input set = ";
		INT_vec_print(cout, set, sz);
		cout << endl;
		
		canonical_set = NEW_INT(sz);
		Elt_transporter = NEW_INT(gen->A->elt_size_in_INT);
		
		gen->trace_set(set, sz, sz /* level */, 
			canonical_set, Elt_transporter, 
			f_implicit_fusion, verbose_level);

		cout << "output set = ";
		INT_vec_print(cout, canonical_set, sz);
		cout << endl;

		
		FREE_INT(canonical_set);
		FREE_INT(Elt_transporter);
		}



	BYTE prefix_iso[1000];
	//BYTE cmd[1000];


	prepare_prefix_iso(prefix_iso);

	sprintf(cmd, "mkdir %s", prefix_iso);
	system(cmd);

	if (f_v) {
		cout << "arc_generator::compute_starter done" << endl;
		}

}

void arc_generator::build_db(int argc, const char **argv, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	BYTE fname_base[1000];

	if (f_v) {
		cout << "arc_generator::build_db" << endl;
		}

	generator *gen;


	gen = prepare_generator(argc, argv, verbose_level);



	gen->f_print_function = TRUE;
	gen->print_function = print_arc;
	gen->print_function_data = this;

	prepare_fname_base(fname_base);

	isomorph_build_db(A, A, gen, 
		target_size, 
		fname_base, (BYTE *)"ISO/", 
		starter_size, verbose_level);

	delete gen;

	if (f_v) {
		cout << "arc_generator::build_db done" << endl;
		}
}

void arc_generator::init(int argc, const char **argv)
{
	F = new finite_field;
	A = new action;
	A_on_lines = new action;
	AG = new action_on_grassmannian;

	read_arguments(argc, argv);
	prepare_prefix_starter();
	
	INT f_semilinear = TRUE;
	
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "arc_generator::init" << endl;
		}
	F->init(q, 0);

	nb_points_total = q * q + q + 1;

	block_size = q + 1;

	sprintf(base_fname, "arc_%ld", q);


	if (is_prime(q)) {
		f_semilinear = FALSE;
		}


	INT f_basis = TRUE;
	if (f_v) {
		cout << "arc_generator::init calling init_matrix_group" << endl;
		}
	A->init_projective_group(3, F, f_semilinear, f_basis, 0 /*verbose_level*/);

	if (f_v) {
		cout << "arc_generator::init calling compute_strong_generators" << endl;
		}


	
	if (f_v) {
		cout << "arc_generator::init group set up" << endl;
		}
	
	if (f_v) {
		cout << "arc_generator::init creating action on lines" << endl;
		}
	Grass = new grassmann;

	Grass->init(3 /*n*/, 2 /*k*/, F, verbose_level - 2);
	AG->init(*A, Grass, verbose_level - 2);
	
	A_on_lines->induced_action_on_grassmannian(A, AG, 
		FALSE /*f_induce_action*/, NULL /*sims *old_G */, 
		MINIMUM(verbose_level - 2, 2));
	
	if (f_v) {
		cout << "action A_on_lines created: ";
		A_on_lines->print_info();
		}

	
	rc.init(F, 3, target_size, 4 /* d */);
	
	if (f_v) {
		cout << "arc_generator::init after rc.init" << endl;
		}

	if (f_v) {
		cout << "arc_generator::init creating projective plane" << endl;
		}


	P2 = new projective_space;

	P2->init(2, F, 
		TRUE /* f_init_incidence_structure */, 
		0 /*verbose_level - 2*/);

	if (f_v) {
		cout << "arc_generator::init after P2->init" << endl;
		}



	if (f_v) {
		cout << "arc_generator::init done" << endl;
		}
}


void arc_generator::early_test_func(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT i, j, rk;
	INT *Pts;
	INT *Candidates;
	INT *Mtx;
		
	if (f_v) {
		cout << "arc_generator::early_test_func checking set ";
		print_set(cout, len, S);
		cout << endl;
		cout << "candidate set of size " << nb_candidates << ":" << endl;
		INT_vec_print(cout, candidates, nb_candidates);
		cout << endl;
		}
	Pts = NEW_INT(len * 3);
	Candidates = NEW_INT(nb_candidates * 3);
	Mtx = NEW_INT(3 * 3);
	for (i = 0; i < len; i++) {
		point_unrank(Pts + i * 3, S[i]);
		}
	for (i = 0; i < nb_candidates; i++) {
		point_unrank(Candidates + i * 3, candidates[i]);
		}

	nb_good_candidates = 0;
	for (j = 0; j < nb_candidates; j++) {
		for (i = 0; i < len - 1; i++) {

			INT_vec_copy(Pts + i * 3, Mtx + 0 * 3, 3);
			INT_vec_copy(Pts + (len - 1) * 3, Mtx + 1 * 3, 3);
			INT_vec_copy(Candidates + j * 3, Mtx + 2 * 3, 3);

			rk = F->rank_of_matrix(Mtx, 3, 0);
			if (rk < 3) {
				if (f_vv) {
					cout << "rank is " << rk << " which is bad" << endl;
					}
				break;
				}
			else {
				if (f_vv) {
					cout << "rank is " << rk << " which is OK" << endl;
					}
				}
			} // next i
		if (i == len - 1) {
			good_candidates[nb_good_candidates++] = candidates[j];
			}
		} // next j
	
	FREE_INT(Pts);
	FREE_INT(Candidates);
	FREE_INT(Mtx);
}

INT arc_generator::check_arc(INT *S, INT len, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	//INT f_vvv = (verbose_level >= 3);
	INT f_OK = TRUE;


	if (f_v) {
		cout << "checking set ";
		print_set(cout, len, S);
		}
	if (!rc.check_rank(len, S, verbose_level - 1)) {
		return FALSE;
		}
	
	if (f_v) {
		cout << "checking set ";
		print_set(cout, len, S);
		}
	if (f_v) {
		cout << endl;
		//print_integer_matrix(cout, S, 1, len);
		print_integer_matrix(cout, rc.M1, rc.m, len);
		if (len > 2) {
			print_set_in_affine_plane(cout, len - 2, S + 2);
			}
		}



	if (f_OK) {
		if (f_v) {
			cout << "accepted" << endl;
			}
		return TRUE;
		}
	else {
		if (f_v) {
			cout << "rejected" << endl;
			}
		return FALSE;
		}
}

void arc_generator::print_set_in_affine_plane(ostream &ost, INT len, INT *S)
{
	::print_set_in_affine_plane(ost, *F, len, S);
}




void arc_generator::point_unrank(INT *v, INT rk)
{
	PG_element_unrank_modified(*F, v, 1 /* stride */, 3 /* len */, rk);
}

INT arc_generator::point_rank(INT *v)
{
	INT rk;
	
	PG_element_rank_modified(*F, v, 1 /* stride */, 3, rk);
	return rk;
}


void arc_generator::lifting_prepare_function_new(exact_cover *E, INT starter_case, 
	INT *candidates, INT nb_candidates, strong_generators *Strong_gens, 
	diophant *&Dio, INT *&col_labels, 
	INT &f_ruled_out, 
	INT verbose_level)
// compute the incidence matrix of tangent lines versus candidate points
// extended by external lines versus candidate points
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT i, j, a, b;
	INT nb_needed;

	if (f_v) {
		cout << "arc_generator::lifting_prepare_function_new nb_candidates=" << nb_candidates << endl;
		}

	nb_needed = target_size - starter_size;
	f_ruled_out = FALSE;


	// compute the line type:

	INT *line_type;

	line_type = NEW_INT(P2->N_lines);
	INT_vec_zero(line_type, P2->N_lines);
	
	for (i = 0; i < starter_size; i++) {
		a = E->starter[i];
		for (j = 0; j < P2->r; j++) {
			b = P2->Lines_on_point[a * P2->r + j];
			line_type[b]++;
			}
		}



	classify C;

	C.init(line_type, P2->N_lines, FALSE, 0);
	if (f_v) {
		cout << "arc_generator::lifting_prepare_function_new line_type:" << endl;
		C.print_naked(TRUE);
		cout << endl;
		}


	// extract the tangent lines:

	INT tangent_lines_fst, nb_tangent_lines;
	INT *tangent_lines;
	INT *tangent_line_idx;
	INT external_lines_fst, nb_external_lines;
	INT *external_lines;
	INT *external_line_idx;
	INT fst, len, idx;


	// find all tangent lines:

	for (i = 0; i < C.nb_types; i++) {
		fst = C.type_first[i];
		len = C.type_len[i];
		idx = C.sorting_perm_inv[fst];
		if (line_type[idx] == 1) {
			break;
			}
		}
	if (i == C.nb_types) {
		cout << "arc_generator::lifting_prepare_function_new there are no tangent lines" << endl;
		exit(1);
		}
	tangent_lines_fst = fst;
	nb_tangent_lines = len;
	tangent_lines = NEW_INT(nb_tangent_lines);
	tangent_line_idx = NEW_INT(P2->N_lines);
	for (i = 0; i < P2->N_lines; i++) {
		tangent_line_idx[i] = -1;
		}
	for (i = 0; i < len; i++) {
		j = C.sorting_perm_inv[tangent_lines_fst + i];
		tangent_lines[i] = j;
		tangent_line_idx[j] = i;
		}


	// find all external lines:
	for (i = 0; i < C.nb_types; i++) {
		fst = C.type_first[i];
		len = C.type_len[i];
		idx = C.sorting_perm_inv[fst];
		if (line_type[idx] == 0) {
			break;
			}
		}
	if (i == C.nb_types) {
		cout << "arc_generator::lifting_prepare_function_new there are no external lines" << endl;
		exit(1);
		}
	external_lines_fst = fst;
	nb_external_lines = len;
	external_lines = NEW_INT(nb_external_lines);
	external_line_idx = NEW_INT(P2->N_lines);
	for (i = 0; i < P2->N_lines; i++) {
		external_line_idx[i] = -1;
		}
	for (i = 0; i < len; i++) {
		j = C.sorting_perm_inv[external_lines_fst + i];
		external_lines[i] = j;
		external_line_idx[j] = i;
		}


	
	col_labels = NEW_INT(nb_candidates);


	INT_vec_copy(candidates, col_labels, nb_candidates);

	if (E->f_lex) {
		E->lexorder_test(col_labels, nb_candidates, Strong_gens->gens, 
			verbose_level - 2);
		}

	if (f_vv) {
		cout << "arc_generator::lifting_prepare_function_new after lexorder test" << endl;
		cout << "arc_generator::lifting_prepare_function_new nb_candidates=" << nb_candidates << endl;
		}

	// compute the incidence matrix between
	// tangent lines and candidate points as well as external lines and candidate points:


	INT nb_rows;
	INT nb_cols;

	nb_rows = nb_tangent_lines + nb_external_lines;
	nb_cols = nb_candidates;

	Dio = new diophant;
	Dio->open(nb_rows, nb_cols);
	Dio->sum = nb_needed;

	for (i = 0; i < nb_tangent_lines; i++) {
		Dio->type[i] = t_EQ;
		Dio->RHS[i] = 1;
		}

	for (i = 0; i < nb_external_lines; i++) {
		Dio->type[nb_tangent_lines + i] = t_ZOR;
		Dio->RHS[nb_tangent_lines + i] = 2;
		}

	Dio->fill_coefficient_matrix_with(0);


	for (i = 0; i < nb_candidates; i++) {
		a = col_labels[i];
		for (j = 0; j < P2->r; j++) {
			b = P2->Lines_on_point[a * P2->r + j];
			if (line_type[b] == 2) {
				cout << "arc_generator::lifting_prepare_function candidate lies on a secant" << endl;
				exit(1);
				}
			idx = tangent_line_idx[b];
			if (idx >= 0) {
				Dio->Aij(idx, i) = 1;
				}
			idx = external_line_idx[b];
			if (idx >= 0) {
				Dio->Aij(nb_tangent_lines + idx, i) = 1;
				}
			}
		}


	FREE_INT(line_type);
	FREE_INT(tangent_lines);
	FREE_INT(tangent_line_idx);
	FREE_INT(external_lines);
	FREE_INT(external_line_idx);
	
	if (f_v) {
		cout << "arc_generator::lifting_prepare_function_new done" << endl;
		}
}

INT arc_generator::arc_test(INT *S, INT len, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT ret = TRUE;
	INT i, a, j, b;

	if (f_v) {
		cout << "arc_generator::arc_test" << endl;
		}

	INT *line_type;

	line_type = NEW_INT(P2->N_lines);
	INT_vec_zero(line_type, P2->N_lines);
	for (i = 0; i < len; i++) {
		a = S[i];
		for (j = 0; j < P2->r; j++) {
			b = P2->Lines_on_point[a * P2->r + j];
			line_type[b]++;
			}
		}
	for (i = 0; i < P2->N_lines; i++) {
		if (line_type[i] > 2) {
			ret = FALSE;
			break;
			}
		}
	FREE_INT(line_type);
	return ret;
}

void arc_generator::report(isomorph &Iso, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	BYTE fname[1000];

	if (f_v) {
		cout << "arc_generator::report" << endl;
		}
	if (target_size == q + 2) {
		sprintf(fname, "hyperovals_%ld.tex", q);
		}
	else {
		sprintf(fname, "arcs_%ld_%ld.tex", q, target_size);
		}

	{
	ofstream f(fname);
	INT f_book = TRUE;
	INT f_title = TRUE;
	BYTE title[1000];
	const BYTE *author = "Orbiter";
	INT f_toc = TRUE;
	INT f_landscape = FALSE;
	INT f_12pt = FALSE;
	INT f_enlarged_page = TRUE;
	INT f_pagenumbers = TRUE;

	if (target_size == q + 2) {
		sprintf(title, "Hyperovals over ${\\mathbb F}_{%ld}$", q);
		}
	else {
		sprintf(title, "Arcs over  ${\\mathbb F}_{%ld}$ of size $%ld$", q, target_size);
		}
	cout << "Writing file " << fname << " with " << Iso.Reps->count << " arcs:" << endl;
	latex_head(f, f_book, f_title, 
		title, author, 
		f_toc, f_landscape, f_12pt, f_enlarged_page, f_pagenumbers);

	f << "\\chapter{Summary}" << endl << endl;
	f << "There are " << Iso.Reps->count << " isomorphism types." << endl << endl;


	Iso.setup_and_open_solution_database(verbose_level - 1);

	INT i, first, c, id;
	INT u, v, h, rep, tt;
	longinteger_object go;
	INT data[1000];



	longinteger_object *Ago, *Ago_induced;
	INT *Ago_INT;

	Ago = new longinteger_object[Iso.Reps->count];
	Ago_induced = new longinteger_object[Iso.Reps->count];
	Ago_INT = NEW_INT(Iso.Reps->count);


	for (h = 0; h < Iso.Reps->count; h++) {
		rep = Iso.Reps->rep[h];
		first = Iso.orbit_fst[rep];
		c = Iso.starter_number[first];
		id = Iso.orbit_perm[first];		
		Iso.load_solution(id, data);

		sims *Stab;
		
		Stab = Iso.Reps->stab[h];

		Iso.Reps->stab[h]->group_order(Ago[h]);
		Ago_INT[h] = Ago[h].as_INT();
		if (f_v) {
			cout << "arc_generator::report computing induced action on the set (in data)" << endl;
			}
		Iso.induced_action_on_set_basic(Stab, data, 0 /*verbose_level*/);
		
			
		Iso.AA->group_order(Ago_induced[h]);
		}


	classify C_ago;

	C_ago.init(Ago_INT, Iso.Reps->count, FALSE, 0);
	cout << "Classification by ago:" << endl;
	C_ago.print(FALSE /*f_backwards*/);



	f << "\\chapter{Invariants}" << endl << endl;

	f << "Classification by automorphism group order: ";
	C_ago.print_naked_tex(f, FALSE /*f_backwards*/);
	f << "\\\\" << endl;

	f << "\\begin{center}" << endl;
	f << "\\begin{tabular}{|c|l|}" << endl;
	f << "\\hline" << endl;
	f << "Ago & Isom. Types \\\\" << endl;
	f << "\\hline" << endl;
	f << "\\hline" << endl;

	INT cnt, length, t, vv, *set;

	cnt = 0;
	for (u = C_ago.nb_types - 1; u >= 0; u--) {
		first = C_ago.type_first[u];
		length = C_ago.type_len[u];
		t = C_ago.data_sorted[first];

		f << t << " & ";

		set = NEW_INT(length);
		for (v = 0; v < length; v++, cnt++) {
			vv = first + v;
			i = C_ago.sorting_perm_inv[vv];
			set[v] = i;
			}

		INT_vec_heapsort(set, length);

		for (v = 0; v < length; v++, cnt++) {

			f << set[v];

			if (v < length - 1) {
				f << ",";
				if ((v + 1) % 10 == 0) {
					f << "\\\\" << endl;
					f << " & " << endl;
					}
				}
			}
		f << "\\\\" << endl;
		if (u > 0) {
			f << "\\hline" << endl;
			}
		FREE_INT(set);
		}
	f << "\\hline" << endl;
	f << "\\end{tabular}" << endl;
	f << "\\end{center}" << endl << endl;


	f << "\\clearpage" << endl << endl;

	f << "\\begin{center}" << endl;
	f << "\\begin{tabular}{|r|r|r|}" << endl;
	f << "\\hline" << endl;
	f << "Isom. Type & $|\\mbox{Aut}|$ & $|\\mbox{Aut}|$ (induced)\\\\" << endl;
	f << "\\hline" << endl;
	f << "\\hline" << endl;

	cnt = 0;
	for (u = 0; u < C_ago.nb_types; u ++) {
		first = C_ago.type_first[u];
		length = C_ago.type_len[u];
		t = C_ago.data_sorted[first];

		set = NEW_INT(length);
		for (v = 0; v < length; v++) {
			vv = first + v;
			i = C_ago.sorting_perm_inv[vv];
			set[v] = i;
			}

		INT_vec_heapsort(set, length);


		for (v = 0; v < length; v++) {
			vv = first + v;
			i = C_ago.sorting_perm_inv[vv];
			h = set[v];
			f << setw(3) << h << " & ";
			Ago[h].print_not_scientific(f);
			f << " & ";
			Ago_induced[h].print_not_scientific(f);
			f << "\\\\" << endl;
			cnt++;
			if ((cnt % 30) == 0) {
				f << "\\hline" << endl;
				f << "\\end{tabular}" << endl;
				f << "\\end{center}" << endl << endl;
				f << "\\begin{center}" << endl;
				f << "\\begin{tabular}{|r|r|r|}" << endl;
				f << "\\hline" << endl;
				f << "Isom. Type & $|\\mbox{Aut}|$ & $|\\mbox{Aut}|$ (induced)\\\\" << endl;
				f << "\\hline" << endl;
				f << "\\hline" << endl;
				}
			}
		FREE_INT(set);
		}

	f << "\\hline" << endl;
	f << "\\end{tabular}" << endl;
	f << "\\end{center}" << endl << endl;


	if (target_size == q + 2) {
		f << "\\chapter{The Hyperovals}" << endl << endl;
		}
	else {
		f << "\\chapter{The Arcs}" << endl << endl;
		}

	f << "\\clearpage" << endl << endl;


	for (h = 0; h < Iso.Reps->count; h++) {
		rep = Iso.Reps->rep[h];
		first = Iso.orbit_fst[rep];
		c = Iso.starter_number[first];
		id = Iso.orbit_perm[first];		
		Iso.load_solution(id, data);


		f << "\\section{Isomorphism type " << h << "}" << endl;
		f << "\\bigskip" << endl;


		if (Iso.Reps->stab[h]) {
			Iso.Reps->stab[h]->group_order(go);
			f << "Stabilizer has order $";
			go.print_not_scientific(f);
			f << "$.\\\\" << endl;
			}
		else {
			//cout << endl;
			}

		sims *Stab;
		
		Stab = Iso.Reps->stab[h];

		if (f_v) {
			cout << "arc_generator::report computing induced action on the set (in data)" << endl;
			}
		Iso.induced_action_on_set_basic(Stab, data, 0 /*verbose_level*/);
		
		longinteger_object go1;
			
		Iso.AA->group_order(go1);
		cout << "action " << Iso.AA->label << " computed, group order is " << go1 << endl;

		f << "Order of the group that is induced on the set is ";
		f << "$";
		go1.print_not_scientific(f);
		f << "$.\\\\" << endl;
		

		schreier Orb;
		//longinteger_object go2;
		
		Iso.AA->compute_all_point_orbits(Orb, Stab->gens, verbose_level - 2);
		f << "With " << Orb.nb_orbits << " orbits on the set.\\\\" << endl;

		classify C_ol;

		C_ol.init(Orb.orbit_len, Orb.nb_orbits, FALSE, 0);

		f << "Orbit lengths: ";
		//INT_vec_print(f, Orb.orbit_len, Orb.nb_orbits);
		C_ol.print_naked_tex(f, FALSE /*f_backwards*/);
		f << " \\\\" << endl;
	
		tt = (target_size + 3) / 4;

		f << "The points by ranks:\\\\" << endl;
		f << "\\begin{center}" << endl;

		for (u = 0; u < 4; u++) {
			f << "\\begin{tabular}[t]{|c|c|c|}" << endl;
			f << "\\hline" << endl;
			f << "$i$ & Rank & Unrank\\\\" << endl;
			f << "\\hline" << endl;
			for (i = 0; i < tt; i++) {
				v = u * tt + i;
				if (v < target_size) {
					INT vec[3];

					point_unrank(vec, data[v]);
					f << "$" << v << "$ & $" << data[v] << "$ & $";
					INT_vec_print(f, vec, 3);
					f << "$\\\\" << endl;
					}
				}
			f << "\\hline" << endl;
			f << "\\end{tabular}" << endl;
			}
		f << "\\end{center}" << endl; 


		report_stabilizer(Iso, f, h /* orbit */, 0 /* verbose_level */);


		report_decompositions(Iso, f, h /* orbit */, 
			data, verbose_level);

		}


	BYTE prefix[1000];
	BYTE label_of_structure_plural[1000];

	sprintf(prefix, "arcs_%ld_%ld", q, target_size);
	sprintf(label_of_structure_plural, "Arcs");
	isomorph_report_data_in_source_code_inside_tex(Iso, 
		prefix, label_of_structure_plural, f, 
		verbose_level);


	Iso.close_solution_database(verbose_level - 1);



	latex_foot(f);
	
	FREE_INT(Ago_INT);
	delete [] Ago;
	delete [] Ago_induced;
	}

	cout << "Written file " << fname << " of size " << file_size(fname) << endl;

}

void arc_generator::report_decompositions(isomorph &Iso, ofstream &f, INT orbit, 
	INT *data, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "arc_generator::report_decompositions" << endl;
		}
	incidence_structure *Inc;
	sims *Stab;
	strong_generators *gens;
	INT *Mtx;
	INT i, j, h;

	Inc = new incidence_structure;
	gens = new strong_generators;

	Stab = Iso.Reps->stab[orbit];
	gens->init_from_sims(Stab, 0 /* verbose_level */);
	//Stab->extract_strong_generators_in_order(*generators, tl, 0 /* verbose_level */);

	Mtx = NEW_INT(P2->N_points * P2->N_lines);
	INT_vec_zero(Mtx, P2->N_points * P2->N_lines);

	for (j = 0; j < P2->N_lines; j++) {
		for (h = 0; h < P2->k; h++) {
			i = P2->Lines[j * P2->k + h];
			Mtx[i * P2->N_lines + j] = 1;
			}
		}

	Inc->init_by_matrix(P2->N_points, P2->N_lines, Mtx, 0 /* verbose_level*/);
	

	partitionstack S;

	INT N;

	if (f_v) {
		cout << "arc_generator::report_decompositions allocating partitionstack" << endl;
		}
	N = Inc->nb_points() + Inc->nb_lines();
	
	S.allocate(N, 0);
	// split off the column class:
	S.subset_continguous(Inc->nb_points(), Inc->nb_lines());
	S.split_cell(0);
	S.split_cell_front_or_back(data, target_size, TRUE /* f_front */, 0 /* verbose_level*/);
				
	INT TDO_depth = N;
	INT TDO_ht;


	if (f_v) {
		cout << "arc_generator::report_decompositions before Inc->compute_TDO_safe" << endl;
		}
	Inc->compute_TDO_safe(S, TDO_depth, verbose_level - 3);
	TDO_ht = S.ht;


	if (S.ht < 50) {
		f << "The TDO decomposition is" << endl;
		Inc->get_and_print_column_tactical_decomposition_scheme_tex(f, TRUE /* f_enter_math */, S);
		}
	else {
		f << "The TDO decomposition is very large (with " << S.ht<< " classes).\\\\" << endl;
		}


	{
		schreier *Sch_points;
		schreier *Sch_lines;
		Sch_points = new schreier;
		Sch_points->init(A /*A_on_points*/);
		Sch_points->initialize_tables();
		Sch_points->init_generators(*gens->gens /* *generators */);
		Sch_points->compute_all_point_orbits(0 /*verbose_level - 2*/);
		
		if (f_v) {
			cout << "found " << Sch_points->nb_orbits << " orbits on points" << endl;
			}
		Sch_lines = new schreier;
		Sch_lines->init(A_on_lines);
		Sch_lines->initialize_tables();
		Sch_lines->init_generators(*gens->gens /* *generators */);
		Sch_lines->compute_all_point_orbits(0 /*verbose_level - 2*/);
		
		if (f_v) {
			cout << "found " << Sch_lines->nb_orbits << " orbits on lines" << endl;
			}
		S.split_by_orbit_partition(Sch_points->nb_orbits, 
			Sch_points->orbit_first, Sch_points->orbit_len, Sch_points->orbit,
			0 /* offset */, 
			verbose_level - 2);
		S.split_by_orbit_partition(Sch_lines->nb_orbits, 
			Sch_lines->orbit_first, Sch_lines->orbit_len, Sch_lines->orbit,
			Inc->nb_points() /* offset */, 
			verbose_level - 2);
		delete Sch_points;
		delete Sch_lines;
	}

	if (S.ht < 50) {
		f << "The TDA decomposition is" << endl;
		Inc->get_and_print_column_tactical_decomposition_scheme_tex(f, TRUE /* f_enter_math */, S);
		}
	else {
		f << "The TDA decomposition is very large (with " << S.ht<< " classes).\\\\" << endl;
		}

	FREE_INT(Mtx);
	delete gens;
	delete Inc;
}

void arc_generator::report_stabilizer(isomorph &Iso, ofstream &f, INT orbit, INT verbose_level)
{
	sims *Stab;
	longinteger_object go;
	INT i;

	Stab = Iso.Reps->stab[orbit];
	Stab->group_order(go);

	f << "The stabilizer of order $" << go << "$ is generated by:\\\\" << endl;
	for (i = 0; i < Stab->gens.len; i++) {
		
		INT *fp, n, ord;
		
		fp = NEW_INT(A->degree);
		n = A->find_fixed_points(Stab->gens.ith(i), fp, 0);
		//cout << "with " << n << " fixed points" << endl;
		FREE_INT(fp);

		ord = A->element_order(Stab->gens.ith(i));

		f << "$$ g_{" << i + 1 << "}=" << endl;
		A->element_print_latex(Stab->gens.ith(i), f);
		f << "$$" << endl << "of order $" << ord << "$ and with " << n << " fixed points." << endl;
		}
	f << endl << "\\bigskip" << endl;
}



// ##################################################################################################
// global functions
// ##################################################################################################

INT callback_arc_test(exact_cover *EC, INT *S, INT len, void *data, INT verbose_level)
{
	arc_generator *Gen = (arc_generator *) data;
	INT f_OK;
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "checking set ";
		print_set(cout, len, S);
		cout << endl;
		}
	f_OK = Gen->arc_test(S, len, verbose_level - 1);
	if (f_OK) {
		if (f_v) {
			cout << "accepted" << endl;
			}
		return TRUE;
		}
	else {
		if (f_v) {
			cout << "rejected" << endl;
			}
		return FALSE;
		}
}


INT check_arc(INT len, INT *S, void *data, INT verbose_level)
{
	arc_generator *Gen = (arc_generator *) data;
	INT f_OK;
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "checking set ";
		print_set(cout, len, S);
		cout << endl;
		}
	f_OK = Gen->check_arc(S, len, verbose_level - 1);
	if (f_OK) {
		if (f_v) {
			cout << "accepted" << endl;
			}
		return TRUE;
		}
	else {
		if (f_v) {
			cout << "rejected" << endl;
			}
		return FALSE;
		}
}

INT placebo_test_function(INT len, INT *S, void *data, INT verbose_level)
{
	//arc_generator *Gen = (arc_generator *) data;
	INT f_OK;
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "checking set ";
		print_set(cout, len, S);
		cout << endl;
		}
	f_OK = TRUE;
	if (f_OK) {
		if (f_v) {
			cout << "accepted" << endl;
			}
		return TRUE;
		}
	else {
		if (f_v) {
			cout << "rejected" << endl;
			}
		return FALSE;
		}
}


void arc_generator_early_test_function(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	void *data, INT verbose_level)
{
	arc_generator *Gen = (arc_generator *) data;
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "arc_generator_early_test_function for set ";
		print_set(cout, len, S);
		cout << endl;
		}
	Gen->early_test_func(S, len, 
		candidates, nb_candidates, 
		good_candidates, nb_good_candidates, 
		verbose_level - 2);
	if (f_v) {
		cout << "arc_generator_early_test_function done" << endl;
		}
}

void placebo_early_test_function(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	void *data, INT verbose_level)
{
	//arc_generator *Gen = (arc_generator *) data;
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "placebo_early_test_function for set ";
		print_set(cout, len, S);
		cout << endl;
		}

	INT_vec_copy(candidates, good_candidates, nb_candidates);
	nb_good_candidates = nb_candidates;

	if (f_v) {
		cout << "placebo_early_test_function done" << endl;
		}
}

void arc_generator_lifting_prepare_function_new(exact_cover *EC, INT starter_case, 
	INT *candidates, INT nb_candidates, strong_generators *Strong_gens, 
	diophant *&Dio, INT *&col_labels, 
	INT &f_ruled_out, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	arc_generator *Gen = (arc_generator *) EC->user_data;

	if (f_v) {
		cout << "arc_generator_lifting_prepare_function_new nb_candidates=" << nb_candidates << endl;
		}

	Gen->lifting_prepare_function_new(EC, starter_case, 
		candidates, nb_candidates, Strong_gens, 
		Dio, col_labels, f_ruled_out, 
		verbose_level - 1);


	if (f_v) {
		cout << "arc_generator_lifting_prepare_function_new nb_rows=" << Dio->m << " nb_cols=" << Dio->n << endl;
		}

	if (f_v) {
		cout << "arc_generator_lifting_prepare_function_new done" << endl;
		}
}

#if 0
void arc_generator_lifting_prepare_function(exact_cover *E, INT starter_case, 
	INT *candidates, INT nb_candidates, strong_generators *Strong_gens, 
	INT *&Inc, INT &nb_free_points, INT &nb_live_blocks2, INT *&live_blocks2, INT &nb_needed, 
	INT &f_has_RHS, INT *&RHS, INT &f_ruled_out, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	arc_generator *Gen = (arc_generator *) E->user_data;

	if (f_v) {
		cout << "arc_generator_lifting_prepare_function nb_candidates=" << nb_candidates << endl;
		}


	f_has_RHS = FALSE;
	RHS = NULL;
	Gen->lifting_prepare_function(E, starter_case, 
		candidates, nb_candidates, Strong_gens, 
		Inc, nb_free_points, nb_live_blocks2, live_blocks2, nb_needed, f_ruled_out, 
		verbose_level - 1);


	if (f_v) {
		cout << "arc_generator_lifting_prepare_function nb_free_points=" << nb_free_points << " nb_candidates=" << nb_candidates << endl;
		}

	if (f_v) {
		cout << "arc_generator_lifting_prepare_function done" << endl;
		}
}

void arc_generator_lifting_cleanup_function(exact_cover *E, INT starter_case, 
	INT *Inc, INT nb_free_points, INT nb_live_blocks2, INT *live_blocks2, 
	INT f_has_RHS, INT *RHS, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "arc_generator_lifting_cleanup_function" << endl;
		}
	FREE_INT(Inc);
	FREE_INT(live_blocks2);
	if (f_v) {
		cout << "arc_generator_lifting_cleanup_function done" << endl;
		}
}
#endif




void print_arc(ostream &ost, INT len, INT *S, void *data)
{
	arc_generator *Gen = (arc_generator *) data;
	
	Gen->print_set_in_affine_plane(ost, len, S);
}

void print_point(ostream &ost, INT pt, void *data)
{
	arc_generator *Gen = (arc_generator *) data;
	INT v[3];
	
	PG_element_unrank_modified(*Gen->F, v, 1 /* stride */, 3 /* len */, pt);
	ost << "(" << v[0] << "," << v[1] << "," << v[2] << ")" << endl;
}

void callback_arc_report(isomorph *Iso, void *data, INT verbose_level)
{
	arc_generator *Gen = (arc_generator *) data;
	
	Gen->report(*Iso, verbose_level);
}



