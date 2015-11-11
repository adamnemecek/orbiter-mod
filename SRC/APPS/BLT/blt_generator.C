// blt_generator.C
// 
// Anton Betten
// moved here from blt.C 5/24/09
//
//
// check_condition() is the test function for Snakes and Ladders. 
// It calls blt_generator::check_condition(), 
// which tests whether a set of points is admissible. 
// This is based on collinearity test and BLT-test.
// Collinearity test only tests the 
// last point in the set, since the subset of all but the last 
// has been checked previously.
// BLT-test tests the whole set, so there is a little bit of waste here.
//
//

#include "orbiter.h"
#include "discreta.h"

#include "blt.h"

void blt_generator::read_arguments(int argc, const char **argv)
{
	INT i;
	
#if 0
	for (i = 1; i < argc; i++) {
		cout << argv[i] << endl;
		}
#endif
	gen->read_arguments(argc, argv, 0);
	
	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-schreier") == 0) {
			f_override_schreier_depth = TRUE;
			override_schreier_depth = atoi(argv[++i]);
			cout << "-schreier " << override_schreier_depth << endl;
			}
		else if (strcmp(argv[i], "-n") == 0) {
			f_override_n = TRUE;
			override_n = atoi(argv[++i]);
			cout << "-override_n " << override_n << endl;
			}
		else if (strcmp(argv[i], "-epsilon") == 0) {
			f_override_epsilon = TRUE;
			override_epsilon = atoi(argv[++i]);
			cout << "-override_epsilon " << override_epsilon << endl;
			}
		else if (strcmp(argv[i], "-BLT") == 0) {
			f_BLT = TRUE;
			cout << "-BLT " << endl;
			}
		else if (strcmp(argv[i], "-ovoid") == 0) {
			f_ovoid = TRUE;
			cout << "-ovoid " << endl;
			}
		else if (strcmp(argv[i], "-semilinear") == 0) {
			f_semilinear = TRUE;
			cout << "-semilinear" << endl;
			}
		}
	if (!f_BLT && !f_ovoid) {
		cout << "please use either -BLT or -ovoid" << endl;
		exit(1);
		}
}

blt_generator::blt_generator()
{
	override_poly = NULL;
	f_semilinear = FALSE;
	gen = NULL;
	F = NULL;
	A = NULL;
	O = NULL;
	f_BLT = FALSE;
	f_ovoid = FALSE;
	f_semilinear = FALSE;
	f_orthogonal_allocated = FALSE;
	nb_sol = 0;
	f_override_schreier_depth = FALSE;
	f_override_n = FALSE;
	override_n = 0;
	f_override_epsilon = FALSE;
	override_epsilon = 0;
}

blt_generator::~blt_generator()
{
	INT f_v = TRUE;

	if (f_v) {
		cout << "blt_generator::~blt_generator before A" << endl;
		}
	if (A) {
		delete A;
		A = NULL;
		}
	if (f_v) {
		cout << "blt_generator::~blt_generator before gen" << endl;
		}
	if (gen) {
		delete gen;
		gen = NULL;
		}
	if (f_orthogonal_allocated) {
		if (f_v) {
			cout << "blt_generator::~blt_generator before O" << endl;
			}
		if (O) {
			delete O;
			}
		f_orthogonal_allocated = FALSE;
		O = NULL;
		}
	if (f_v) {
		cout << "blt_generator::~blt_generator done" << endl;
		}
	
}



void blt_generator::init_basic(INT q, 
	const BYTE *starter_directory_name, 
	const BYTE *starter_prefix, 
	int argc, const char **argv, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);

	if (f_v) {
		cout << "blt_generator::init_basic" << endl;
		}

	F = new finite_field;
	A = new action;
	gen = new generator;
	
	if (f_vv) {
		cout << "blt_generator::init_basic before read_arguments" << endl;
		}

	read_arguments(argc, argv);
	blt_generator::q = q;
	blt_generator::starter_directory_name = starter_directory_name;
	blt_generator::starter_prefix = starter_prefix;

	target_size = q + 1;
		

	if (f_vv) {
		cout << "blt_generator::init_basic q=" << q << " target_size = " << target_size << endl;
		}
	
	override_poly = override_polynomial_subfield(q);
	
	if (f_vv && override_poly) {
		cout << "blt_generator::init_basic init override_poly=" << override_poly << endl;
		}
	F->init_override_polynomial(q, override_poly, 0);


	n = 0;
	epsilon = 0;
	
	sprintf(gen->fname_base, "%s%s", starter_directory_name, starter_prefix);
	
	if (f_BLT) {
		epsilon = 0;
		n = 5;
		}
	else if (f_ovoid) {
		if (f_override_n) {
			n = override_n;
			if (f_vv) {
				cout << "blt_generator::init_basic override value of n=" << n << endl;
				}
			}
		if (f_override_epsilon) {
			epsilon = override_epsilon;
			if (f_vv) {
				cout << "blt_generator::init_basic override value of epsilon=" << epsilon << endl;
				}
			}
		//sprintf(gen->fname_base, "%sovoid_O%s_%ld_%ld", gen->prefix, plus_minus_letter(epsilon), n, q);
		}
	else {
		cout << "neither f_BLT nor f_ovoid is TRUE" << endl;
		exit(1);
		}
	
	f_semilinear = TRUE;
	if (is_prime(q)) {
		f_semilinear = FALSE;
		}
	if (f_vv) {
		cout << "blt_generator::init_basic f_semilinear=" << f_semilinear << endl;
		}
	if (f_v) {
		cout << "blt_generator::init_basic finished" << endl;
		}
}

void blt_generator::init_group(INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT f_basis = TRUE;

	if (f_v) {
		cout << "blt_generator::init_group" << endl;
		}
	if (f_vv) {
		cout << "blt_generator::init_group epsilon=" << epsilon << endl;
		cout << "blt_generator::init_group n=" << n << endl;
		cout << "blt_generator::init_group q=" << q << endl;
		cout << "blt_generator::init_group f_semilinear=" << f_semilinear << endl;
		}
	if (f_vv) {
		cout << "blt_generator::init_group before A->init_orthogonal_group" << endl;
		}
	A->init_orthogonal_group(epsilon, n, F, 
		TRUE /* f_on_points */, 
		FALSE /* f_on_lines */, 
		FALSE /* f_on_points_and_lines */, 
		f_semilinear, f_basis, verbose_level - 1);
	if (f_vv) {
		cout << "blt_generator::init_group after A->init_orthogonal_group" << endl;
		}
	
	if (f_vv) {
		cout << "blt_generator::init_group computing lex least base" << endl;
		}
	A->lex_least_base_in_place(0 /*verbose_level - 2*/);
	if (f_vv) {
		cout << "blt_generator::init_group computing lex least base done" << endl;
		cout << "blt_generator::init_group base: ";
		INT_vec_print(cout, A->base, A->base_len);
		cout << endl;
		}
	
	action_on_orthogonal *AO;

	AO = A->G.AO;
	O = AO->O;

	if (f_v) {
		cout << "blt_generator::init_group degree = " << A->degree << endl;
		}
		
	//init_orthogonal_hash(verbose_level);

	if (A->degree < 200) {
		if (f_v) {
			cout << "blt_generator::init_group before test_Orthogonal" << endl;
			}
		test_Orthogonal(epsilon, n - 1, q);
		}
	//A->Sims->print_all_group_elements();

	if (FALSE) {
		cout << "blt_generator::init_group before A->Sims->print_all_transversal_elements" << endl;
		A->Sims->print_all_transversal_elements();
		cout << "blt_generator::init_group after A->Sims->print_all_transversal_elements" << endl;
		}


	if (FALSE /*f_vv*/) {
		O->F->print(FALSE);
		}


	
	if (f_v) {
		cout << "blt_generator::init_group finished" << endl;
		}
}

void blt_generator::init_orthogonal(INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	
	if (f_v) {
		cout << "blt_generator::init_orthogonal" << endl;
		}
	if (f_vv) {
		cout << "epsilon=" << epsilon << endl;
		cout << "n=" << n << endl;
		cout << "q=" << q << endl;
		cout << "override_poly=" << override_poly << endl;
		cout << "f_semilinear=" << f_semilinear << endl;
		}

	f_orthogonal_allocated = TRUE;
	O = new orthogonal;
	O->init(epsilon, n, F, verbose_level - 3);
	if (f_vv) {
		cout << "created O^" << plus_minus_string(epsilon) << "(" << n << "," << q << ") with " 
			<< O->nb_points << " points and " << O->nb_lines << " lines" << endl << endl;
		}

	init_orthogonal_hash(verbose_level);

	if (f_vv) {
		O->F->print(FALSE);
		}

	if (f_v) {
		cout << "blt_generator::init_orthogonal finished" << endl;
		}
}

void blt_generator::init_orthogonal_hash(INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "blt_generator::init_orthogonal_hash" << endl;
		}

	init_hash_table_parabolic(*O->F, 4, 0/*verbose_level*/);

	if (f_v) {
		cout << "blt_generator::init_orthogonal finished" << endl;
		}
}

void blt_generator::init2(INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	//INT f_vvv = (verbose_level >= 3);

	if (f_v) {
		cout << "blt_generator::init2" << endl;
		}
	if (gen->f_max_depth) {
		gen->depth = gen->max_depth;
		}
	else {
		gen->depth = q + 1;
		}
	
	if (f_v) {
		cout << "blt_generator::init2 depth = " << gen->depth << endl;
		}

	
	gen->init(A, A, A->Strong_gens, 
		gen->depth /* sz */, verbose_level);
	
#if 0
	// not needed since we have an early_test_func:
	gen->init_check_func(::check_conditions, 
		(void *)this /* candidate_check_data */);
#endif

	// we have an early test function:

	gen->init_early_test_func(
		early_test_func_callback, 
		this,  
		verbose_level);

	// We also have an incremental check function. 
	// This is only used by the clique finder:
	gen->init_incremental_check_func(
		check_function_incremental_callback, 
		this /* candidate_check_data */);


	gen->f_print_function = TRUE;
	gen->print_function = print_set;
	gen->print_function_data = (void *) this;
	
	
	INT nb_oracle_nodes = ONE_MILLION;
	
	if (f_vv) {
		cout << "blt_generator::init2 calling init_oracle with " << nb_oracle_nodes << " nodes" << endl;
		}
	
	gen->init_oracle(nb_oracle_nodes, verbose_level - 1);

	if (f_vv) {
		cout << "blt_generator::init2 after init_root_node" << endl;
		}
	
	//cout << "verbose_level = " << verbose_level << endl;
	//cout << "verbose_level_group_theory = " << verbose_level_group_theory << endl;
	
	gen->root[0].init_root_node(gen, 0/*verbose_level - 2*/);
	if (f_v) {
		cout << "blt_generator::init2 done" << endl;
		}
}

INT blt_generator::pair_test(INT a, INT x, INT y, INT verbose_level)
// We assume that a is an element of a set S of size at least two such that 
// S \cup \{ x \} is BLT and 
// S \cup \{ y \} is BLT.
// In order to test if S \cup \{ x, y \} is BLT, we only need to test 
// the triple \{ x,y,a\}
{
	INT v1[5], v2[5], v3[5];
	INT f12, f13, f23;
	INT d;

	O->unrank_point(v1, 1, a, 0);
	O->unrank_point(v2, 1, x, 0);
	O->unrank_point(v3, 1, y, 0);
	f12 = O->evaluate_bilinear_form(v1, v2, 1);
	f13 = O->evaluate_bilinear_form(v1, v3, 1);
	f23 = O->evaluate_bilinear_form(v2, v3, 1);
	d = O->F->product3(f12, f13, f23);
	if (d == 0) {
		return FALSE;
		}
	if (O->f_is_minus_square[d]) {
		return FALSE;
		}
	else {
		return TRUE;
		}
	
}

INT blt_generator::check_conditions(INT len, INT *S, INT verbose_level)
{
	INT f_OK = TRUE;
	INT f_BLT_test = FALSE;
	INT f_collinearity_test = FALSE;
	//INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	
	//f_v = TRUE;
	//f_vv = TRUE;
	
	if (f_vv) {
		cout << "checking set ";
		print_set(cout, len, S);
		}
	if (!collinearity_test(S, len, verbose_level)) {
		f_OK = FALSE;
		f_collinearity_test = TRUE;
		}
	if (f_BLT) {
		if (!O->BLT_test(len, S, verbose_level)) {
			f_OK = FALSE;
			f_BLT_test = TRUE;
			}
		}


	if (f_OK) {
		if (f_vv) {
			cout << "OK" << endl;
			}
		return TRUE;
		}
	else {
		if (f_vv) {
			cout << "not OK because of ";
			if (f_BLT_test) {
				cout << "BLT test";
				}
			if (f_collinearity_test) {
				cout << "collinearity test";
				}
			cout << endl;
			}
		return FALSE;
		}
}

INT blt_generator::collinearity_test(INT *line, INT len, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, x, y;
	INT f_OK = TRUE;
	INT fxy;
	
	if (f_v) {
		cout << "\ncollinearity test for" << endl;
		for (i = 0; i < len; i++) {
			O->unrank_point(O->v1, 1, line[i], 0);
			INT_vec_print(cout, O->v1, n);
			cout << endl;
			}
		}
	y = line[len - 1];
	O->unrank_point(O->v1, 1, y, 0);
	
	for (i = 0; i < len - 1; i++) {
		x = line[i];
		O->unrank_point(O->v2, 1, x, 0);
		fxy = O->evaluate_bilinear_form(O->v1, O->v2, 1);
		
		if (fxy == 0) {
			f_OK = FALSE;
			if (f_v) {
				cout << "not OK; ";
				cout << "{x,y}={" << x << "," << y << "} are collinear" << endl;
				INT_vec_print(cout, O->v1, n);
				cout << endl;
				INT_vec_print(cout, O->v2, n);
				cout << endl;
				cout << "fxy=" << fxy << endl;
				}
			break;
			}
		}
	
	if (f_v) {
		if (!f_OK) {
			cout << "collinearity test fails" << endl;
			}
		}
	return f_OK;
}

void blt_generator::print(ostream &ost, INT *S, INT len)
{
	INT i;
	
	for (i = 0; i < len; i++) {
		O->unrank_point(O->v1, 1, S[i], 0);
		INT_vec_print(ost, O->v1, n);
		ost << endl;
		}
}






void blt_generator::create_system(INT level, 
	INT orbit_at_level, 
	INT level_of_candidates_file, 
	const BYTE *output_prefix, INT width, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);


	if (f_v) {
		cout << "blt_generator::create_system" << endl;
		}


	BYTE prefix[1000];
	orbit_rep *R;
	INT starter_sz = level;



	INT max_starter;
	INT nb;

	INT *lines_on_pt;
	INT *Perp;
	INT i, j, a;
	classify C;
	INT f, h;
	INT fst, len, b, pt;
	INT *free_pts;
	INT nb_free_pts;
	INT *free_pt_idx;
	INT *Incma;
	INT *pts_on_line;

	INT *Pts1;
	INT *Pts2;

	sprintf(prefix, "BLT_%ld", q);

	R = new orbit_rep;
	R->init_from_file(A, prefix, 
		level, orbit_at_level, level_of_candidates_file, 
		early_test_func_callback, 
		this /* early_test_func_callback_data */, 
		verbose_level - 1
		);

	nb = q + 1 - level;


	if (f_vv) {
		cout << "blt_generator::create_system Read starter " << orbit_at_level << " : ";
		INT_vec_print(cout, R->rep, starter_sz);
		cout << endl;
		}

	max_starter = R->rep[level - 1];

	if (f_vv) {
		cout << "blt_generator::create_system max_starter=" << max_starter << endl;
		cout << "blt_generator::create_system Group order=" << R->stab_go << endl;
		cout << "blt_generator::create_system nb_candidates=" << R->nb_candidates 
			<< " at level " << starter_sz << endl;
		}


	INT nb_candidates2;
	
	if (f_vv) {
		cout << "blt_generator::create_system Before lexorder_test" << endl;
		}
	A->lexorder_test(R->candidates, R->nb_candidates, nb_candidates2, 
		R->Strong_gens->gens, max_starter, verbose_level - 2);
	if (f_vv) {
		cout << "blt_generator::create_system After lexorder_test nb_candidates=" << nb_candidates2 
			<< " eliminated " << R->nb_candidates - nb_candidates2 << " candidates" << endl;
		}
	R->nb_candidates = nb_candidates2;


	if (R->nb_candidates < nb) {
		cout << "blt_generator::create_system nb_candidates < nb, the case is eliminated" << endl;
		delete R;
		return;
		}
	

	lines_on_pt = NEW_INT(starter_sz * (q + 1));
	for (i = 0; i < starter_sz; i++) {
		O->lines_on_point_by_line_rank(R->rep[i], lines_on_pt + i * (q + 1), 0 /* verbose_level */);
		}

	if (f_vv) {
		cout << "blt_generator::create_system Lines on partial BLT set:" << endl;
		INT_matrix_print(lines_on_pt, starter_sz, q + 1);
		}

	Perp = NEW_INT(starter_sz * (q + 1) * (q + 1));
	for (i = 0; i < starter_sz; i++) {
		for (j = 0; j < q + 1; j++) {
			a = lines_on_pt[i * (q + 1) + j];
			O->points_on_line_by_line_rank(a, Perp + i * (q + 1) * (q + 1) + j * (q + 1), 0 /* verbose_level */);
			}
		}
	if (f_vv) {
		cout << "blt_generator::create_system Perp:" << endl;
		INT_matrix_print(Perp, starter_sz * (q + 1), q + 1);
		}
	

	C.init(Perp, starter_sz * (q + 1) * (q + 1), TRUE, 0);

	C.print(FALSE /* f_reverse */);


	f = C.second_type_first[0];
	nb_free_pts = C.second_type_len[0];
	if (f_v) {
		cout << "blt_generator::create_system nb_free_pts=" << nb_free_pts << endl;
		}
	free_pts = NEW_INT(nb_free_pts);
	free_pt_idx = NEW_INT(O->nb_points);
	for (h = 0; h < O->nb_points; h++) {
		free_pt_idx[h] = -1;
		}
	
	for (h = 0; h < nb_free_pts; h++) {
		b = C.second_sorting_perm_inv[f + h];
		fst = C.type_first[b];
		len = C.type_len[b];
		if (len != 1) {
			cout << "len != 1" << endl;
			exit(1);
			}
		pt = C.data_sorted[fst];
		//cout << "h=" << h << " b=" << b << " len=" << len << " pt=" << pt << endl;
		free_pts[h] = pt;
		free_pt_idx[pt] = h;
		}

	if (f_v) {
		cout << "blt_generator::create_system There are " << nb_free_pts << " free points" << endl;
		}


	Pts1 = NEW_INT(nb_free_pts * 5);
	Pts2 = NEW_INT(R->nb_candidates * 5);
	for (i = 0; i < nb_free_pts; i++) {
		O->unrank_point(Pts1 + i * 5, 1, free_pts[i], 0 /*verbose_level - 1*/);
		}
	for (i = 0; i < R->nb_candidates; i++) {
		O->unrank_point(Pts2 + i * 5, 1, R->candidates[i], 0 /*verbose_level - 1*/);
		}

	Incma = NEW_INT(nb_free_pts * R->nb_candidates);
	pts_on_line = NEW_INT(q + 1);
	INT_vec_zero(Incma, nb_free_pts * R->nb_candidates);
	for (i = 0; i < nb_free_pts; i++) {
		for (j = 0; j < R->nb_candidates; j++) {
			a = O->evaluate_bilinear_form(Pts1 + i * 5, Pts2 + j * 5, 1);
			if (a == 0) {
				Incma[i * R->nb_candidates + j] = 1;
				}
			}
		}

	if (f_v) {
		cout << "blt_generator::create_system Incma computed" << endl;
		}
	//INT_matrix_print(Incma, nb_free_pts, R->nb_candidates);


	write_problem_to_file_wassermann(
		output_prefix, 
		Incma, 
		nb_free_pts, R->nb_candidates, 
		starter_sz, orbit_at_level, 
		width);

	delete R;

	FREE_INT(lines_on_pt);
	FREE_INT(Perp);
	FREE_INT(Incma);
	FREE_INT(pts_on_line);
	FREE_INT(Pts1);
	FREE_INT(Pts2);

	if (f_v) {
		cout << "blt_generator::create_system done" << endl;
		}
}

void blt_generator::create_graphs(INT starter_depth, 
	INT orbit_at_level_r, INT orbit_at_level_m, 
	INT level_of_candidates_file, 
	const BYTE *output_prefix, 
	INT f_lexorder_test, INT f_eliminate_graphs_if_possible, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT f_v3 = (verbose_level >= 3);


	if (f_v) {
		cout << "blt_generator::create_graphs" << endl;
		cout << "blt_generator::create_graphs starter_depth = " << starter_depth << endl;
		cout << "blt_generator::create_graphs f_lexorder_test=" << f_lexorder_test << endl;
		}


	//f_memory_debug = TRUE;


	BYTE fname[1000];
	BYTE fname_list_of_cases[1000];
	BYTE graph_fname_base[1000];
	INT orbit;
	INT nb_orbits;
	//INT width;
	INT *list_of_cases;
	INT nb_of_cases;
	//INT nb_vertices;




	sprintf(fname, "%s%s_lvl_%ld", starter_directory_name, starter_prefix, starter_depth);
	sprintf(fname_list_of_cases, "%slist_of_cases_%s_%ld_%ld_%ld.txt", output_prefix, starter_prefix, starter_depth, orbit_at_level_r, orbit_at_level_m);

	nb_orbits = count_number_of_orbits_in_file(fname, 0);
	if (f_v) {
		cout << "blt_generator::create_graphs There are " << nb_orbits << " starters" << endl;
		}
	if (nb_orbits < 0) {
		cout << "Something is wrong, nb_orbits is negative" << endl;
		exit(1);
		}

#if 0
	//width = log10(nb_orbits) + 1;
	if (f_v) {
		cout << "blt_generator::create_graphs width=" << width << endl;
		}
#endif

	nb_of_cases = 0;
	list_of_cases = NEW_INT(nb_orbits);
	for (orbit = 0; orbit < nb_orbits; orbit++) {
		if ((orbit % orbit_at_level_m) != orbit_at_level_r) {
			continue;
			}
		if (f_v3) {
			cout << "blt_generator::create_graphs creating graph associated with orbit " << orbit << " / " << nb_orbits << ":" << endl;
			}

		
		colored_graph *CG;
		INT nb_vertices = -1;


		if (create_graph(starter_depth, orbit, level_of_candidates_file, 
			output_prefix, 
			f_lexorder_test, f_eliminate_graphs_if_possible, 
			nb_vertices, graph_fname_base,
			CG,  
			verbose_level - 2)) {
			list_of_cases[nb_of_cases++] = orbit;

			BYTE fname[1000];

			sprintf(fname, "%s%s.bin", output_prefix, CG->fname_base);
			CG->save(fname, verbose_level - 2);
			
			nb_vertices = CG->nb_points;
			//delete CG;
			}

		if (CG) {
			delete CG;
			}
		if (f_vv) {
			if (nb_vertices >= 0) {
				cout << "blt_generator::create_graphs creating graph associated with orbit " << orbit << " / " << nb_orbits << " with " << nb_vertices << " vertices created" << endl;
				}
			else {
				cout << "blt_generator::create_graphs creating graph associated with orbit " << orbit << " / " << nb_orbits << " is ruled out" << endl;
				}
			}
		}

	write_set_to_file(fname_list_of_cases, list_of_cases, nb_of_cases, 0 /*verbose_level */);
	if (f_v) {
		cout << "blt_generator::create_graphs Written file " << fname_list_of_cases << " of size " << file_size(fname_list_of_cases) << endl;
		}

	FREE_INT(list_of_cases);

	//registry_dump_sorted();
}


INT blt_generator::create_graph(INT level, 
	INT orbit_at_level, INT level_of_candidates_file, 
	const BYTE *output_prefix, 
	INT f_lexorder_test, INT f_eliminate_graphs_if_possible, 
	INT &nb_vertices, BYTE *graph_fname_base,
	colored_graph *&CG,  
	INT verbose_level)
// returns TRUE if a graph was written, FALSE otherwise
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT f_v3 = (verbose_level >= 3);


	if (f_v) {
		cout << "blt_generator::create_graph" << endl;
		cout << "blt_generator::create_graph level = " << level << endl;
		cout << "blt_generator::create_graph f_lexorder_test=" << f_lexorder_test << endl;
		}

	CG = NULL;
	
	INT ret;

	BYTE prefix[1000];
	orbit_rep *R;

	//INT i;


	INT starter_sz = level;
	INT max_starter;
	INT nb;

	nb_vertices = 0;

	sprintf(prefix, "%s%s", starter_directory_name, starter_prefix);
	

	R = new orbit_rep;
	R->init_from_file(A, prefix, 
		level, orbit_at_level, level_of_candidates_file, 
		early_test_func_callback, 
		this /* early_test_func_callback_data */, 
		verbose_level - 1
		);
	nb = q + 1 - starter_sz;


	if (f_vv) {
		cout << "blt_generator::create_graph Case " << orbit_at_level << " / " << R->nb_cases << " Read starter : ";
		INT_vec_print(cout, R->rep, starter_sz);
		cout << endl;
		}

	max_starter = R->rep[starter_sz - 1];

	if (f_vv) {
		cout << "blt_generator::create_graph Case " << orbit_at_level << " / " << R->nb_cases << " max_starter=" << max_starter << endl;
		cout << "blt_generator::create_graph Case " << orbit_at_level << " / " << R->nb_cases << " Group order=" << R->stab_go << endl;
		cout << "blt_generator::create_graph Case " << orbit_at_level << " / " << R->nb_cases << " nb_candidates=" << R->nb_candidates << " at level " << level << endl;
		}



	if (f_lexorder_test) {
		INT nb_candidates2;
	
		if (f_v3) {
			cout << "blt_generator::create_graph Case " << orbit_at_level << " / " << R->nb_cases << " Before lexorder_test" << endl;
			}
		A->lexorder_test(R->candidates, R->nb_candidates, nb_candidates2, 
			R->Strong_gens->gens, max_starter, verbose_level - 3);
		if (f_vv) {
			cout << "blt_generator::create_graph After lexorder_test nb_candidates=" << nb_candidates2 << " eliminated " << R->nb_candidates - nb_candidates2 << " candidates" << endl;
			}
		R->nb_candidates = nb_candidates2;
		}


	// we must do this. 
	// For instance, what of we have no points left, then the minimal color stuff break down.
	//if (f_eliminate_graphs_if_possible) {
		if (R->nb_candidates < nb) {
			if (f_v) {
				cout << "blt_generator::create_graph Case " << orbit_at_level << " / " << R->nb_cases << " nb_candidates < nb, the case is eliminated" << endl;
				}
			delete R;
			return FALSE;
			}
		//}


	nb_vertices = R->nb_candidates;


	INT *point_color;
	INT nb_colors;

	INT *lines_on_pt;
	
	lines_on_pt = NEW_INT(1 /*starter_sz*/ * (q + 1));
	O->lines_on_point_by_line_rank(R->rep[0], lines_on_pt, 0 /* verbose_level */);

	if (f_v3) {
		cout << "Case " << orbit_at_level << " Lines on partial BLT set:" << endl;
		INT_matrix_print(lines_on_pt, 1 /*starter_sz*/, q + 1);
		}

	INT special_line;

	special_line = lines_on_pt[0];

	compute_colors(orbit_at_level, 
		R->rep, starter_sz, 
		special_line, 
		R->candidates, R->nb_candidates, 
		point_color, nb_colors, 
		verbose_level - 3);


	classify C;

	C.init(point_color, R->nb_candidates, FALSE, 0);
	if (f_v3) {
		cout << "blt_generator::create_graph Case " << orbit_at_level << " / " << R->nb_cases << " point colors (1st classification): ";
		C.print(FALSE /* f_reverse */);
		cout << endl;
		}


	classify C2;

	C2.init(point_color, R->nb_candidates, TRUE, 0);
	if (f_vv) {
		cout << "blt_generator::create_graph Case " << orbit_at_level << " / " << R->nb_cases << " point colors (2nd classification): ";
		C2.print(FALSE /* f_reverse */);
		cout << endl;
		}



	INT f, l, idx;

	f = C2.second_type_first[0];
	l = C2.second_type_len[0];
	idx = C2.second_sorting_perm_inv[f + 0];
#if 0
	if (C.type_len[idx] != minimal_type_multiplicity) {
		cout << "idx != minimal_type" << endl;
		cout << "idx=" << idx << endl;
		cout << "minimal_type=" << minimal_type << endl;
		cout << "C.type_len[idx]=" << C.type_len[idx] << endl;
		cout << "minimal_type_multiplicity=" << minimal_type_multiplicity << endl;
		exit(1);
		}
#endif
	INT minimal_type, minimal_type_multiplicity;
	
	minimal_type = idx;
	minimal_type_multiplicity = C2.type_len[idx];

	if (f_vv) {
		cout << "blt_generator::create_graph Case " << orbit_at_level << " / " << R->nb_cases << " minimal type is " << minimal_type << endl;
		cout << "blt_generator::create_graph Case " << orbit_at_level << " / " << R->nb_cases << " minimal_type_multiplicity " << minimal_type_multiplicity << endl;
		}

	if (f_eliminate_graphs_if_possible) {
		if (minimal_type_multiplicity == 0) {
			cout << "blt_generator::create_graph Case " << orbit_at_level << " / " << R->nb_cases << " Color class " << minimal_type << " is empty, the case is eliminated" << endl;
			ret = FALSE;
			goto finish;
			}
		}



	if (f_vv) {
		cout << "blt_generator::create_graph Case " << orbit_at_level << " / " << R->nb_cases << " Computing adjacency list, nb_points=" << R->nb_candidates << endl;
		}

	UBYTE *bitvector_adjacency;
	INT bitvector_length_in_bits;
	INT bitvector_length;

	compute_adjacency_list_fast(R->rep[0], 
		R->candidates, R->nb_candidates, point_color, 
		bitvector_adjacency, bitvector_length_in_bits, bitvector_length, 
		verbose_level - 2);

	if (f_vv) {
		cout << "blt_generator::create_graph Case " << orbit_at_level << " / " << R->nb_cases << " Computing adjacency list done" << endl;
		cout << "blt_generator::create_graph Case " << orbit_at_level << " / " << R->nb_cases << " bitvector_length=" << bitvector_length << endl;
		}


	if (f_v) {
		cout << "blt_generator::create_graph creating colored_graph" << endl;
		}

	//colored_graph *CG;

	CG = new colored_graph;

	CG->init(R->nb_candidates /* nb_points */, nb_colors, 
		point_color, bitvector_adjacency, verbose_level - 2);
		// the adjacency becomes part of the colored_graph object
	
	INT i;
	for (i = 0; i < R->nb_candidates; i++) {
		CG->points[i] = R->candidates[i];
		}
	CG->init_user_data(R->rep, level, verbose_level - 2);
	sprintf(CG->fname_base, "graph_BLT_%ld_%ld_%ld", q, level, orbit_at_level);
		

	if (f_v) {
		cout << "blt_generator::create_graph colored_graph created" << endl;
		}

#if 0
	sprintf(graph_fname_base, "graph_BLT_%ld_%ld_%ld", q, level, orbit_at_level);
	BYTE fname[1000];

	sprintf(fname, "%s%s.bin", output_prefix, graph_fname_base);

	save_colored_graph(fname, R->nb_candidates, nb_colors, 
		R->candidates, point_color, 
		R->rep, level, 
		bitvector_adjacency, bitvector_length,
		verbose_level - 1);
		// GALOIS/galois_global.C


	if (f_vv) {
		cout << "blt_generator::create_graph Case " << orbit_at_level << " Written file " << fname << " of size " << file_size(fname) << endl;
		}

#endif

	//FREE_UBYTE(bitvector_adjacency);
	FREE_INT(lines_on_pt);
	FREE_INT(point_color);


	ret = TRUE;

finish:
	delete R;
	return ret;
}



void blt_generator::compute_adjacency_list_fast(INT first_point_of_starter, 
	INT *points, INT nb_points, INT *point_color, 
	UBYTE *&bitvector_adjacency, INT &bitvector_length_in_bits, INT &bitvector_length, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT L;
	INT i, j, k, c1, c2;
	INT *Pts;
	INT *form_value;
	INT v1[5];
	INT m[5];
	INT f12, f13, f23, d;
	UINT cnt;
	INT two;
	finite_field *F;
	INT *Pi, *Pj;

	if (f_v) {
		cout << "blt_generator::compute_adjacency_list_fast" << endl;
		}
	L = (nb_points * (nb_points - 1)) >> 1;

	bitvector_length_in_bits = L;
	bitvector_length = (L + 7) >> 3;
	bitvector_adjacency = NEW_UBYTE(bitvector_length);
	for (i = 0; i < bitvector_length; i++) {
		bitvector_adjacency[i] = 0;
		}
	
	Pts = NEW_INT(nb_points * 5);
	form_value = NEW_INT(nb_points);
	O->unrank_point(v1, 1, first_point_of_starter, 0);
	if (f_v) {
		cout << "blt_generator::compute_adjacency_list_fast unranking points" << endl;
		}
	for (i = 0; i < nb_points; i++) {
		O->unrank_point(Pts + i * 5, 1, points[i], 0);
		form_value[i] = O->evaluate_bilinear_form(v1, Pts + i * 5, 1);
		}

	if (f_v) {
		cout << "blt_generator::compute_adjacency_list_fast computing adjacencies" << endl;
		}

	cnt = 0;
	F = O->F;
	two = F->add(1, 1);
	
	for (i = 0; i < nb_points; i++) {
		f12 = form_value[i];
		c1 = point_color[i];
		Pi = Pts + i * 5;
		m[0] = F->mult(Pi[0], two);
		m[1] = Pi[2];
		m[2] = Pi[1];
		m[3] = Pi[4];
		m[4] = Pi[3];
		
		for (j = i + 1; j < nb_points; j++, cnt++) {
			k = ij2k(i, j, nb_points);
		
			if ((cnt & ((1 << 25) - 1)) == 0 && cnt) {
				cout << "adjacency " << cnt << " / " << L << endl;
				}
			c2 = point_color[j];
			if (c1 == c2) {
				bitvector_m_ii(bitvector_adjacency, k, 0);
				continue;
				}
			f13 = form_value[j];
			Pj = Pts + j * 5;
			f23 = F->add5(
				F->mult(m[0], Pj[0]), 
				F->mult(m[1], Pj[1]), 
				F->mult(m[2], Pj[2]), 
				F->mult(m[3], Pj[3]), 
				F->mult(m[4], Pj[4])
				);
			d = F->product3(f12, f13, f23);
			if (d == 0) {
				bitvector_m_ii(bitvector_adjacency, k, 0);
				}
			else {
				if (O->f_is_minus_square[d]) {
					bitvector_m_ii(bitvector_adjacency, k, 0);
					}
				else {
					bitvector_m_ii(bitvector_adjacency, k, 1);
					}
				}
			
			} // next j
		} // next i



	FREE_INT(Pts);
	FREE_INT(form_value);
	if (f_v) {
		cout << "blt_generator::compute_adjacency_list_fast done" << endl;
		}
}



void blt_generator::compute_colors(INT orbit_at_level, 
	INT *starter, INT starter_sz, 
	INT special_line, 
	INT *candidates, INT nb_candidates, 
	INT *&point_color, INT &nb_colors, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT p1, p2;
	INT v1[5];
	INT v2[5];
	INT v3[5];
	INT *pts_on_special_line;
	INT idx, i;


	if (f_v) {
		cout << "blt_generator::compute_colors" << endl;
		}
	O->unrank_line(p1, p2, special_line, 0/*verbose_level*/);
	if (f_vv) {
		cout << "after unrank_line " << special_line << ":" << endl;
		cout << "p1=" << p1 << " p2=" << p2 << endl;
		}
	O->unrank_point(v1, 1, p1, 0);
	O->unrank_point(v2, 1, p2, 0);
	if (f_vv) {
		cout << "p1=" << p1 << " ";
		INT_vec_print(cout, v1, 5);
		cout << endl;
		cout << "p2=" << p2 << " ";
		INT_vec_print(cout, v2, 5);
		cout << endl;
		}
	if (p1 != starter[0]) {
		cout << "p1 != starter[0]" << endl;
		exit(1);
		}
	
	pts_on_special_line = NEW_INT(q + 1);
	O->points_on_line(p1, p2, pts_on_special_line, 0/*verbose_level*/);
	
	if (f_vv) {
		cout << "pts_on_special_line:" << endl;
		INT_vec_print(cout, pts_on_special_line, q + 1);
		cout << endl;
		}

	if (!INT_vec_search(pts_on_special_line, q + 1, starter[0], idx)) {
		cout << "cannot find the first point on the line" << endl;
		exit(1);
		}
	for (i = idx; i < q + 1; i++) {
		pts_on_special_line[i] = pts_on_special_line[i + 1];
		}
	if (f_vv) {
		cout << "pts_on_special_line without the first starter point:" << endl;
		INT_vec_print(cout, pts_on_special_line, q);
		cout << endl;
		}
	
	INT a, b, t, c, j, h;
	INT *starter_t;
	
	starter_t = NEW_INT(starter_sz);
	starter_t[0] = -1;
	for (i = 1; i < starter_sz; i++) {
		O->unrank_point(v3, 1, starter[i], 0);
		a = O->evaluate_bilinear_form(v1, v3, 1);
		b = O->evaluate_bilinear_form(v2, v3, 1);
		if (a == 0) {
			cout << "a == 0, this should not be" << endl;
			exit(1);
			}
		// <v3,t*v1+v2> = t*<v3,v1>+<v3,v2> = t*a+b = 0
		// Thus, t = -b/a
		t = O->F->mult(O->F->negate(b), O->F->inverse(a));
		starter_t[i] = t;
		}

	if (f_vv) {
		cout << "starter_t:" << endl;
		INT_vec_print(cout, starter_t, starter_sz);
		cout << endl;
		}

	INT *free_pts;
	INT *open_colors;
	INT *open_colors_inv;

	free_pts = NEW_INT(q);
	open_colors = NEW_INT(q);
	open_colors_inv = NEW_INT(q);

	point_color = NEW_INT(nb_candidates);

	nb_colors = q - starter_sz + 1;
	j = 0;
	for (i = 0; i < q; i++) {
		for (h = 1; h < starter_sz; h++) {
			if (starter_t[h] == i)
				break;
			}
		if (h == starter_sz) {
			free_pts[j] = pts_on_special_line[i];
			open_colors[j] = i;
			j++;
			}
		}
	if (j != nb_colors) {
		cout << "extension_data::setup error: j != nb_colors" << endl;
		exit(1);
		}
	if (f_vv) {
		cout << "The " << nb_colors << " free points are :" << endl;
		INT_vec_print(cout, free_pts, nb_colors);
		cout << endl;
		cout << "The " << nb_colors << " open colors are :" << endl;
		INT_vec_print(cout, open_colors, nb_colors);
		cout << endl;
		}
	for ( ; j < q; j++) {
		open_colors[j] = starter_t[j - nb_colors + 1];
		}
	if (f_vv) {
		cout << "open_colors :" << endl;
		INT_vec_print(cout, open_colors, q);
		cout << endl;
		}
	for (i = 0; i < q; i++) {
		j = open_colors[i];
		open_colors_inv[j] = i;
		}
	if (f_vv) {
		cout << "open_colors_inv :" << endl;
		INT_vec_print(cout, open_colors_inv, q);
		cout << endl;
		}


	for (i = 0; i < nb_candidates; i++) {
		O->unrank_point(v3, 1, candidates[i], 0);
		a = O->evaluate_bilinear_form(v1, v3, 1);
		b = O->evaluate_bilinear_form(v2, v3, 1);
		if (a == 0) {
			cout << "a == 0, this should not be" << endl;
			exit(1);
			}
		// <v3,t*v1+v2> = t*<v3,v1>+<v3,v2> = t*a+b = 0
		// Thus, t = -b/a
		t = O->F->mult(O->F->negate(b), O->F->inverse(a));
		c = open_colors_inv[t];
		if (c >= nb_colors) {
			cout << "c >= nb_colors" << endl;
			cout << "i=" << i << endl;
			cout << "candidates[i]=" << candidates[i] << endl;
			cout << "as vector: ";
			INT_vec_print(cout, v3, 5);
			cout << endl;
			cout << "a=" << a << endl;
			cout << "b=" << b << endl;
			cout << "t=" << t << endl;
			cout << "c=" << c << endl;
			cout << "nb_colors=" << nb_colors << endl;
			
			exit(1);
			}
		point_color[i] = c;
		}

	if (f_vv) {
		cout << "point colors:" << endl;
		INT_vec_print(cout, point_color, nb_candidates);
		cout << endl;
		}

	FREE_INT(pts_on_special_line);
	FREE_INT(starter_t);
	FREE_INT(free_pts);
	FREE_INT(open_colors);
	FREE_INT(open_colors_inv);
	if (f_v) {
		cout << "blt_generator::compute_colors done" << endl;
		}
}



void blt_generator::write_problem_to_file_wassermann(
	const BYTE *output_prefix, 
	INT *Incma, 
	INT nb_free_pts, INT nb_candidates, 
	INT starter_level, INT starter_case, 
	INT width)
{
	BYTE str[1000];
	BYTE fname[1000];
	INT i, j, d, nb;
	
	sprintf(str, "%sBLT_%ld_level_%ld_case_%s0%ldld.txt", output_prefix, q, starter_level, "%", width);
	cout << "str=" << str << endl;
	sprintf(fname, str, starter_case);

	nb = q + 1 - starter_level;
	{
	ofstream fp(fname);
	fp << nb_free_pts << " " << nb_candidates << " " << nb << endl;

	for (i = 0; i < nb_free_pts; i++) {

		fp << "EQ 1 ";
		d = 0;
		for (j = 0; j < nb_candidates; j++) {
			if (Incma[i * nb_candidates + j]) {
				d++;
				}
			}

		fp << d;
		for (j = 0; j < nb_candidates; j++) {
			if (Incma[i * nb_candidates + j]) {
				fp << " " << j;
				}
			}
		fp << endl;
		}
	}
	cout << "written file " << fname << " of size " << file_size(fname) << endl;
}

void blt_generator::early_test_func(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT i, j, a;
	INT f_OK;
	INT v[5];
	INT *v1, *v2, *v3;
	INT m1[5];
	INT m3[5];
	INT *Pts;
	INT *Candidates;
	INT two;
	INT fxy, fxz, fyz;
		
	if (f_v) {
		cout << "blt_generator::early_test_func checking set ";
		print_set(cout, len, S);
		cout << endl;
		cout << "candidate set of size " << nb_candidates << ":" << endl;
		INT_vec_print(cout, candidates, nb_candidates);
		cout << endl;
		if (f_vv) {
			for (i = 0; i < nb_candidates; i++) {
				O->unrank_point(v, 1, candidates[i], 2/*verbose_level - 4*/);
				cout << "candidate " << i << "=" << candidates[i] << ": ";
				INT_vec_print(cout, v, 5);
				cout << endl;
				}
			}
		}
	Pts = NEW_INT(len * 5);
	Candidates = NEW_INT(nb_candidates * 5);
	for (i = 0; i < len; i++) {
		O->unrank_point(Pts + i * 5, 1, S[i], 0/*verbose_level - 4*/);
		}
	for (i = 0; i < nb_candidates; i++) {
		O->unrank_point(Candidates + i * 5, 1, candidates[i], 0/*verbose_level - 4*/);
		}
	
	two = O->F->add(1, 1);


	nb_good_candidates = 0;
	
	for (j = 0; j < nb_candidates; j++) {
		v1 = Pts;
		v3 = Candidates + j * 5;

		m1[0] = O->F->mult(two, v1[0]);
		m1[1] = v1[2];
		m1[2] = v1[1];
		m1[3] = v1[4];
		m1[4] = v1[3];

		//fxz = evaluate_bilinear_form(v1, v3, 1);
		// too slow !!!
		fxz = O->F->add5(
				O->F->mult(m1[0], v3[0]), 
				O->F->mult(m1[1], v3[1]), 
				O->F->mult(m1[2], v3[2]), 
				O->F->mult(m1[3], v3[3]), 
				O->F->mult(m1[4], v3[4]) 
			);

		m3[0] = O->F->mult(two, v3[0]);
		m3[1] = v3[2];
		m3[2] = v3[1];
		m3[3] = v3[4];
		m3[4] = v3[3];

		f_OK = TRUE;
		for (i = 1; i < len; i++) {
			//fxy = evaluate_bilinear_form(v1, v2, 1);

			v2 = Pts + i * 5;
		
			fxy = O->F->add5(
				O->F->mult(m1[0], v2[0]), 
				O->F->mult(m1[1], v2[1]), 
				O->F->mult(m1[2], v2[2]), 
				O->F->mult(m1[3], v2[3]), 
				O->F->mult(m1[4], v2[4]) 
				);
		
			//fyz = evaluate_bilinear_form(v2, v3, 1);
			fyz = O->F->add5(
					O->F->mult(m3[0], v2[0]), 
					O->F->mult(m3[1], v2[1]), 
					O->F->mult(m3[2], v2[2]), 
					O->F->mult(m3[3], v2[3]), 
					O->F->mult(m3[4], v2[4]) 
				);

			a = O->F->product3(fxy, fxz, fyz);

			if (a == 0) {
				f_OK = FALSE;
				break;
				}
			if (O->f_is_minus_square[a]) {
				f_OK = FALSE;
				break;
				}

			}
		if (f_OK) {
			good_candidates[nb_good_candidates++] = candidates[j];
			}
		}
	FREE_INT(Pts);
	FREE_INT(Candidates);
}

INT blt_generator::check_function_incremental(INT len, INT *S, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, a;
	INT f_OK;
	INT *v1, *v2, *v3;
	INT m1[5];
	INT m3[5];
	INT *Pts;
	INT two;
	INT fxy, fxz, fyz;
		
	if (f_v) {
		cout << "blt_generator::check_function_incremental checking set ";
		print_set(cout, len, S);
		cout << endl;
		}

	Pts = NEW_INT(len * 5);
	for (i = 0; i < len; i++) {
		O->unrank_point(Pts + i * 5, 1, S[i], 0/*verbose_level - 4*/);
		}

	two = O->F->add(1, 1);

	v1 = Pts;
	v3 = Pts + (len - 1) * 5;

	m1[0] = O->F->mult(two, v1[0]);
	m1[1] = v1[2];
	m1[2] = v1[1];
	m1[3] = v1[4];
	m1[4] = v1[3];

	//fxz = evaluate_bilinear_form(v1, v3, 1);
	// too slow !!!
	fxz = O->F->add5(
			O->F->mult(m1[0], v3[0]), 
			O->F->mult(m1[1], v3[1]), 
			O->F->mult(m1[2], v3[2]), 
			O->F->mult(m1[3], v3[3]), 
			O->F->mult(m1[4], v3[4]) 
		);

	m3[0] = O->F->mult(two, v3[0]);
	m3[1] = v3[2];
	m3[2] = v3[1];
	m3[3] = v3[4];
	m3[4] = v3[3];

	f_OK = TRUE;
	for (i = 1; i < len - 1; i++) {
		//fxy = evaluate_bilinear_form(v1, v2, 1);

		v2 = Pts + i * 5;
		
		fxy = O->F->add5(
			O->F->mult(m1[0], v2[0]), 
			O->F->mult(m1[1], v2[1]), 
			O->F->mult(m1[2], v2[2]), 
			O->F->mult(m1[3], v2[3]), 
			O->F->mult(m1[4], v2[4]) 
			);
		
		//fyz = evaluate_bilinear_form(v2, v3, 1);
		fyz = O->F->add5(
				O->F->mult(m3[0], v2[0]), 
				O->F->mult(m3[1], v2[1]), 
				O->F->mult(m3[2], v2[2]), 
				O->F->mult(m3[3], v2[3]), 
				O->F->mult(m3[4], v2[4]) 
			);

		a = O->F->product3(fxy, fxz, fyz);

		if (a == 0) {
			f_OK = FALSE;
			break;
			}
		
		if (O->f_is_minus_square[a]) {
			f_OK = FALSE;
			break;
			}

		}
	FREE_INT(Pts);
	return f_OK;
}

void blt_generator::Law_71(INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT set[100];
	INT q = 71;
	//matrix_group *M;
	//orthogonal *O;

	if (f_v) {
		cout << "Law_71" << endl;
		}

	action_on_orthogonal *AO = A->G.AO;
	orthogonal *O;

	//M = A->subaction->G.matrix_grp;
	//O = M->O;
	O = AO->O;
	//M = A->subaction->G.matrix_grp;
	//O = M->O;

	create_Law_71_BLT_set(O, set, verbose_level);
#if 0
	if (!G->check_conditions(cout, q + 1, set, verbose_level)) {
		cout << "the set is not a BLT set" << endl;
		exit(1);
		}
	cout << "BLT test passed" << endl;
#endif

	
	write_set_to_file("Law71.txt", set, q + 1, verbose_level);

#if 0
	r = G->open_database_and_identify_object(set, G->transporter, 
		G->f_use_implicit_fusion, verbose_level);
		
	cout << "Law_71 identified as r=" << r << endl;
#endif
}





void blt_generator::report(isomorph &Iso, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	BYTE fname[1000];

	if (f_v) {
		cout << "blt_generator::report" << endl;
		}
	sprintf(fname, "report_BLT_%ld.tex", q);

	{
	ofstream f(fname);
	INT f_book = TRUE;
	INT f_title = TRUE;
	BYTE title[1000];
	const BYTE *author = "Anton Betten";
	INT f_toc = TRUE;
	INT f_landscape = FALSE;
	INT f_12pt = FALSE;
	INT f_enlarged_page = TRUE;
	INT f_pagenumbers = TRUE;

	sprintf(title, "BLT-sets of Q$(4,%ld)$", q);
	cout << "Writing file " << fname << " with " << Iso.Reps->count << " BLT-sets:" << endl;
	latex_head(f, f_book, f_title, 
		title, author, 
		f_toc, f_landscape, f_12pt, f_enlarged_page, f_pagenumbers);

	f << "\\chapter{Summary}" << endl << endl;
	f << "There are " << Iso.Reps->count << " BLT-sets." << endl << endl;


	//Iso.setup_and_open_solution_database(verbose_level - 1);

	INT i, first, c, id;
	INT u, v, h, rep, tt;
	longinteger_object go;
	INT data[1000];
	INT data2[1000];


	longinteger_object *Ago, *Ago_induced;

	Ago = new longinteger_object[Iso.Reps->count];
	Ago_induced = new longinteger_object[Iso.Reps->count];


	for (h = 0; h < Iso.Reps->count; h++) {
		if (f_v) {
			cout << "blt_generator::report looking at representative h=" << h << endl;
			}
		rep = Iso.Reps->rep[h];
		first = Iso.orbit_fst[rep];
		c = Iso.starter_number[first];
		id = Iso.orbit_perm[first];		
		Iso.load_solution(id, data);

		sims *Stab;
		
		Stab = Iso.Reps->stab[h];

		Iso.Reps->stab[h]->group_order(Ago[h]);
		//f << "Stabilizer has order $";
		//go.print_not_scientific(f);
		if (f_v) {
			cout << "blt_generator::report computing induced action on the set (in data)" << endl;
			}
		Iso.induced_action_on_set(Stab, data, 2 /*verbose_level*/);
		if (f_v) {
			cout << "blt_generator::report induced action on the set (in data) computed" << endl;
			}
		
			
		Iso.AA->group_order(Ago_induced[h]);
		}


	cout << "Computing intersection and plane invariants" << endl;
	INT **intersection_type;
	INT *highest_intersection_number;
	INT **intersection_matrix;
	INT *nb_planes;

	set_of_sets *Sos;
	set_of_sets *Sos2;
	set_of_sets *Sos3;

	decomposition *D2;
	decomposition *D3;

	grassmann *G;
	projective_space *P;
	//INT f_semilinear = TRUE;
	INT set_size = q + 1;

	P = new projective_space;
	
	if (f_v) {
		cout << "before P->init" << endl;
		}

#if 0
	if (is_prime(q)) {
		f_semilinear = FALSE;
		}
#endif


	P->init(4, F, 
		FALSE /* f_init_incidence_structure */, 
		verbose_level);

	if (f_v) {
		cout << "after P->init" << endl;
		}

	G = new grassmann;

	G->init(5, 3, F, 0 /*verbose_level - 2*/);


	longinteger_object **R;
	INT **Sos2_idx;
	INT **Sos3_idx;

	Sos = new set_of_sets[Iso.Reps->count];
	Sos2 = new set_of_sets[Iso.Reps->count];
	Sos3 = new set_of_sets[Iso.Reps->count];
	D2 = new decomposition[Iso.Reps->count];
	D3 = new decomposition[Iso.Reps->count];
	R = new plonginteger_object[Iso.Reps->count];
	Sos2_idx = NEW_PINT(Iso.Reps->count);
	Sos3_idx = NEW_PINT(Iso.Reps->count);

	if (f_v) {
		cout << "blt_generator::report computing invariants" << endl;
		}
	intersection_type = NEW_PINT(Iso.Reps->count);
	highest_intersection_number = NEW_INT(Iso.Reps->count);
	intersection_matrix = NEW_PINT(Iso.Reps->count);
	nb_planes = NEW_INT(Iso.Reps->count);
	for (h = 0; h < Iso.Reps->count; h++) {
		if (f_v) {
			cout << "blt_generator::report looking at representative h=" << h << endl;
			}
		rep = Iso.Reps->rep[h];
		first = Iso.orbit_fst[rep];
		c = Iso.starter_number[first];
		id = Iso.orbit_perm[first];		
		Iso.load_solution(id, data);


		INT v5[5];

		for (i = 0; i < set_size; i++) {
			O->unrank_point(v5, 1, data[i], 0 /* verbose_level */);
			data2[i] = P->rank_point(v5);
			}


		if (f_v) {
			cout << "blt_generator::report before P->plane_intersections" << endl;
			}
		P->plane_intersections(G, 
			data2, set_size, R[h], Sos[h], verbose_level);


		if (f_v) {
			cout << "blt_generator::report before intersection_matrix" << endl;
			}
		Sos[h].intersection_matrix(
			intersection_type[h], highest_intersection_number[h], 
			intersection_matrix[h], nb_planes[h], 
			verbose_level);
		
		if (f_v) {
			cout << "blt_generator::report before extract_largest_sets" << endl;
			}
		Sos[h].extract_largest_sets(Sos2[h], Sos2_idx[h], verbose_level);

		if (f_v) {
			cout << "blt_generator::report before remove_sets_of_given_size" << endl;
			}
		Sos[h].remove_sets_of_given_size(3, Sos3[h], Sos3_idx[h], verbose_level);

		if (f_v) {
			cout << "blt_generator::report before Sos2[h].compute_tdo_decomposition" << endl;
			}
		Sos2[h].compute_tdo_decomposition(D2[h], verbose_level);
		

		D2[h].get_row_scheme(verbose_level);
		D2[h].get_col_scheme(verbose_level);
		if (Sos3[h].nb_sets) {
			if (f_v) {
				cout << "blt_generator::report before Sos3[h].compute_tdo_decomposition" << endl;
				}
			Sos3[h].compute_tdo_decomposition(D3[h], verbose_level);
			D3[h].get_row_scheme(verbose_level);
			D3[h].get_col_scheme(verbose_level);
			}
#if 0
		P->plane_intersection_invariant(G, 
			data2, set_size, 
			intersection_type[h], highest_intersection_number[h], 
			intersection_matrix[h], nb_planes[h], 
			verbose_level);
#endif
		
		}


	cout << "Computing intersection and plane invariants done" << endl;

	f << "\\chapter{Invariants}" << endl << endl;

	f << "\\chapter{The BLT-Sets}" << endl << endl;

	f << "\\clearpage" << endl << endl;


	for (h = 0; h < Iso.Reps->count; h++) {
		rep = Iso.Reps->rep[h];
		first = Iso.orbit_fst[rep];
		c = Iso.starter_number[first];
		id = Iso.orbit_perm[first];		
		Iso.load_solution(id, data);



		f << "\\section{Isomorphism Type " << h << "}" << endl;
		f << "\\bigskip" << endl;

		if (Iso.Reps->stab[h]) {
			Iso.Reps->stab[h]->group_order(go);
			f << "Stabilizer has order $";
			go.print_not_scientific(f);
			f << "$\\\\" << endl;
			}
		else {
			//cout << endl;
			}

		INT a, j;
		f << "Plane intersection type is ";
		for (i = highest_intersection_number[h]; i >= 0; i--) {

			a = intersection_type[h][i];
			if (a == 0) 
				continue;
			f << "$" << i;
			if (a > 9) {
				f << "^{" << a << "}";
				}
			else if (a > 1) {
				f << "^" << a;
				}
#if 0
			if (i < nb_types - 1)
				f << ",\\,";
#endif
			f << "$ ";
			}
		f << "\\\\" << endl;
		f << "Plane invariant is ";

		if (nb_planes[h] < 10) {
			f << "$$";
			f << "\\left[" << endl;
			f << "\\begin{array}{*{" << nb_planes[h] << "}{c}}" << endl;
			for (i = 0; i < nb_planes[h]; i++) {
				for (j = 0; j < nb_planes[h]; j++) {
					f << intersection_matrix[h][i * nb_planes[h] + j];
					if (j < nb_planes[h] - 1) {
						f << " & ";
						}
					}
				f << "\\\\" << endl;
				}
			f << "\\end{array}" << endl;
			f << "\\right]" << endl;
			f << "$$" << endl;
			}
		else {
			f << "too big (" << nb_planes[h] << " planes)\\\\" << endl;
			}

		INT f_enter_math = FALSE;
		INT f_print_subscripts = TRUE;
		
		f << "$$" << endl;
		D2[h].print_row_decomposition_tex(
			f, f_enter_math, f_print_subscripts, verbose_level - 1);
		f << "\\quad" << endl;
		D2[h].print_column_decomposition_tex(
			f, f_enter_math, f_print_subscripts, verbose_level - 1);
		f << "$$" << endl;
		D2[h].Stack->print_classes_tex(f);
		
		if (Sos3[h].nb_sets) {
			f << "$$" << endl;

			D3[h].print_row_decomposition_tex(
				f, f_enter_math, f_print_subscripts, verbose_level - 1);
			f << "$$" << endl;
			f << "$$" << endl;
			D3[h].print_column_decomposition_tex(
				f, f_enter_math, f_print_subscripts, verbose_level - 1);
			f << "$$" << endl;
			D3[h].Stack->print_classes_tex(f);

			INT t, fst_col, fst, len, u, a;
			
			fst_col = D3[h].Stack->startCell[1];
			for (t = 0; t < D3[h].Stack->ht; t++) {
				if (!D3[h].Stack->is_col_class(t)) {
					continue;
					}
				f << "Column cell " << t << ":\\\\" << endl;
				len = D3[h].Stack->cellSize[t];
				fst = D3[h].Stack->startCell[t];
				INT *Cell;
				Cell = NEW_INT(len);
				for (u = 0; u < len; u++) {
					a = D3[h].Stack->pointList[fst + u] - fst_col;
					Cell[u] = a;
					}
				INT_vec_heapsort(Cell, len);
#if 0
				for (u = 0; u < len; u++) {
					a = Cell[u];
					b = Sos3_idx[h][a];
					f << a << " (rank = ";
					R[h][b].print_not_scientific(f);
					f << ") = ";
					G->unrank_longinteger(R[h][b], 0 /* verbose_level */);
					f << "$\\left[" << endl;
					f << "\\begin{array}{*{" << 5 << "}{c}}" << endl;
					for (i = 0; i < 3; i++) {
						for (j = 0; j < 5; j++) {
							c = G->M[i * 5 + j];
							f << c;
							if (j < 4) {
								f << "&";
								}
							}
						f << "\\\\" << endl;
						}
					f << "\\end{array}" << endl;
					f << "\\right]$\\\\" << endl;
					}
#endif
				FREE_INT(Cell);
				}
			}
		


		sims *Stab;
		
		Stab = Iso.Reps->stab[h];

		if (f_v) {
			cout << "blt_generator::report computing induced action on the set (in data)" << endl;
			}
		Iso.induced_action_on_set(Stab, data, 0 /*verbose_level*/);
		
		longinteger_object go1;
			
		Iso.AA->group_order(go1);
		cout << "action " << Iso.AA->label << " computed, group order is " << go1 << endl;

		f << "Order of the group that is induced on the object is ";
		f << "$";
		go1.print_not_scientific(f);
		f << "$\\\\" << endl;
		
		{
		INT nb_ancestors;
		nb_ancestors = Iso.UF->count_ancestors();
		
		f << "Number of ancestors on $" << Iso.level << "$-sets is " << nb_ancestors << ".\\\\" << endl;

		INT *orbit_reps;
		INT nb_orbits;
		strong_generators *Strong_gens;
		//vector_ge SG;
		//INT *tl;
			
		Strong_gens = new strong_generators;
		//tl = NEW_INT(Iso.AA->base_len);
		Strong_gens->init_from_sims(Iso.AA->Sims, 0);
		//Iso.AA->Sims->extract_strong_generators_in_order(SG, tl, verbose_level);
		orbits_on_k_sets(Iso.AA, Iso.AA, Strong_gens /* SG, tl */, 
			Iso.level, orbit_reps, nb_orbits, verbose_level);

		f << "Number of orbits on $" << Iso.level << "$-sets is " << nb_orbits << ".\\\\" << endl;
		FREE_INT(orbit_reps);
		//FREE_INT(tl);
		delete Strong_gens;
		}

		schreier Orb;
		//longinteger_object go2;
		
		Iso.AA->compute_all_point_orbits(Orb, Stab->gens, verbose_level - 2);
		f << "With " << Orb.nb_orbits << " orbits on the object\\\\" << endl;

		classify C_ol;

		C_ol.init(Orb.orbit_len, Orb.nb_orbits, FALSE, 0);

		f << "Orbit lengths: ";
		//INT_vec_print(f, Orb.orbit_len, Orb.nb_orbits);
		C_ol.print_naked_tex(f, FALSE /* f_reverse */);
		f << " \\\\" << endl;
	
		tt = (target_size + 3) / 4;

		f << "The points by ranks:\\\\" << endl;
		f << "\\begin{center}" << endl;

		for (u = 0; u < 4; u++) {
			f << "\\begin{tabular}[t]{|c|c|}" << endl;
			f << "\\hline" << endl;
			f << "$i$ & Rank \\\\" << endl;
			f << "\\hline" << endl;
			for (i = 0; i < tt; i++) {
				v = u * tt + i;
				if (v < target_size) {
					f << "$" << v << "$ & $" << data[v] << "$ \\\\" << endl;
					}
				}
			f << "\\hline" << endl;
			f << "\\end{tabular}" << endl;
			}
		f << "\\end{center}" << endl; 

		f << "The points:\\\\" << endl;
		INT v5[5];
		for (i = 0; i < target_size; i++) {
			O->unrank_point(v5, 1, data[i], 0 /* verbose_level */);
			//Grass->unrank_INT(data[i], 0/*verbose_level - 4*/);
			if ((i % 4) == 0) {
				if (i) {
					f << "$$" << endl;
					}
				f << "$$" << endl;
				}
			//f << "\\left[" << endl;
			//f << "\\begin{array}{c}" << endl;
			f << "P_{" << i /*data[i]*/ << "}=";
			INT_vec_print(f, v5, 5);
#if 0
			for (u = 0; u < 5; u++) {
				for (v = 0; v < n; v++) {
					f << Grass->M[u * n + v];
					}
				f << "\\\\" << endl;
				}
#endif
			//f << "\\end{array}" << endl;
			//f << "\\right]" << endl;
			}
		f << "$$" << endl;


		longinteger_object so;

		Stab->group_order(so);
		f << "Stabilizer of order ";
		so.print_not_scientific(f);
		f << " is generated by:\\\\" << endl;
		for (i = 0; i < Stab->gens.len; i++) {
		
			INT *fp, n;
		
			fp = NEW_INT(A->degree);
			n = A->find_fixed_points(Stab->gens.ith(i), fp, 0);
			//cout << "with " << n << " fixed points" << endl;
			FREE_INT(fp);

			f << "$$ g_{" << i + 1 << "}=" << endl;
			A->element_print_latex(Stab->gens.ith(i), f);
			f << "$$" << endl << "with " << n << " fixed points" << endl;
			}



		//report_stabilizer(Iso, f, h /* orbit */, 0 /* verbose_level */);


		}


	BYTE prefix[1000];
	BYTE label_of_structure_plural[1000];

	sprintf(prefix, "BLT_%ld", q);
	sprintf(label_of_structure_plural, "BLT-Sets");
	isomorph_report_data_in_source_code_inside_tex(Iso, 
		prefix, label_of_structure_plural, f, verbose_level);



	//Iso.close_solution_database(verbose_level - 1);



	latex_foot(f);
	//FREE_INT(Rk_of_span);
	delete G;
	delete P;
	delete [] Ago;
	delete [] Ago_induced;
	}

	cout << "Written file " << fname << " of size " << file_size(fname) << endl;
	if (f_v) {
		cout << "blt_generator::report done" << endl;
		}

}

void blt_generator::subset_orbits(isomorph &Iso, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	BYTE fname[1000];

	if (f_v) {
		cout << "blt_generator::subset_orbits" << endl;
		cout << "A->elt_size_in_INT=" << A->elt_size_in_INT << endl;
		}
	sprintf(fname, "report_BLT_%ld_subset_orbits.tex", q);


	Iso.load_table_of_solutions(verbose_level);
	
	Iso.depth_completed = Iso.level /*- 2*/;

	Iso.gen->recreate_schreier_vectors_up_to_level(Iso.level - 1, TRUE /* f_compact */, verbose_level);

	INT i;
	
	if (f_v) {
		for (i = 0; i <= Iso.level + 1; i++) {
			cout << "gen->first_oracle_node_at_level[" << i << "]=" << Iso.gen->first_oracle_node_at_level[i] << endl;
			}
		cout << "Iso.depth_completed=" << Iso.depth_completed << endl;
		}
	Iso.iso_test_init2(verbose_level);


	{
	ofstream f(fname);
	INT f_book = TRUE;
	INT f_title = TRUE;
	BYTE title[1000];
	const BYTE *author = "Anton Betten";
	INT f_toc = TRUE;
	INT f_landscape = FALSE;
	INT f_12pt = FALSE;
	INT f_enlarged_page = TRUE;
	INT f_pagenumbers = TRUE;

	sprintf(title, "BLT-sets of Q$(4,%ld)$", q);
	cout << "Writing file " << fname << " with " << Iso.Reps->count << " BLT-sets:" << endl;
	latex_head(f, f_book, f_title, 
		title, author, 
		f_toc, f_landscape, f_12pt, f_enlarged_page, f_pagenumbers);

	f << "\\chapter{Summary}" << endl << endl;
	f << "There are " << Iso.Reps->count << " BLT-sets." << endl << endl;


	Iso.setup_and_open_solution_database(verbose_level - 1);

	INT h, rep, first, c, id;
	longinteger_object go;
	INT data[1000];
	//INT data2[1000];

	for (h = 0; h < Iso.Reps->count; h++) {
		rep = Iso.Reps->rep[h];
		first = Iso.orbit_fst[rep];
		c = Iso.starter_number[first];
		id = Iso.orbit_perm[first];		
		Iso.load_solution(id, data);



		f << "\\section{Isomorphism Type " << h << "}" << endl;
		f << "\\bigskip" << endl;

		INT_vec_print(cout, data, Iso.size);
		cout << endl;

		sims *Stab;
		
		Stab = Iso.Reps->stab[h];

		if (f_v) {
			cout << "blt_generator::subset_orbits computing induced action on the set (in data)" << endl;
			}
		Iso.induced_action_on_set(Stab, data, 0 /*verbose_level*/);
		
		cout << "data after induced_action_on_set:" << endl;
		INT_vec_print(cout, data, Iso.size);
		cout << endl;
		
		longinteger_object go1;
			
		Iso.AA->group_order(go1);
		cout << "action " << Iso.AA->label << " computed, group order is " << go1 << endl;

		f << "Order of the group that is induced on the object is ";
		f << "$";
		go1.print_not_scientific(f);
		f << "$\\\\" << endl;

		{
		INT *orbit_reps;
		INT nb_orbits;
		//vector_ge SG;
		//INT *tl;
		strong_generators *Strong_gens;
		
		Strong_gens = new strong_generators;
		Strong_gens->init_from_sims(Iso.AA->Sims, 0);
		//tl = NEW_INT(Iso.AA->base_len);
		//Iso.AA->Sims->extract_strong_generators_in_order(SG, tl, verbose_level);
		orbits_on_k_sets(Iso.AA, Iso.AA, Strong_gens /* SG, tl */, 
			Iso.level, orbit_reps, nb_orbits, verbose_level);

		cout << "Orbit reps: nb_orbits=" << nb_orbits << endl;
		INT_matrix_print(orbit_reps, nb_orbits, Iso.level);

		f << "Number of orbits on $" << Iso.level << "$-sets is " << nb_orbits << ".\\\\" << endl;

		INT *rearranged_set;
		INT *transporter;
		INT u;
		INT case_nb;
		INT f_implicit_fusion = FALSE;
		INT cnt_special_orbits;
		INT f_vv = FALSE;
		INT idx;
		
		rearranged_set = NEW_INT(Iso.size);
		transporter = NEW_INT(A->elt_size_in_INT);

		cnt_special_orbits = 0;
		for (u = 0; u < nb_orbits; u++) {
			cout << "orbit " << u << ":" << endl;
			INT_vec_print(cout, orbit_reps + u * Iso.level, Iso.level);
			cout << endl;



			rearrange_subset(Iso.size, Iso.level, data, orbit_reps + u * Iso.level, rearranged_set, 0/*verbose_level - 3*/);
				// in GALOIS/sorting.C


			//INT_vec_print(cout, rearranged_set, Iso.size);
			//cout << endl;
			INT f_failure_to_find_point, f_found;

			A->element_one(transporter, 0);
			case_nb = Iso.trace_set(rearranged_set, transporter, 
				f_implicit_fusion, f_failure_to_find_point, 0 /*verbose_level - 2*/);


			f_found = Iso.find_extension_easy_new(rearranged_set, case_nb, idx, 0 /* verbose_level */);
#if 0
			f_found = Iso.identify_solution_relaxed(prefix, transporter, 
				f_implicit_fusion, orbit_no0, f_failure_to_find_point, 3 /*verbose_level*/);
#endif

			cout << "case_nb=" << case_nb << endl;
			if (f_failure_to_find_point) {
				cout << "blt_generator::subset_orbits f_failure_to_find_point" << endl;
				exit(1);
				}	
			if (!f_found) {
				if (f_vv) {
					cout << "blt_generator::subset_orbits not found" << endl;
					}
				continue;
				}
			cnt_special_orbits++;
			} // next u

		f << "Number of special orbits on $" << Iso.level << "$-sets is " << cnt_special_orbits << ".\\\\" << endl;

		FREE_INT(rearranged_set);
		FREE_INT(transporter);
		FREE_INT(orbit_reps);
		//FREE_INT(tl);
		delete Strong_gens;
		}

		}

	Iso.close_solution_database(verbose_level - 1);



	latex_foot(f);
	//FREE_INT(Rk_of_span);
	}

	cout << "Written file " << fname << " of size " << file_size(fname) << endl;
	if (f_v) {
		cout << "blt_generator::subset_orbits done" << endl;
		}
}


// ####################################################################################
// global functions:
// ####################################################################################



void print_set(ostream &ost, INT len, INT *S, void *data)
{
	blt_generator *Gen = (blt_generator *) data;
	
	//print_vector(ost, S, len);
	Gen->print(ost, S, len);
}

INT check_conditions(INT len, INT *S, void *data, INT verbose_level)
{
	blt_generator *Gen = (blt_generator *) data;
	return Gen->check_conditions(len, S, verbose_level);
}



void early_test_func_callback(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	void *data, INT verbose_level)
{
	blt_generator *BLT = (blt_generator *) data;
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "early_test_func for set ";
		print_set(cout, len, S);
		cout << endl;
		}
	BLT->early_test_func(S, len, 
		candidates, nb_candidates, 
		good_candidates, nb_good_candidates, 
		verbose_level - 2);
	if (f_v) {
		cout << "early_test_func done" << endl;
		}
}

INT check_function_incremental_callback(INT len, INT *S, void *data, INT verbose_level)
{
	blt_generator *BLT = (blt_generator *) data;
	INT f_OK;
	
	f_OK = BLT->check_function_incremental(len, S, verbose_level);
	return f_OK; 
}



void callback_report(isomorph *Iso, void *data, INT verbose_level)
{
	blt_generator *Gen = (blt_generator *) data;
	
	Gen->report(*Iso, verbose_level);
}

void callback_subset_orbits(isomorph *Iso, void *data, INT verbose_level)
{
	blt_generator *Gen = (blt_generator *) data;
	
	Gen->subset_orbits(*Iso, verbose_level);
}






