// regular_ls_generator.C
// 
// Anton Betten
// 1/1/13

#include "orbiter.h"
#include "discreta.h"

#include "regular_ls.h"

void regular_ls_generator::init_basic(int argc, const char **argv, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);

	if (f_v) {
		cout << "regular_ls_generator::init_basic" << endl;
		}

	gen = new generator;
	
	if (f_vv) {
		cout << "regular_ls_generator::init_basic before read_arguments" << endl;
		}

	read_arguments(argc, argv);

	sprintf(base_fname, "rls_%ld_%ld", m, k);


	m2 = (m * (m - 1)) >> 1;
	v1 = NEW_INT(m);

	row_sum = NEW_INT(m);
	pairs = NEW_INT(m2);
	open_rows = NEW_INT(m);
	open_row_idx = NEW_INT(m);
	open_pairs = NEW_INT(m2);
	open_pair_idx = NEW_INT(m2);

}

void regular_ls_generator::read_arguments(int argc, const char **argv)
{
	INT i;
	INT f_m = FALSE;
	INT f_n = FALSE;
	INT f_k = FALSE;
	INT f_r = FALSE;
	
#if 0
	for (i = 1; i < argc; i++) {
		cout << argv[i] << endl;
		}
#endif
	gen->read_arguments(argc, argv, 0);
	
	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-m") == 0) {
			f_m = TRUE;
			m = atoi(argv[++i]);
			cout << "-m " << m << endl;
			}
		else if (strcmp(argv[i], "-n") == 0) {
			f_n = TRUE;
			n = atoi(argv[++i]);
			cout << "-n " << n << endl;
			}
		else if (strcmp(argv[i], "-k") == 0) {
			f_k = TRUE;
			k = atoi(argv[++i]);
			cout << "-k " << k << endl;
			}
		else if (strcmp(argv[i], "-r") == 0) {
			f_r = TRUE;
			r = atoi(argv[++i]);
			cout << "-r " << r << endl;
			}
		}
	if (!f_m) {
		cout << "regular_ls_generator::read_arguments Please use option -m <m>" << endl;
		exit(1);
		}
	if (!f_n) {
		cout << "regular_ls_generator::read_arguments Please use option -n <n>" << endl;
		exit(1);
		}
	if (!f_k) {
		cout << "regular_ls_generator::read_arguments Please use option -k <k>" << endl;
		exit(1);
		}
	if (!f_r) {
		cout << "regular_ls_generator::read_arguments Please use option -r <r>" << endl;
		exit(1);
		}
}

regular_ls_generator::regular_ls_generator()
{
	null();
}

regular_ls_generator::~regular_ls_generator()
{
	freeself();
}

void regular_ls_generator::null()
{
	gen = NULL;
	A = NULL;
	A2 = NULL;
	initial_pair_covering = NULL;
	row_sum = NULL;
	pairs = NULL;
	open_rows = NULL;
	open_row_idx = NULL;
	open_pairs = NULL;
	open_pair_idx = NULL;
}


void regular_ls_generator::freeself()
{
	if (initial_pair_covering) {
		FREE_INT(initial_pair_covering);
		}
	if (row_sum) {
		FREE_INT(row_sum);
		}
	if (pairs) {
		FREE_INT(pairs);
		}
	if (open_rows) {
		FREE_INT(open_rows);
		}
	if (open_row_idx) {
		FREE_INT(open_row_idx);
		}
	if (open_pairs) {
		FREE_INT(open_pairs);
		}
	if (open_pair_idx) {
		FREE_INT(open_pair_idx);
		}
	null();
}

void regular_ls_generator::init_group(INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "regular_ls_generator::init_group" << endl;
		}

	if (f_v) {
		cout << "regular_ls_generator::init_group creating symmetric group of degree " << m << endl;
		}
	A = new action;
	A->init_symmetric_group(m /* degree */, 0 /* verbose_level - 2*/);
	

	if (f_v) {
		cout << "regular_ls_generator::init_group done" << endl;
		}
}

void regular_ls_generator::init_action_on_k_subsets(INT k, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "regular_ls_generator::init_action_on_k_subsets" << endl;
		}

	//regular_ls_generator::onk = onk;

	if (f_v) {
		cout << "regular_ls_generator::init_action_on_k_subsets creating action on k-subsets for k=" << k << endl;
		}
	A2 = new action;
	A2->induced_action_on_k_subsets(*A, k, verbose_level - 2);

	Aonk = A2->G.on_k_subsets;
	
	if (f_v) {
		cout << "regular_ls_generator::init_action_on_k_subsets before A2->induced_action_override_sims" << endl;
		}

	if (f_v) {
		cout << "regular_ls_generator::init_action_on_k_subsets done" << endl;
		}
}

void regular_ls_generator::init_generator(INT target_depth, 
	INT f_has_initial_pair_covering, INT *initial_pair_covering,
	strong_generators *Strong_gens, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	//INT f_vvv = (verbose_level >= 3);
	INT i;

	if (f_v) {
		cout << "regular_ls_generator::init_generator" << endl;
		}
	if (regular_ls_generator::initial_pair_covering) {
		FREE_INT(regular_ls_generator::initial_pair_covering);
		}
	if (gen->f_max_depth) {
		gen->depth = gen->max_depth;
		}
	else {
		gen->depth = target_depth;
		}
	regular_ls_generator::target_depth = target_depth;
	regular_ls_generator::initial_pair_covering = NEW_INT(m2);
	if (f_has_initial_pair_covering) {
		for (i = 0; i < m2; i++) {
			regular_ls_generator::initial_pair_covering[i] = initial_pair_covering[i];
			}
		}
	else {
		for (i = 0; i < m2; i++) {
			regular_ls_generator::initial_pair_covering[i] = FALSE;
			}
		}
	
	if (f_v) {
		cout << "regular_ls_generator::init_generator depth = " << gen->depth << endl;
		}

	strcpy(gen->fname_base, base_fname);
	
	gen->init(A, A2, Strong_gens, gen->depth, 0/*verbose_level - 3*/);
	
#if 0
	// not needed since we have an early_test_func:
	gen->init_check_func(::check_conditions, 
		(void *)this /* candidate_check_data */);
#endif

	// we have an early test function:

	gen->init_early_test_func(
		rls_generator_early_test_function, 
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
		cout << "regular_ls_generator::init_generator calling init_oracle with " << nb_oracle_nodes << " nodes" << endl;
		}
	
	gen->init_oracle(nb_oracle_nodes, verbose_level - 1);

	if (f_vv) {
		cout << "regular_ls_generator::init_generator after init_root_node" << endl;
		}
	
	//cout << "verbose_level = " << verbose_level << endl;
	//cout << "verbose_level_group_theory = " << verbose_level_group_theory << endl;
	
	gen->root[0].init_root_node(gen, 0/*verbose_level - 2*/);
	if (f_v) {
		cout << "regular_ls_generator::init_generator done" << endl;
		}
}

void regular_ls_generator::compute_starter(INT starter_depth, 
	const BYTE *prefix_starter, 
	INT f_lex, INT f_write_candidate_file, 
	INT f_draw_poset, INT f_embedded, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	BYTE cmd[1000];

	if (f_v) {
		cout << "regular_ls_generator::compute_starter" << endl;
		}

	cout << "regular_ls_generator::compute_starter prefix_starter=" << prefix_starter << endl;
	
	sprintf(cmd, "mkdir %s", prefix_starter);
	system(cmd);

	sprintf(gen->fname_base, "%s%s", prefix_starter, base_fname);


	cout << "gen->fname_base=" << gen->fname_base << endl;
	
	gen->f_W = TRUE;
	gen->classify(0 /* from_level */, starter_depth /* to_level */, 
		f_lex, f_write_candidate_file, 
		verbose_level);


	if (f_draw_poset) {
		if (f_v) {
			cout << "arc_generator::compute_starter before gen->draw_poset" << endl;
			}

		gen->draw_poset(base_fname, starter_depth, 0 /* data1 */, f_embedded, 0 /* gen->verbose_level */);
		
		}
	if (f_v) {
		cout << "regular_ls_generator::compute_starter done" << endl;
		}

}

void regular_ls_generator::early_test_func(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	INT verbose_level)
{
	//verbose_level = 10;
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT i, j, a, b, p;
	INT f_OK;

	if (f_v) {
		cout << "regular_ls_generator::early_test_func checking set ";
		print_set(cout, len, S);
		cout << endl;
		cout << "candidate set of size " << nb_candidates << ":" << endl;
		INT_vec_print(cout, candidates, nb_candidates);
		cout << endl;
		}
	INT_vec_zero(row_sum, m);
	INT_vec_copy(initial_pair_covering, pairs, m2);

	if (f_vv) {
		cout << "pairs initially:" << endl;
		INT_vec_print(cout, pairs, m2);
		cout << endl;
		}
	for (i = 0; i < len; i++) {

		unrank_k_subset(S[i], v1, m, k);
		for (a = 0; a < k; a++) {
			row_sum[v1[a]]++;
			for (b = a + 1; b < k; b++) {
				p = ij2k(v1[a], v1[b], m);
				pairs[p] = TRUE;
				}
			}
		
		}
	if (f_vv) {
		cout << "pairs after adding in the chosen sets:" << endl;
		INT_vec_print(cout, pairs, m2);
		cout << endl;
		}
	

	nb_good_candidates = 0;
	
	for (j = 0; j < nb_candidates; j++) {
		f_OK = TRUE;

		if (f_vv) {
			cout << "Testing candidate " << j << " = " << candidates[j] << endl;
			}

		// do some testing
		unrank_k_subset(candidates[j], v1, m, k);
		if (f_vv) {
			cout << "Testing candidate " << j << " = " << candidates[j] << " = ";
			INT_vec_print(cout, v1, k);
			cout << endl;
			}
		for (a = 0; a < k; a++) {
			if (row_sum[v1[a]] == r) {
				f_OK = FALSE;
				break;
				}
			for (b = a + 1; b < k; b++) {
				p = ij2k(v1[a], v1[b], m);
				if (pairs[p]) {
					f_OK = FALSE;
					break;
					}
				}
			if (!f_OK) {
				break;
				}
			}


		if (f_OK) {
			if (f_vv) {
				cout << "Testing candidate " << j << " = " << candidates[j] << " is good" << endl;
				}
			good_candidates[nb_good_candidates++] = candidates[j];
			}
		}
}

INT regular_ls_generator::check_function_incremental(INT len, INT *S, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, a, b, p;
	INT f_OK;
		
	if (f_v) {
		cout << "regular_ls_generator::check_function_incremental checking set ";
		print_set(cout, len, S);
		cout << endl;
		}

	INT_vec_zero(row_sum, m);
	INT_vec_copy(initial_pair_covering, pairs, m2);

	for (i = 0; i < len - 1; i++) {

		unrank_k_subset(S[i], v1, m, k);
		for (a = 0; a < k; a++) {
			row_sum[v1[a]]++;
			for (b = a + 1; b < k; b++) {
				p = ij2k(v1[a], v1[b], m);
				pairs[p] = TRUE;
				}
			}
		
		}

	f_OK = TRUE;
	unrank_k_subset(S[len - 1], v1, m, k);
	for (a = 0; a < k; a++) {
		if (row_sum[v1[a]] == r) {
			f_OK = FALSE;
			break;
			}
		for (b = a + 1; b < k; b++) {
			p = ij2k(v1[a], v1[b], m);
			if (pairs[p]) {
				f_OK = FALSE;
				break;
				}
			}
		if (!f_OK) {
			break;
			}
		}

	return f_OK;
}

void regular_ls_generator::print(ostream &ost, INT *S, INT len)
{
	INT i;
	
	for (i = 0; i < len; i++) {
		ost << S[i] << " ";
		ost << endl;
		}
}

void regular_ls_generator::lifting_prepare_function_new(exact_cover *E, INT starter_case, 
	INT *candidates, INT nb_candidates, strong_generators *Strong_gens, 
	diophant *&Dio, INT *&col_labels, 
	INT &f_ruled_out, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT i, a, h1, h2, p, idx;
	INT nb_needed;
	INT nb_open_rows, nb_open_pairs;

	if (f_v) {
		cout << "regular_ls_generator::lifting_prepare_function_new nb_candidates=" << nb_candidates << endl;
		}

	nb_needed = n /* target_size */ - E->starter_size;
	f_ruled_out = FALSE;


	INT_vec_zero(row_sum, m);
	INT_vec_copy(initial_pair_covering, pairs, m2);

	if (f_vv) {
		cout << "pairs initially:" << endl;
		INT_vec_print(cout, pairs, m2);
		cout << endl;
		}
	for (i = 0; i < E->starter_size; i++) {

		unrank_k_subset(E->starter[i], v1, m, k);
		for (h1 = 0; h1 < k; h1++) {
			row_sum[v1[h1]]++;
			for (h2 = h1 + 1; h2 < k; h2++) {
				p = ij2k(v1[h1], v1[h2], m);
				pairs[p] = TRUE;
				}
			}
		
		}

	nb_open_rows = 0;
	INT_vec_mone(open_row_idx, m);
	for (i = 0; i < m; i++) {
		if (row_sum[i] < r) {
			open_rows[nb_open_rows] = i;
			open_row_idx[i] = nb_open_rows;
			nb_open_rows++;
			}
		}

	nb_open_pairs = 0;
	INT_vec_mone(open_pair_idx, m2);

	for (i = 0; i < m2; i++) {
		if (pairs[i] == FALSE) {
			open_pairs[nb_open_pairs] = i;
			open_pair_idx[i] = nb_open_pairs;
			nb_open_pairs++;
			}
		}

	
	col_labels = NEW_INT(nb_candidates);


	INT_vec_copy(candidates, col_labels, nb_candidates);

	if (E->f_lex) {
		E->lexorder_test(col_labels, nb_candidates, Strong_gens->gens, 
			verbose_level - 2);
		}

	if (f_vv) {
		cout << "regular_ls_generator::lifting_prepare_function_new after lexorder test" << endl;
		cout << "regular_ls_generator::lifting_prepare_function_new nb_candidates=" << nb_candidates << endl;
		}

	// compute the incidence matrix between
	// open rows and open pairs versus candidate blocks:


	INT nb_rows;
	INT nb_cols;

	nb_rows = nb_open_rows + nb_open_pairs;
	nb_cols = nb_candidates;

	Dio = new diophant;
	Dio->open(nb_rows, nb_cols);
	Dio->sum = nb_needed;

	for (i = 0; i < nb_open_rows; i++) {
		Dio->type[i] = t_EQ;
		Dio->RHS[i] = r - row_sum[open_rows[i]];
		}

	for (i = 0; i < nb_open_pairs; i++) {
		Dio->type[nb_open_rows + i] = t_LE;
		Dio->RHS[nb_open_rows + i] = 1;
		}

	Dio->fill_coefficient_matrix_with(0);


	for (i = 0; i < nb_candidates; i++) {
		a = col_labels[i];


		unrank_k_subset(a, v1, m, k);

		for (h1 = 0; h1 < k; h1++) {

			if (row_sum[v1[h1]] == r) {
				cout << "regular_ls_generator::lifting_prepare_function_new row_sum[v1[h1]] == onr" << endl;
				exit(1);
				}
			idx = open_row_idx[v1[h1]];
			Dio->Aij(idx, i) = 1;
			
			for (h2 = h1 + 1; h2 < k; h2++) {
				p = ij2k(v1[h1], v1[h2], m);
				if (pairs[p]) {
					cout << "regular_ls_generator::lifting_prepare_function_new pairs[p]" << endl;
					exit(1);
					}
				idx = open_pair_idx[p];
				Dio->Aij(nb_open_rows + idx, i) = 1;
				}
			}

		}


	
	if (f_v) {
		cout << "regular_ls_generator::lifting_prepare_function_new done" << endl;
		}
}


#if 0
void regular_ls_generator::extend(const BYTE *fname, 
	INT f_single_case, INT single_case, 
	INT N, INT K, INT R, INT f_lambda_reached, INT depth, 
	INT f_lexorder_test, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	INT orbit_at_level;
	INT nb_orbits;

	if (f_v) {
		cout << "regular_ls_generator::extend" << endl;
		cout << "regular_ls_generator::extend N=" << N << " K=" << K << " R=" << R << " f_lambda_reached=" << f_lambda_reached << endl;
		}

	nb_orbits = count_number_of_orbits_in_file(fname, verbose_level - 2);

	if (f_v) {
		cout << "regular_ls_generator::extend nb_orbits = " << nb_orbits << endl;
		}



	if (f_single_case) {
			extend_a_single_case(fname, 
				N, K, R, f_lambda_reached, 
				f_lexorder_test, 
				single_case, nb_orbits, depth, 
				verbose_level);
		}
	else {

		for (orbit_at_level = 0; orbit_at_level < nb_orbits; orbit_at_level++) {
			extend_a_single_case(fname, 
				N, K, R, f_lambda_reached, 
				f_lexorder_test, 
				orbit_at_level, nb_orbits, depth, 
				verbose_level);

			}
		}

	
}

void regular_ls_generator::extend_a_single_case(const BYTE *fname, 
	INT N, INT K, INT R, INT f_lambda_reached, 
	INT f_lexorder_test, 
	INT orbit_at_level, INT nb_orbits, INT depth, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT *starter;
	INT starter_sz;
	sims *Stab;
	strong_generators *Strong_gens;

	
	INT *Mtx;
	//INT *open_pairs;
	//INT *open_pair_idx;
	//INT nb_open_pairs;
	INT i, a, b, p;
	longinteger_object ago;

	if (f_v) {
		cout << "regular_ls_generator::extend_a_single_case reading representative of orbit " << orbit_at_level << " / " << nb_orbits << endl;
		}
		
	INT *set;
	INT *candidates;
	INT nb_candidates;
	INT nb_cases;
	INT nb, u, f_OK;
	

	nb = INT_n_choose_k(m, K);


	Mtx = NEW_INT(m * n);
	//open_pairs = NEW_INT(m2);
	//open_pair_idx = NEW_INT(m2);
	set = NEW_INT(K);
	candidates = NEW_INT(nb);


	A->read_set_and_stabilizer(fname, 
		orbit_at_level, starter, starter_sz, Stab, 
		Strong_gens, 
		nb_cases, 
		verbose_level);
		// ACTION/action.C

	Stab->group_order(ago);
	if (f_v) {
		cout << "regular_ls_generator::extend_a_single_case read representative of orbit " << orbit_at_level << " / " << nb_orbits << " ago=" << ago << endl;
		INT_vec_print(cout, starter, starter_sz);
		cout << endl;
		}

	INT_vec_zero(row_sum, m);
	INT_vec_zero(pairs, m2);
	INT_vec_zero(Mtx, m * n);

	for (i = 0; i < n; i++) {

		unrank_k_subset(starter[i], v1, m, k);
		for (a = 0; a < k; a++) {
			row_sum[v1[a]]++;
			Mtx[v1[a] * n + i] = 1;
			for (b = a + 1; b < k; b++) {
				p = ij2k(v1[a], v1[b], m);
				pairs[p] = TRUE;
				}
			}
		
		}
	if (f_vv) {
		cout << "Mtx:" << endl;
		INT_matrix_print(Mtx, m, n);
		cout << "ago=" << ago << endl;
		cout << "row_sum:";
		INT_vec_print(cout, row_sum, m);
		cout << endl;
		cout << "pairs:";
		INT_vec_print(cout, pairs, m2);
		cout << endl;
		}

	nb_candidates = 0;
	for (u = 0; u < nb; u++) {
		unrank_k_subset(u, set, m, K);
		f_OK = TRUE;
		for (a = 0; a < K; a++) {
			for (b = a + 1; b < K; b++) {
				p = ij2k(set[a], set[b], m);
				if (pairs[p]) {
					f_OK = FALSE;
					break;
					}
				}
			if (!f_OK) {
				break;
				}
			}
		if (f_OK) {
			candidates[nb_candidates++] = u;
			}
		}
	if (f_vv) {
		cout << "There are " << nb_candidates << " candidates" << endl;
		for (i = 0; i < nb_candidates; i++) {
			cout << i << " : " << candidates[i] << " : ";
			unrank_k_subset(candidates[i], set, m, K);
			INT_vec_print(cout, set, K);
			cout << endl;
			}
		}


	init_generator(N, TRUE, pairs, 
		Strong_gens, 
		verbose_level - 1);

	INT t0 = os_ticks();
	INT f_use_invariant_subset_if_available = TRUE;
	INT f_debug = FALSE;
	INT f_implicit_fusion = FALSE;
	
	gen->depth = depth;

	gen->f_W = TRUE;
	//gen->f_prefix = TRUE;
	sprintf(gen->fname_base, "case_%ld", orbit_at_level);

	gen->init_root_node_invariant_subset(
		candidates, nb_candidates, verbose_level - 1);

	gen->main(t0, 500 /*schreier_depth*/, 
		f_use_invariant_subset_if_available, 
		f_implicit_fusion, 
		f_debug, 
		verbose_level + 1);

	INT nb_starters;
	INT orbit_at_depth;
	INT nb_sol, nbs;
	INT *Nb_sol;
	INT **Solutions;

	nb_starters = gen->nb_orbits_at_level(depth);

	cout << "regular_ls_generator::extend_a_single_case after generator_main, there are " << nb_starters << " starters" << endl;

	Nb_sol = NEW_INT(nb_starters);
	Solutions = NEW_PINT(nb_starters);
	nb_sol = 0;
	for (orbit_at_depth = 0; orbit_at_depth < nb_starters; orbit_at_depth++) {

		if (f_v) {
			cout << "regular_ls_generator::extend_a_single_case " << orbit_at_level << " / " << nb_orbits << ", before handle_starter " << orbit_at_depth << " / " << nb_starters << endl;
			}
		handle_starter(fname, 
			N, K, R, f_lambda_reached, 
			f_lexorder_test, 
			orbit_at_level, nb_orbits, 
			orbit_at_depth, nb_starters, depth, 
			pairs,
			Solutions[orbit_at_depth], nbs,  
			verbose_level - 2);
		Nb_sol[orbit_at_depth] = nbs;
		nb_sol += Nb_sol[orbit_at_depth];

		if (f_v) {
			cout << "regular_ls_generator::extend_a_single_case " << orbit_at_level << " / " << nb_orbits << ", after handle_starter " << orbit_at_depth << " / " << nb_starters << " #sol=" << Nb_sol[orbit_at_depth] << " total = " << nb_sol << endl;
			}
		}
	if (f_v) {
		cout << "The " << nb_sol << " solutions are:" << endl;
		for (orbit_at_depth = 0; orbit_at_depth < nb_starters; orbit_at_depth++) {
			if (Nb_sol[orbit_at_depth]) {
				cout << "case " << orbit_at_depth << " has " << Nb_sol[orbit_at_depth] << " solutions:" << endl;
				INT_matrix_print(Solutions[orbit_at_depth], Nb_sol[orbit_at_depth], N);
				}
			}
		
		}

	if (nb_sol) {
		INT print_mod = 1000;

		isomorph_build_db(A, A2, gen, 
			N /* size */, gen->fname_base, (BYTE *) "ISO/", depth /* level */, verbose_level - 2);

		isomorph_init_solutions_from_memory(A, A2, gen, 
			N /* size */, gen->fname_base, (BYTE *) "ISO/", depth /* level */, 
			Solutions, Nb_sol, verbose_level);
	
		isomorph_compute_orbits(A, A2, gen, 
			N /* size */, gen->fname_base, (BYTE *) "ISO/", depth /* level */, 
			verbose_level);
	
		isomorph_testing(A, A2, gen, 
			N /* size */, gen->fname_base, (BYTE *) "ISO/", depth /* level */, 
			FALSE /*f_play_back*/, NULL /*old_event_file*/, print_mod, verbose_level);
		}



#if 0	
	exit(1);
		
	nb_open_pairs = 0;
	for (p = 0; p < m2; p++) {
		if (pairs[p]) {
			open_pair_idx[p] = -1;
			}
		else {
			open_pair_idx[p] = nb_open_pairs;
			open_pairs[nb_open_pairs++] = p;
			}
		}
	if (f_vv) {
		cout << "nb_open_pairs=" << nb_open_pairs << endl;
		cout << "open pairs" << endl;
		INT_vec_print(cout, open_pairs, nb_open_pairs);
		cout << endl;
		}


	INT nb_candidates;
	INT *candidates;

	nb_candidates = 0;
	for (u = 0; u < nb; u++) {
		unrank_k_subset(u, set, m, K);
		f_OK = TRUE;
		for (a = 0; a < K; a++) {
			for (b = a + 1; b < K; b++) {
				p = ij2k(set[a], set[b], m);
				if (pairs[p]) {
					f_OK = FALSE;
					break;
					}
				}
			if (!f_OK) {
				break;
				}
			}
		if (f_OK) {
			candidates[nb_candidates++] = u;
			}
		}
	if (f_vv) {
		cout << "There are " << nb_candidates << " candidates" << endl;
		for (i = 0; i < nb_candidates; i++) {
			cout << i << " : " << candidates[i] << " : ";
			unrank_k_subset(candidates[i], set, m, K);
			INT_vec_print(cout, set, K);
			cout << endl;
			}
		}


		
	INT *C;
	INT *System;
	INT *RHS;
	INT nb_rows, p_idx;

	C = NEW_INT(nb_candidates * K);
	for (u = 0; u < nb_candidates; u++) {
		unrank_k_subset(candidates[u], set, m, K);
		for (i = 0; i < K; i++) {
			C[u * K + i] = set[i];
			}
		}
	if (f_vv) {
		cout << "candidates:" << endl;
		INT_matrix_print(C, nb_candidates, K);
		}
	nb_rows = m + nb_open_pairs;
	System = NEW_INT(nb_rows * nb_candidates);
	RHS = NEW_INT(nb_rows);
	for (i = 0; i < nb_rows; i++) {
		for (j = 0; j < nb_candidates; j++) {
			System[i * nb_candidates + j] = 0;
			}
		RHS[i] = 0;
		}

	for (u = 0; u < nb_candidates; u++) {
		for (i = 0; i < K; i++) {
			a = C[u * K + i];
			System[a * nb_candidates + u] = 1;
			}
		for (a = 0; a < K; a++) {
			for (b = a + 1; b < K; b++) {
				p = ij2k(C[u * K + a], C[u * K + b], m);
				p_idx = open_pair_idx[p];
				if (p_idx < 0) {
					cout << "p_idx < 0" << endl;
					exit(1);
					}
				System[(m + p_idx) * nb_candidates + u] = 1;
				}
			}
		}

	for (i = 0; i < m; i++) {
		RHS[i] = R;
		}
	for (i = 0; i < nb_open_pairs; i++) {
		RHS[m + i] = 1;
		}
	diophant *D;

	D = new diophant;
	D->open(nb_rows, nb_candidates);
	for (i = 0; i < nb_rows; i++) {
		for (j = 0; j < nb_candidates; j++) {
			D->Aij(i, j) = System[i * nb_candidates + j];
			}
		D->RHSi(i) = RHS[i];
		D->f_le[i] = FALSE;
		}
	for (j = 0; j < nb_candidates; j++) {
		D->x_max[j] = 1;
		}
	D->f_x_max = TRUE;
	D->sum = N;

	if (!f_lambda_reached) {
		for (i = 0; i < nb_open_pairs; i++) {
			D->f_le[m + i] = TRUE;
			}
		}

	if (f_vv) {
		cout << "System created:" << endl;
		D->print();
		cout << "before D->solve_all_mckay" << endl;
		}
	D->solve_all_mckay(0 /*verbose_level*/);
	if (f_vv) {
		cout << "after D->solve_all_mckay number of solutions = " << D->_resultanz << endl;
		}


	delete D;


	FREE_INT(RHS);
	FREE_INT(System);
	FREE_INT(C);
#endif



	FREE_INT(starter);
	delete Stab;
	delete Strong_gens;




	FREE_INT(Mtx);
	//FREE_INT(open_pairs);
	//FREE_INT(open_pair_idx);
	FREE_INT(set);
	FREE_INT(candidates);
	FREE_INT(Nb_sol);

	if (f_v) {
		cout << "regular_ls_generator::extend_a_single_case orbit " << orbit_at_level << " / " << nb_orbits << " ago=" << ago << " done, total number of solutions = " << nb_sol << endl;
		}

}

void regular_ls_generator::handle_starter(const BYTE *fname, 
	INT N, INT K, INT R, INT f_lambda_reached, 
	INT f_lexorder_test, 
	INT orbit_at_level, INT nb_orbits, 
	INT orbit_at_depth, INT nb_starters, INT depth, 
	INT *pairs, 
	INT *&Solutions, INT &nb_sol, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT f_vvv = (verbose_level >= 3);
	INT *Mtx;
	INT nb_open_pairs;
	INT *open_pairs;
	INT *open_pair_idx;
	INT *pairs2;
	INT *starter;
	INT *set;
	INT i, j, u, a, b, p, nb, f_OK;
	INT nb_candidates;
	INT *candidates;
	INT max_starter;


	if (f_v) {
		cout << "regular_ls_generator::handle_starter " << orbit_at_level << " / " << nb_orbits << ", " << orbit_at_depth << " / " << nb_starters << endl;
		cout << "verbose_level = " << verbose_level << endl;
		}
	nb = INT_n_choose_k(m, K);
	
	Mtx = NEW_INT(m * depth);
	starter = NEW_INT(m);
	set = NEW_INT(m);
	pairs2 = NEW_INT(m2);
	open_pairs = NEW_INT(m2);
	open_pair_idx = NEW_INT(m2);
	candidates = NEW_INT(nb);

	INT_vec_copy(pairs, pairs2, m2);
	INT_vec_zero(row_sum, m);
	INT_vec_zero(Mtx, m * depth);

	gen->get_set_by_level(depth, orbit_at_depth, starter);
	max_starter = starter[depth - 1];

	if (f_v) {
		cout << "regular_ls_generator::handle_starter " << orbit_at_level << " / " << nb_orbits << ", " << orbit_at_depth << " / " << nb_starters << " starter=";
		INT_vec_print(cout, starter, depth);
		cout << endl;
		if (f_vv) {
			for (i = 0; i < depth; i++) {
				unrank_k_subset(starter[i], v1, m, onk);
				cout << i << " : " << starter[i] << " : ";
				INT_vec_print(cout, v1, onk);
				cout << endl;
				}
			}
		}

	for (i = 0; i < depth; i++) {

		unrank_k_subset(starter[i], v1, m, onk);
		for (a = 0; a < onk; a++) {
			row_sum[v1[a]]++;
			Mtx[v1[a] * depth + i] = 1;
			for (b = a + 1; b < onk; b++) {
				p = ij2k(v1[a], v1[b], m);
				if (pairs2[p]) {
					cout << "regular_ls_generator::handle_starter pairs2[p]" << endl;
					exit(1);
					}
				pairs2[p] = TRUE;
				}
			}
		
		}

	if (f_v) {
		cout << "Mtx:" << endl;
		INT_matrix_print(Mtx, m, depth);
		cout << "row_sum=";
		INT_vec_print(cout, row_sum, m);
		cout << endl;
		cout << "pairs2=";
		INT_vec_print(cout, pairs2, m2);
		cout << endl;
		}

	nb_candidates = 0;
	for (u = 0; u < nb; u++) {
		unrank_k_subset(u, set, m, onk);
		f_OK = TRUE;
		for (a = 0; a < onk; a++) {
			if (row_sum[set[a]] == onr) {
				f_OK = FALSE;
				break;
				}
			for (b = a + 1; b < onk; b++) {
				p = ij2k(set[a], set[b], m);
				if (pairs2[p]) {
					f_OK = FALSE;
					break;
					}
				}
			if (!f_OK) {
				break;
				}
			}
		if (f_OK) {
			candidates[nb_candidates++] = u;
			}
		}
	if (f_vv) {
		cout << "There are " << nb_candidates << " candidates" << endl;
		}
	if (f_vvv) {
		for (i = 0; i < nb_candidates; i++) {
			cout << i << " : " << candidates[i] << " : ";
			unrank_k_subset(candidates[i], set, m, onk);
			INT_vec_print(cout, set, onk);
			cout << endl;
			}
		}
	if (f_lexorder_test) {
		INT nb_candidates2;
		strong_generators *Strong_gens;
		longinteger_object go;
	
		if (f_v) {
			cout << "Doing a lexorder_test" << endl;
			}

		gen->get_stabilizer_generators(Strong_gens,  
			depth, orbit_at_depth, 0 /* verbose_level */);
		Strong_gens->group_order(go);
		if (f_v) {
			cout << "Stabilizer order = " << go << endl;
			cout << "Stabilizer generators:" << endl;
			Strong_gens->print_generators();
			}
		if (f_vv) {
			cout << "Before lexorder_test" << endl;
			}
		A2->lexorder_test(candidates, nb_candidates, nb_candidates2, 
			Strong_gens->gens, max_starter, verbose_level - 3);
		if (f_v) {
			cout << "After lexorder_test nb_candidates=" << nb_candidates2 << " eliminated " << nb_candidates - nb_candidates2 << " candidates" << endl;
			}
		nb_candidates = nb_candidates2;
		delete Strong_gens;
		}


	nb_open_pairs = 0;
	for (p = 0; p < m2; p++) {
		if (pairs2[p]) {
			open_pair_idx[p] = -1;
			}
		else {
			open_pair_idx[p] = nb_open_pairs;
			open_pairs[nb_open_pairs++] = p;
			}
		}
	if (f_vvv) {
		cout << "nb_open_pairs=" << nb_open_pairs << endl;
		cout << "open pairs" << endl;
		INT_vec_print(cout, open_pairs, nb_open_pairs);
		cout << endl;
		cout << "open_pair_idx" << endl;
		INT_vec_print(cout, open_pair_idx, m2);
		cout << endl;
		}


		
	INT *C;
	INT *System;
	INT *RHS;
	INT nb_rows, p_idx;

	C = NEW_INT(nb_candidates * onk);
	for (u = 0; u < nb_candidates; u++) {
		unrank_k_subset(candidates[u], set, m, onk);
		for (i = 0; i < onk; i++) {
			C[u * onk + i] = set[i];
			}
		}
#if 0
	if (f_vv) {
		cout << "candidates:" << endl;
		INT_matrix_print(C, nb_candidates, onk);
		}
#endif
	nb_rows = m + nb_open_pairs;
	System = NEW_INT(nb_rows * nb_candidates);
	RHS = NEW_INT(nb_rows);
	INT_vec_zero(System, nb_rows * nb_candidates);
	INT_vec_zero(RHS, nb_rows);

	for (u = 0; u < nb_candidates; u++) {
		for (i = 0; i < onk; i++) {
			a = C[u * onk + i];
			System[a * nb_candidates + u] = 1;
			}
		for (a = 0; a < onk; a++) {
			for (b = a + 1; b < onk; b++) {
				p = ij2k(C[u * onk + a], C[u * onk + b], m);
				p_idx = open_pair_idx[p];
				if (p_idx < 0) {
					cout << "p_idx < 0" << endl;
					exit(1);
					}
				System[(m + p_idx) * nb_candidates + u] = 1;
				}
			}
		}

	for (i = 0; i < m; i++) {
		RHS[i] = onr - row_sum[i];
		}
	for (i = 0; i < nb_open_pairs; i++) {
		RHS[m + i] = 1;
		}
	diophant *D;

	D = new diophant;
	D->open(nb_rows, nb_candidates);
	for (i = 0; i < nb_rows; i++) {
		for (j = 0; j < nb_candidates; j++) {
			D->Aij(i, j) = System[i * nb_candidates + j];
			}
		D->RHSi(i) = RHS[i];
		D->type[i] = t_EQ;
		}
	for (j = 0; j < nb_candidates; j++) {
		D->x_max[j] = 1;
		}
	D->f_x_max = TRUE;
	D->sum = N - depth;

	if (!f_lambda_reached) {
		for (i = 0; i < nb_open_pairs; i++) {
			D->type[m + i] = t_LE;
			}
		}

	if (f_vvv) {
		cout << "System created:" << endl;
		//D->print();
		D->print_dense();
		cout << "before D->solve_all_mckay" << endl;
		}
	if (nb_candidates) {
		INT nb_backtrack;
		diophant_solve_all_mckay(D, nb_backtrack, 0 /*verbose_level*/);
		nb_sol = D->_resultanz;
		}
	else {
		nb_sol = 0;
		}
	if (f_vv) {
		cout << "after D->solve_all_mckay number of solutions = " << nb_sol << endl;
		}
	if (nb_sol) {
		INT *Sol;
		INT nb_sol1;
		
		Solutions = NEW_INT(nb_sol * N);
		D->get_solutions(Sol, nb_sol1, 0 /* verbose_level */);
		if (depth + D->sum != N) {
			cout << "depth + S->sum != N" << endl;
			exit(1);
			}
		for (i = 0; i < nb_sol; i++) {
			if (FALSE && i == 0) {
				cout << "Solution 0, starter:" << endl;
				INT_vec_print(cout, starter, depth);
				cout << endl;
				}
			for (j = 0; j < depth; j++) {
				Solutions[i * N + j] = starter[j];
				}
			if (FALSE && i == 0) {
				cout << "Solution 0:" << endl;
				INT_vec_print(cout, Solutions + i * N, N);
				cout << endl;
				}
			for (j = 0; j < D->sum; j++) {
				u = Sol[i * D->sum + j];
				Solutions[i * N + depth + j] = candidates[u];
				}
			if (FALSE && i == 0) {
				cout << "Solution 0:" << endl;
				INT_vec_print(cout, Solutions + i * N, N);
				cout << endl;
				}
			}
		FREE_INT(Sol);
		}
	else {
		Solutions = 0;
		}


	delete D;


	FREE_INT(RHS);
	FREE_INT(System);
	FREE_INT(C);




	FREE_INT(candidates);


	FREE_INT(Mtx);
	FREE_INT(starter);
	FREE_INT(set);
	FREE_INT(pairs2);
	FREE_INT(open_pairs);
	FREE_INT(open_pair_idx);
	if (f_v) {
		cout << "regular_ls_generator::handle_starter " << orbit_at_level << " / " << nb_orbits << ", " << orbit_at_depth << " / " << nb_starters << " done, with " << nb_sol << " solutions" << endl;
		}
}
#endif



// ####################################################################################
// global functions:
// ####################################################################################



void print_set(ostream &ost, INT len, INT *S, void *data)
{
	regular_ls_generator *Gen = (regular_ls_generator *) data;
	
	//print_vector(ost, S, len);
	Gen->print(ost, S, len);
}

void rls_generator_early_test_function(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	void *data, INT verbose_level)
{
	regular_ls_generator *Gen = (regular_ls_generator *) data;
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "rls_generator_early_test_function for set ";
		print_set(cout, len, S);
		cout << endl;
		}
	Gen->early_test_func(S, len, 
		candidates, nb_candidates, 
		good_candidates, nb_good_candidates, 
		verbose_level - 2);
	if (f_v) {
		cout << "rls_generator_early_test_function done" << endl;
		}
}

INT check_function_incremental_callback(INT len, INT *S, void *data, INT verbose_level)
{
	regular_ls_generator *Gen = (regular_ls_generator *) data;
	INT f_OK;
	
	f_OK = Gen->check_function_incremental(len, S, verbose_level);
	return f_OK; 
}


void rls_generator_lifting_prepare_function_new(exact_cover *EC, INT starter_case, 
	INT *candidates, INT nb_candidates, strong_generators *Strong_gens, 
	diophant *&Dio, INT *&col_labels, 
	INT &f_ruled_out, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	regular_ls_generator *Gen = (regular_ls_generator *) EC->user_data;

	if (f_v) {
		cout << "rls_generator_lifting_prepare_function_new nb_candidates=" << nb_candidates << endl;
		}

	Gen->lifting_prepare_function_new(EC, starter_case, 
		candidates, nb_candidates, Strong_gens, 
		Dio, col_labels, f_ruled_out, 
		verbose_level - 1);


	if (f_v) {
		cout << "rls_generator_lifting_prepare_function_new nb_rows=" << Dio->m << " nb_cols=" << Dio->n << endl;
		}

	if (f_v) {
		cout << "rls_generator_lifting_prepare_function_new done" << endl;
		}
}



