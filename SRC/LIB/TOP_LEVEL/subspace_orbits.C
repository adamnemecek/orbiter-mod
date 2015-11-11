// subspace_orbits.C
// 
// Anton Betten
//
// started:    January 25, 2010
// moved here: March 29, 2012
// 
//
//

#include "orbiter.h"


void lift_generators_to_subfield_structure(action *Aq, INT n, INT s, finite_field *Fq, strong_generators *&SG, INT verbose_level);

subspace_orbits::subspace_orbits()
{
	P = NULL;
	A_linear = NULL;
	//Mtx = NULL;
	tmp_M = NULL;
	base_cols = NULL;
	Gen = NULL;
	f_wedge_action = FALSE;
	A_wedge = NULL;
	f_PGL2_OnConic = FALSE;
	A_OnConic = NULL;
	f_subfield_structure_action = FALSE;
	f_has_strong_generators = FALSE;
	Strong_gens = NULL;
	f_print_generators = FALSE;
	f_has_extra_test_func = FALSE;
}

subspace_orbits::~subspace_orbits()
{
	if (P) {
		delete P;
		}
	if (A_linear) {
		delete A_linear;
		}
	if (A_OnConic) {
		delete A_OnConic;
		}
	if (A_wedge) {
		delete A_wedge;
		}
	if (tmp_M) {
		FREE_INT(tmp_M);
		}
	if (base_cols) {
		FREE_INT(base_cols);
		}
	if (Gen) {
		delete Gen;
		}
	if (Strong_gens) {
		delete Strong_gens;
		}
}

void subspace_orbits::init_group(INT *group_generator_data, INT group_generator_size, 
	INT f_group_order_target, const BYTE *group_order_target, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	vector_ge gens;

	if (f_v) {
		cout << "subspace_orbits::init_group" << endl;
		}
	A_linear->init_group_from_generators(group_generator_data, group_generator_size, 
		f_group_order_target, group_order_target, 
		&gens, Strong_gens, 
		verbose_level);
	f_has_strong_generators = TRUE;
	if (f_v) {
		cout << "subspace_orbits::init_group done" << endl;
		}
}

void subspace_orbits::init_PGL2q_OnConic(BYTE *prefix, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "subspace_orbits::init_PGL2q_OnConic initializing action of PGL(2,q) on conic" << endl;
		}
	if (!A_linear->f_has_sims) {
		cout << "subspace_orbits::init_PGL2q_OnConic A_linear does not have sims, so we create it" << endl;
		A_linear->create_sims(verbose_level);
		}
	if (!A_linear->f_has_strong_generators) {
		cout << "subspace_orbits::init_PGL2q_OnConic A_linear does not have strong generators" << endl;
		//A_linear->create_sims(verbose_level);
		exit(1);
		}
	A_OnConic = new action;
	A_OnConic->induced_action_by_representation_on_conic(A_linear, 
		FALSE /* f_induce_action */, NULL, 
		verbose_level);
	vector_space_dimension = A_OnConic->G.Rep->dimension;
	if (f_v) {
		cout << "subspace_orbits::init_PGL2q_OnConic vector_space_dimension=" << vector_space_dimension << endl;
		}
	if (f_v) {
		cout << "subspace_orbits::init created action of PGL2_on conic:" << endl;
		A_OnConic->print_info();
		}
	sprintf(prefix, "PGL2_OnConic_%ld_%ld", n, q);
}

void subspace_orbits::init_wedge_action(BYTE *prefix, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "subspace_orbits::init_wedge_action initializing wedge action" << endl;
		}
	if (!A_linear->f_has_sims) {
		cout << "subspace_orbits::init_wedge_action A_linear does not have sims, so we create it" << endl;
		A_linear->create_sims(verbose_level);
		}
	if (!A_linear->f_has_strong_generators) {
		cout << "subspace_orbits::init_wedge_action A_linear does not have strong generators" << endl;
		//>create_sims(verbose_level);
		exit(1);
		}
	A_wedge = new action;
	action_on_wedge_product *AW;


	if (f_v) {
		cout << "subspace_orbits::init_wedge_action before induced_wedge_action:" << endl;
		}
	AW = new action_on_wedge_product;

	AW->init(*A_linear, verbose_level);
	
	vector_space_dimension = AW->wedge_dimension;
	if (f_v) {
		cout << "subspace_orbits::init_wedge_action vector_space_dimension=" << vector_space_dimension << endl;
		}
		
	A_wedge->induced_action_on_wedge_product(A_linear, 
		AW, 
		FALSE /* f_induce_action */, NULL, 
		verbose_level);
	if (f_v) {
		cout << "subspace_orbits::init_wedge_action created wedge action:" << endl;
		A_wedge->print_info();
		}
	sprintf(prefix, "Wedge_%ld_%ld", n, q);
}

void subspace_orbits::init_monomial_group(BYTE *prefix, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT *Elt1;
	sims *S;
	action *A;
	finite_field *F;
	longinteger_domain D;
	longinteger_object target_go;
	INT *go_factored;
	INT n, q;
	vector_ge *gens;
	INT *data;
	INT i, h, hh;
	
	if (f_v) {
		cout << "subspace_orbits::init_monomial_group initializing monomial group" << endl;
		}
	A = A_linear;
	F = P->F;
	q = F->q;
	n = vector_space_dimension;
	if (f_v) {
		cout << "n=" << n << " q=" << q << endl;
		}
	if (A->degree < 1000) {
		cout << "The degree of the action is small, so we list the points:" << endl;
		display_all_PG_elements(n - 1, *F);
		}
	Elt1 = NEW_INT(A->elt_size_in_INT);
	go_factored = NEW_INT(2 * n + 1);
	data = NEW_INT(n * n + 1);

	// group order 
	// = n! * q^n * e if not projective
	// = n! * q^(n-1) * e if projective
	// where e is the degree of the field if f_semilinear is TRUE
	// and e = 1 otherwise
	
	for (i = 0; i < n; i++) {
		go_factored[i] = n - i;
		}
	for (i = 0; i < n; i++) {
		if (i == n - 1) {
			go_factored[n + i] = 1; // because it is projective
			}
		else {
			go_factored[n + i] = q - 1;
			}
		}
	if (Mtx->f_semilinear) {
		go_factored[2 * n] = F->e;
		}
	else {
		go_factored[2 * n] = 1;
		}
	D.multiply_up(target_go, go_factored, 2 * n + 1);
	if (f_v) {
		cout << "group ordered factored: ";
		INT_vec_print(cout, go_factored, 2 * n + 1);
		cout << endl;
		cout << "target_go=" << target_go << endl;
		}
	gens = new vector_ge;
	gens->init(A);
	gens->allocate(n - 1 + 1 + 1);
	for (h = 0; h < n - 1 + 1 + 1; h++) {

		INT_vec_zero(data, n * n + 1);
		for (i = 0; i < n; i++) {
			data[i * n + i] = 1;
			}

		if (h < n - 1) {
			// swap basis vector h and h + 1:
			hh = h + 1;
			data[h * n + h] = 0;
			data[hh * n + hh] = 0;
			data[h * n + hh] = 1;
			data[hh * n + h] = 1;
			}
		else if (h == n - 1) {
			data[0] = F->alpha_power(1);
			}
		else if (h == n) {
			if (Mtx->f_semilinear) {
				data[n * n] = 1;
				}
			}
		A->make_element(Elt1, data, 0 /*verbose_level - 1*/);
		if (f_v) {
			cout << "generator " << h << ":" << endl;
			A->element_print_quick(Elt1, cout);
			}
		gens->copy_in(h, Elt1);
		}
	if (f_v) {
		cout << "subspace_orbits::init_monomial_group creating group" << endl;
		}
	S = create_sims_from_generators_randomized(A, 
		gens, TRUE /* f_target_go */, 
		target_go, 0 /*verbose_level - 1*/);
	if (f_v) {
		cout << "subspace_orbits::init_monomial_group after creating group" << endl;
		}
	Strong_gens = new strong_generators;
	Strong_gens->init_from_sims(S, 0);
	if (f_v) {
		cout << "subspace_orbits::init_monomial_group after extracting strong generators" << endl;
		}
	if (f_vv) {
		INT f_print_as_permutation = FALSE;
		INT f_offset = FALSE;
		INT offset = 0;
		INT f_do_it_anyway_even_for_big_degree = FALSE;
		INT f_print_cycles_of_length_one = FALSE;
		
		cout << "strong generators are:" << endl;
		Strong_gens->gens->print(cout, f_print_as_permutation, 
			f_offset, offset, f_do_it_anyway_even_for_big_degree, 
			f_print_cycles_of_length_one);
		}
	f_has_strong_generators = TRUE;
	delete S;
	delete gens;
	FREE_INT(data);
	FREE_INT(go_factored);
	FREE_INT(Elt1);
	sprintf(prefix, "Monomial_%ld_%ld", n, q);
}

void subspace_orbits::init_subfield_structure_action(BYTE *prefix, INT s, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	
	if (f_v) {
		cout << "subspace_orbits::init_subfield_structure_action" << endl;
		cout << "s=" << s << endl;
		}
	f_subfield_structure_action = TRUE;
	subspace_orbits::s = s;


	if (f_v) {
		cout << "subspace_orbits::init_subfield_structure_action before lift_generators_to_subfield_structure" << endl;
		}
	lift_generators_to_subfield_structure(A_linear, P->n + 1, s, P->F, SGens, verbose_level - 1);

	sprintf(prefix, "Subfield_%ld_%ld_%ld", n, q, s);
	
	if (f_v) {
		cout << "subspace_orbits::init_subfield_structure_action done" << endl;
		}
}

void subspace_orbits::init(int argc, const char **argv, 
	INT n, finite_field *F, INT depth, 
	INT f_wedge_action, INT f_PGL2OnConic, INT f_monomial_group, 
	INT f_semilinear, 
	INT f_subfield_structure_action, INT s,
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_basis = TRUE;
	BYTE prefix[1000];

	if (f_v) {
		cout << "subspace_orbits::init" << endl;
		}
	subspace_orbits::n = n;
	subspace_orbits::F = F;
	subspace_orbits::q = F->q;
	subspace_orbits::depth = depth;
	subspace_orbits::f_wedge_action = f_wedge_action;
	subspace_orbits::f_PGL2_OnConic = f_PGL2OnConic;
	subspace_orbits::f_monomial_group = f_monomial_group;
	
	if (f_monomial_group) {
		f_basis = FALSE;
		}
	
	P = new projective_space;

	if (f_v) {
		cout << "subspace_orbits::init before P->init" << endl;
		}
	P->init(n - 1, F, 
		FALSE /* f_init_incidence_structure */, 
		MINIMUM(2, verbose_level));

	if (f_v) {
		cout << "subspace_orbits::init after P->init" << endl;
		}
	

	A_linear = new action;
	A_linear->init_projective_group(n, 
		F, 
		f_semilinear, 
		f_basis, 
		verbose_level - 1);

	Mtx = A_linear->G.matrix_grp;

	vector_space_dimension = n;

	if (f_PGL2_OnConic) {
		init_PGL2q_OnConic(prefix, verbose_level);
		}
	else if (f_wedge_action) {
		init_wedge_action(prefix, verbose_level);
		}
	else if (f_monomial_group) {
		init_monomial_group(prefix, verbose_level);
		}
	else if (f_subfield_structure_action) {
		init_subfield_structure_action(prefix, s, verbose_level);
		}
	else {
		sprintf(prefix, "PGL_%ld_%ld", n, q);
		}

	if (f_v) {
		cout << "subspace_orbits::init vector_space_dimension=" << vector_space_dimension << endl;
		}

	tmp_M = NEW_INT(vector_space_dimension * vector_space_dimension);
	base_cols = NEW_INT(vector_space_dimension);
	Gen = new generator;

	Gen->read_arguments(argc, argv, 0);

	Gen->prefix[0] = 0;
	sprintf(Gen->fname_base, "%s%s", Gen->prefix, prefix);
	
	
	Gen->depth = depth;
}


void subspace_orbits::init2(INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	strong_generators *gens;
	
	if (f_v) {
		cout << "subspace_orbits->init2" << endl;
		}
	if (f_has_strong_generators) {
		// init_group or init_monomial has been called
		if (f_v) {
			cout << "subspace_orbits->init2 initializing with strong generators" << endl;
			}
		gens = Strong_gens;
		}
	else {
		// init_group or init_monomial has not been called
		if (f_v) {
			cout << "subspace_orbits->init2 initializing full group" << endl;
			}
		gens = A_linear->Strong_gens;
		}



	if (f_print_generators) {
		INT f_print_as_permutation = FALSE;
		INT f_offset = TRUE;
		INT offset = 1;
		INT f_do_it_anyway_even_for_big_degree = TRUE;
		INT f_print_cycles_of_length_one = TRUE;
		
		cout << "subspace_orbits->init2 printing generators for the group:" << endl;
		gens->gens->print(cout, f_print_as_permutation, 
			f_offset, offset, 
			f_do_it_anyway_even_for_big_degree, 
			f_print_cycles_of_length_one);
		}



	if (f_PGL2_OnConic) {
		Gen->init(A_linear, A_OnConic, gens, Gen->depth, verbose_level);
		}
	else if (f_wedge_action) {
		Gen->init(A_linear, A_wedge, gens, Gen->depth, verbose_level);
		}
	else if (f_monomial_group) {
		Gen->init(A_linear, A_linear, gens, Gen->depth, verbose_level);
		}
	else if (f_subfield_structure_action) {
		Gen->init(A_linear, A_linear, SGens, Gen->depth, verbose_level);
		}
	else {
		Gen->init(A_linear, A_linear, gens, Gen->depth, verbose_level);
		}

#if 0
	Gen->init_check_func(
		subspace_orbits_test_func, 
		this /* candidate_check_data */);
#endif
	Gen->init_early_test_func(
		subspace_orbits_early_test_func, 
		this /*void *data */,  
		verbose_level);

	//Gen->init_incremental_check_func(
		//check_mindist_incremental, 
		//this /* candidate_check_data */);

	Gen->init_vector_space_action(vector_space_dimension, 
		P->F, 
		subspace_orbits_rank_point_func, 
		subspace_orbits_unrank_point_func, 
		this, 
		verbose_level);
#if 0
	Gen->f_print_function = TRUE;
	Gen->print_function = print_set;
	Gen->print_function_data = this;
#endif	

	INT nb_oracle_nodes = 1000;
	
	if (f_v) {
		cout << "subspace_orbits->init2 before Gen->init_oracle" << endl;
		}
	Gen->init_oracle(nb_oracle_nodes, verbose_level - 1);
	if (f_v) {
		cout << "subspace_orbits->init2 calling Gen->init_root_node" << endl;
		}
	Gen->root[0].init_root_node(Gen, verbose_level - 1);
	
	schreier_depth = Gen->depth;
	f_use_invariant_subset_if_available = FALSE;
	f_implicit_fusion = FALSE;
	f_debug = FALSE;
	if (f_v) {
		cout << "subspace_orbits->init2 done" << endl;
		}
}

void subspace_orbits::read_data_file(INT depth_completed, const BYTE *fname_data_file, 
	INT f_exportmagma, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_recompute_schreier = TRUE;
	
	if (f_v) {
		cout << "subspace_orbits::read_data_file" << endl;
		}
	Gen->read_data_file(depth_completed, 
		fname_data_file, 
		verbose_level - 1);
	
	// ignore the last level: Schreier vectors have not yet been computed
	depth_completed--;
	
	
	if (f_v) {
		cout << "read_data_file: after reading file " << fname_data_file << endl;
		cout << "depth_completed = " << depth_completed << endl;
		}

	if (f_recompute_schreier) {
		if (f_v) {
			cout << "recomputing Schreier vectors" << endl;
			}
		Gen->recreate_schreier_vectors_up_to_level(
			depth_completed, 
			TRUE /* f_compact */, 
			MINIMUM(verbose_level, 1));
		}
	if (f_v) {
		cout << "read_data_file: recreated Schreier vectors" << endl;
		}

#if 0
	INT level;
	INT f, l, i;

	for (level = 0; level <= depth_completed; level++) {
		if (f_v) {
			cout << "level=" << level << endl;
			}
		f = Gen->first_oracle_node_at_level[level];
		l = Gen->first_oracle_node_at_level[level + 1] - f;
		if (f_v) {
			cout << "f=" << f << " l=" << l << endl;
			}
#if 0
		oracle *O;
		

		for (i = 0; i < l; i++) {
			if (f_v) {
				cout << "analyzing node " << i << "/" << l << ":" << endl;
				}
			O = Gen->root + f + i;
			analyze_schreier_vector(O->sv, verbose_level - 1);
			}

#endif

		}
#endif


	if (f_exportmagma) {
		INT level;

		level = depth_completed + 1;

		wedge_product_export_magma(Gen, 
			n, q, vector_space_dimension, 
			level, verbose_level);

		

		} // exportmagma
}

void subspace_orbits::compute_orbits(INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT t0 = os_ticks();
	
	if (f_v) {
		cout << "subspace_orbits::compute_orbits calling generator_main" << endl;
		cout << "A=";
		Gen->A->print_info();
		cout << "A2=";
		Gen->A2->print_info();
		}
	Gen->main(t0, 
		schreier_depth, 
		f_use_invariant_subset_if_available, 
		f_implicit_fusion, 
		f_debug, 
		verbose_level - 1);
	
	INT nb_orbits;
	
	if (f_v) {
		cout << "subspace_orbits::compute_orbits done with generator_main" << endl;
		}
	nb_orbits = Gen->nb_orbits_at_level(depth);
	if (f_v) {
		cout << "subspace_orbits::compute_orbits we found " << nb_orbits << " orbits at depth " << depth<< endl;
		}
}

void subspace_orbits::compute_Kramer_Mesner_matrix(INT t, INT k, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "subspace_orbits::compute_Kramer_Mesner_matrix t=" << t << " k=" << k << ":" << endl;
		}

	// compute Kramer Mesner matrices
	Vector V;
	INT i, j, a;
	
	V.m_l(k);
	for (i = 0; i < k; i++) {
		V[i].change_to_matrix();
		calc_Kramer_Mesner_matrix_neighboring(Gen, i, V[i].as_matrix(), verbose_level - 2);
			// DISCRETA/snakesandladders.C
		
		if (f_v) {
			cout << "matrix level " << i << ":" << endl;
			V[i].as_matrix().print(cout);
			}
		}
	
	matrix Mtk, Mtk_inf;
	
	Mtk_from_MM(V, Mtk, t, k, TRUE, q, verbose_level - 2);
	cout << "M_{" << t << "," << k << "} sup:" << endl;
	Mtk.print(cout);
	
	Mtk_sup_to_inf(Gen, t, k, Mtk, Mtk_inf, 0/*verbose_level - 2*/);	
	cout << "M_{" << t << "," << k << "} inf:" << endl;
	Mtk_inf.print(cout);
	
	cout << endl;
	cout << endl;
	
	INT nb_t_orbits;
	INT nb_k_orbits;
	INT first_t, first_k;
	INT len, rep, size;
	INT set1[1000];
	//INT set2[1000];
	
	first_t = Gen->first_oracle_node_at_level[t];
	first_k = Gen->first_oracle_node_at_level[k];
	nb_t_orbits = Mtk_inf.s_m();
	nb_k_orbits = Mtk_inf.s_n();
	for (i = 0; i < nb_t_orbits; i++) {
		len = Gen->orbit_length_as_INT(i, t);
		//cout << "i=" << i << " len=" << len << endl;
		Gen->get_set(first_t + i, set1, size);
		if (size != t) {
			cout << "size != t" << endl;
			exit(1);
			}
		for (j = 0; j < nb_k_orbits; j++) {
			a = Mtk_inf.s_iji(i, j);
			cout << setw(2) << a << " ";
			}
		cout << "| ";
		cout << setw(3) << i << " " << setw(3) << len << " ";
		if (t == 1) {
			rep = set1[0];
			schreier Schreier;

			Schreier.init(Gen->A2);
			Schreier.init_generators_by_hdl(Gen->root[0].nb_strong_generators, Gen->root[0].hdl_strong_generators, verbose_level - 1);
			Schreier.compute_point_orbit(rep, 0 /* verbose_level */);
			if (Schreier.orbit_len[0] != len) {
				cout << "Schreier.orbit_len[0] != len" << endl;
				exit(1);
				}
			INT *pts;
			INT len1;

			pts = NEW_INT(len);
			Schreier.get_orbit(0 /* orbit_idx */, pts, len1, 0 /* verbose_level */);
			
			//cout << "{";
			INT_vec_print(cout, pts, len);
			//cout << "}";
			FREE_INT(pts);
			}
		cout << endl;
		}
	
}

INT subspace_orbits::test_set(INT len, INT *S, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT ret = TRUE;
	INT i, rk;
	
	if (f_v) {
		cout << "subspace_orbits::test_set" << endl;
		cout << "Testing set ";
		INT_vec_print(cout, S, len);
		cout << endl;
		}
	for (i = 0; i < len; i++) {
		PG_element_unrank_modified(*P->F, tmp_M + i * vector_space_dimension, 1, vector_space_dimension, S[i]);
		}
	if (f_vv) {
		cout << "coordinate matrix:" << endl;
		print_integer_matrix_width(cout, tmp_M, len, vector_space_dimension, vector_space_dimension, P->F->log10_of_q);
		}
	rk = P->F->Gauss_simple(tmp_M, len, vector_space_dimension, base_cols, 0 /*verbose_level - 2*/);
	if (f_v) {
		cout << "the matrix has rank " << rk << endl;
		}
	if (rk < len) {
		ret = FALSE;
		}
	if (ret) {
		if (f_has_extra_test_func) {
			ret = (*extra_test_func)(this, len, S, extra_test_func_data, verbose_level);
			}
		}

	if (ret) {
		if (f_v) {
			cout << "OK" << endl;
			}
		}
	else {
		if (f_v) {
			cout << "not OK" << endl;
			}
		}
	return ret;
}


// ####################################################################################
// global functions:
// ####################################################################################


INT subspace_orbits_rank_point_func(INT *v, void *data)
{
	subspace_orbits *G;
	generator *gen;
	INT rk;
	
	G = (subspace_orbits *) data;
	gen = G->Gen;
	PG_element_rank_modified(*gen->F, v, 1, gen->vector_space_dimension, rk);
	return rk;
}

void subspace_orbits_unrank_point_func(INT *v, INT rk, void *data)
{
	subspace_orbits *G;
	generator *gen;
	
	G = (subspace_orbits *) data;
	gen = G->Gen;
	PG_element_unrank_modified(*gen->F, v, 1, gen->vector_space_dimension, rk);
}

#if 0
INT subspace_orbits_test_func(ostream &ost, INT len, INT *S, void *data, INT verbose_level)
{
	subspace_orbits *SubOrb = (subspace_orbits *) data;
	INT i, rk;
	INT f_OK = TRUE;
	INT f_v = TRUE; //(verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	
	if (f_v) {
		ost << "subspace_orbits_test_func: checking set ";
		print_set(ost, len, S);
		}
	for (i = 0; i < len; i++) {
		PG_element_unrank_modified(*SubOrb->P->F, SubOrb->tmp_M, 1, SubOrb->n, S[i]);
		}
	if (f_vv) {
		ost << "coordinate matrix:" << endl;
		print_integer_matrix_width(ost, SubOrb->tmp_M, len, SubOrb->n, SubOrb->n, SubOrb->P->F->log10_of_q);
		}
	rk = SubOrb->P->F->Gauss_simple(SubOrb->tmp_M, len, SubOrb->n, SubOrb->base_cols, verbose_level - 2);
	if (f_v) {
		ost << "the matrix has rank " << rk << endl;
		}
	if (rk < len) {
		f_OK = FALSE;
		}
	if (f_OK) {
		if (SubOrb->f_has_extra_test_func) {
			f_OK = (*SubOrb->extra_test_func)(SubOrb, len, S, SubOrb->extra_test_func_data, verbose_level);
			}
		}

	if (f_OK) {
		if (f_v) {
			ost << "OK" << endl;
			}
		return TRUE;
		}
	else {
		return FALSE;
		}
}
#endif

void subspace_orbits_early_test_func(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	void *data, INT verbose_level)
{
	//verbose_level = 1;

	subspace_orbits *SubOrb;
	INT f_v = (verbose_level >= 1);
	INT i;

	SubOrb = (subspace_orbits *) data;

	if (f_v) {
		cout << "subspace_orbits_early_test_func" << endl;
		cout << "testing " << nb_candidates << " candidates" << endl;
		}
	nb_good_candidates = 0;
	for (i = 0; i < nb_candidates; i++) {
		S[len] = candidates[i];
		if (SubOrb->test_set(len + 1, S, verbose_level - 1)) {
			good_candidates[nb_good_candidates++] = candidates[i];
			}
		}
	if (f_v) {
		cout << "subspace_orbits_early_test_func" << endl;
		cout << "Out of " << nb_candidates << " candidates, " << nb_good_candidates << " survive" << endl;
		}
}


void lift_generators_to_subfield_structure(action *Aq, INT n, INT s, finite_field *Fq, strong_generators *&SG, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT q, Q, m, t;
	finite_field *FQ;
	action *AQ;
	subfield_structure *S;
	sims *Sims;
	INT *EltQ;
	INT *Eltq;
	INT *Mtx;

	if (f_v) {
		cout << "lift_generators_to_subfield_structure" << endl;
		}
	q = Fq->q;
	Q = i_power_j(q, s);
	m = n / s;
	if (m * s != n) {
		cout << "lift_generators_to_subfield_structure s must divide n" << endl;
		exit(1);
		}
	if (f_v) {
		cout << "lift_generators_to_subfield_structure creating subfield structure" << endl;
		}
	if (f_v) {
		cout << "n=" << n << endl;
		cout << "s=" << s << endl;
		cout << "m=" << m << endl;
		cout << "q=" << q << endl;
		cout << "Q=" << Q << endl;
		}
	FQ = new finite_field;
	FQ->init(Q, 0);

	AQ = new action;
	
	if (f_v) {
		cout << "lift_generators_to_subfield_structure creating AQ" << endl;
		}
	AQ->init_general_linear_group(m, FQ, FALSE /* f_semilinear */, TRUE /* f_basis */, verbose_level - 2);
	if (f_v) {
		cout << "lift_generators_to_subfield_structure creating AQ done" << endl;
		}

	longinteger_object order_GLmQ;
	longinteger_object target_go;
	longinteger_domain D;
	INT r;

	AQ->group_order(order_GLmQ);
	

	cout << "lift_generators_to_subfield_structure order of GL(m,Q) = " << order_GLmQ << endl;
	D.integral_division_by_INT(order_GLmQ, 
		q - 1, target_go, r);
	cout << "lift_generators_to_subfield_structure target_go = " << target_go << endl;

	S = new subfield_structure;
	S->init(FQ, Fq, verbose_level);

	cout << "lift_generators_to_subfield_structure creating subfield structure done" << endl;
		

	vector_ge *gens;
	vector_ge *gens1;
	INT nb_gens;

	gens = AQ->Strong_gens->gens;
	nb_gens = gens->len;

	gens1 = new vector_ge;

	Eltq = NEW_INT(Aq->elt_size_in_INT);
	Mtx = NEW_INT(n * n);

	cout << "lift_generators_to_subfield_structure lifting generators" << endl;
	gens1->init(Aq);
	gens1->allocate(nb_gens);
	for (t = 0; t < nb_gens; t++) {
		cout << "lift_generators_to_subfield_structure " << t << 
" / " << nb_gens << endl;
		EltQ = gens->ith(t);
		S->lift_matrix(EltQ, m, Mtx, 0 /* verbose_level */);
		if (f_v) {
			cout << "lifted matrix:" << endl;
			INT_matrix_print(Mtx, n, n);
			}
		Aq->make_element(Eltq, Mtx, verbose_level - 1);
		if (f_v) {
			cout << "after make_element:" << endl;
			Aq->element_print_quick(Eltq, cout);
			}
		Aq->element_move(Eltq, gens1->ith(t), 0);
		cout << "lift_generators_to_subfield_structure " << t << 
" / " << nb_gens << " done" << endl;
		}

	if (f_v) {
		cout << "lift_generators_to_subfield_structure creating lifted group:" << endl;
		}
	Sims = create_sims_from_generators_with_target_group_order(Aq, 
		gens1, target_go, 0 /* verbose_level */);

#if 0
	Sims = create_sims_from_generators_without_target_group_order(Aq, 
		gens1, MINIMUM(2, verbose_level - 3));
#endif

	if (f_v) {
		cout << "lift_generators_to_subfield_structure creating lifted group done" << endl;
		}

	longinteger_object go;

	Sims->group_order(go);

	if (f_v) {
		cout << "go=" << go << endl;
		}


	SG = new strong_generators;

	SG->init_from_sims(Sims, 0 /* verbose_level */);
	if (f_v) {
		cout << "lift_generators_to_subfield_structure strong generators are:" << endl;
		SG->print_generators();
		}


	delete gens1;
	FREE_INT(Eltq);
	FREE_INT(Mtx);
	delete Sims;
	delete S;
	delete AQ;
	delete FQ;
	if (f_v) {
		cout << "lift_generators_to_subfield_structure done" << endl;
		}

}


void wedge_product_export_magma(generator *Gen, INT n, INT q, INT vector_space_dimension, INT level, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);


	if (f_v) {
		cout << "wedge_product_export_magma" << endl;
		}
	
	//INT level;
	INT *the_set;
	INT *v;
	INT a, i, j, h, fst, len, ii, jj;
	longinteger_object go;
	INT *Elt;
	
	//level = depth_completed + 1;


	the_set = NEW_INT(level);
	v = NEW_INT(vector_space_dimension);
	Elt = NEW_INT(Gen->A->elt_size_in_INT);
	
	fst = Gen->first_oracle_node_at_level[level];
	len = Gen->first_oracle_node_at_level[level + 1] - fst;
	if (f_v) {
		cout << "exporting to magma" << endl;
		cout << "fst=" << fst << " len=" << len << endl;
		}
	oracle *O;
	BYTE fname[1000];

	sprintf(fname, "Wedge_n%ld_q%ld_d%ld.magma", n, q, level);

	{
	ofstream f(fname);

	f << "// file " << fname << endl;
	f << "n := " << n << ";" << endl;
	f << "q := " << q << ";" << endl;
	f << "d := " << level << ";" << endl;
	f << "n2 := " << vector_space_dimension << ";" << endl;
	f << "V := VectorSpace (GF (q), n2);" << endl;
	f << endl;
	f << "/* list of orbit reps */" << endl;
	f << "L := [" << endl;
	f << endl;

	for (i = 0; i < len; i++) {
		O = Gen->root + fst + i;
	
		f << "// orbit rep " << i << endl;
		f << "[" << endl;
		O->store_set_to(Gen, level - 1, the_set);
	 	for (j = 0; j < level; j++) {
			a = the_set[j];
			(*Gen->unrank_point_func)(v, a, Gen->rank_point_data);
			f << "[ ";
			for (h = 0; h < vector_space_dimension; h++) {
				f << v[h];
				if (h < vector_space_dimension - 1)
					f << ", ";
				}
			f << " ]";
			if (j < level - 1) {
				f << "," << endl;
				}
			else {
				f << "]" << endl;
				}
			}
		if (i < len - 1) {
			f << "," << endl << endl;
			}
		else {
			f << endl << "];" << endl << endl;
			}
		} // next i

	f << "// list of orbit lengths " << endl;
	f << "len := \[";
	
	for (i = 0; i < len; i++) {

		if ((i % 20) == 0) {
			f << endl;
			f << "// orbits " << i << " and following:" << endl;
			}

		Gen->orbit_length(i, level, go);
		f << go;
		if (i < len - 1) {
			f << ", ";
			}
		}
	f << "];" << endl << endl;


	f << "// subspaces of vector space " << endl;
	f << "L := [sub< V | L[i]>: i in [1..#L]];" << endl;

	f << "// stabilisers " << endl;
	f << "P := GL(n, q);" << endl;
	f << "E := ExteriorSquare (P);" << endl;


	f << "// base:" << endl;
	f << "BV := VectorSpace (GF (q), n);" << endl;
	f << "B := [ BV | " << endl;
	for (i = 0; i < Gen->A->base_len; i++) {
		a = Gen->A->base[i];
		PG_element_unrank_modified(*Gen->F, v, 1, n, a);
		//(*Gen->unrank_point_func)(v, a, Gen->rank_point_data);
		f << "[ ";
		for (h = 0; h < n; h++) {
			f << v[h];
			if (h < n - 1)
				f << ", ";
			}
        	if (i < Gen->A->base_len - 1)
				f << "], " << endl;
		else f << " ]" << endl;
		}
	f << "];" << endl;
	f << endl;
	f << "P`Base := B;" << endl;

	f << "// list of stabilizer generators" << endl;
	f << "S := [" << endl;
	f << endl;

	for (i = 0; i < len; i++) {
		O = Gen->root + fst + i;
	
		f << "// orbit rep " << i << " has " << O->nb_strong_generators << " strong generators";
		if (O->nb_strong_generators) {
			f << ", transversal lengths: ";
			INT_vec_print(f, O->tl, Gen->A->base_len);
			}
		f << endl;
		f << "[" << endl;

	 	for (j = 0; j < O->nb_strong_generators; j++) {

			Gen->A->element_retrieve(O->hdl_strong_generators[j], Elt, 0);
			
				f << "[";
			//Gen->A->element_print_quick(Elt, f);
			for (ii = 0; ii < n; ii++) {
				f << "[";
				for (jj = 0; jj < n; jj++) {
					a = Elt[ii * n + jj];
					f << a;
					if (jj < n - 1) {
						f << ", ";
						}
					else {
						f << "]";
						}
					}
				if (ii < n - 1) {
					f << "," << endl;
					}
				else {
					f << "]";
					}
				}
			
			if (j < O->nb_strong_generators - 1) {
				f << "," << endl;
				}
			}
			f << "]" << endl;
		if (i < len - 1) {
			f << "," << endl << endl;
			}
		else {
			f << endl << "];" << endl << endl;
			}
		} // next i

         f << endl << 
"T := [sub<GL(n, q) | [&cat (s): s in S[i]]> : i in [1..#S]];" << endl << endl;
	} // file f

	FREE_INT(the_set);
	FREE_INT(v);
	FREE_INT(Elt);

	if (f_v) {
		cout << "written file " << fname << " of size " << file_size(fname) << endl;
		}
	if (f_v) {
		cout << "wedge_product_export_magma done" << endl;
		}
}



