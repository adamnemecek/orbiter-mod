// discreta_global.C
//
// Anton Betten
// Nov 19, 2007

#include "orbiter.h"

void free_global_data()
{
	orthogonal_points_free_global_data();
	longinteger_free_global_data();
}

void the_end(INT t0)
{
	cout << "***************** The End **********************" << endl;
	if (f_memory_debug) {
		cout << "freeing global data" << endl;
		free_global_data();
		registry_dump();
		registry_dump_sorted();

		dump_object_memory();
		}
	time_check(cout, t0);
	cout << endl;
}

void dump_object_memory()
{
	cout << "GALOIS:" << endl;
	if (finite_field::cntr_new)
		cout << "finite_field::cntr_new=" << finite_field::cntr_new << endl;
	if (finite_field::cntr_objects)
		cout << "finite_field::cntr_objects=" << finite_field::cntr_objects << endl;
	if (longinteger_object::cntr_new)
		cout << "longinteger_object::cntr_new=" << longinteger_object::cntr_new << endl;
	if (longinteger_object::cntr_objects)
		cout << "longinteger_object::cntr_objects=" << longinteger_object::cntr_objects << endl;
	if (longinteger_domain::cntr_new)
		cout << "longinteger_domain::cntr_new=" << longinteger_domain::cntr_new << endl;
	if (longinteger_domain::cntr_objects)
		cout << "longinteger_domain::cntr_objects=" << longinteger_domain::cntr_objects << endl;
	if (unipoly_domain::cntr_new)
		cout << "unipoly_domain::cntr_new=" << unipoly_domain::cntr_new << endl;
	if (unipoly_domain::cntr_objects)
		cout << "unipoly_domain::cntr_objects=" << unipoly_domain::cntr_objects << endl;
	if (rank_checker::cntr_new)
		cout << "rank_checker::cntr_new=" << rank_checker::cntr_new << endl;
	if (rank_checker::cntr_objects)
		cout << "rank_checker::cntr_objects=" << rank_checker::cntr_objects << endl;
	if (mp_graphics::cntr_new)
		cout << "mp_graphics::cntr_new=" << mp_graphics::cntr_new << endl;
	if (mp_graphics::cntr_objects)
		cout << "mp_graphics::cntr_objects=" << mp_graphics::cntr_objects << endl;
	if (partitionstack::cntr_new)
		cout << "partitionstack::cntr_new=" << partitionstack::cntr_new << endl;
	if (partitionstack::cntr_objects)
		cout << "partitionstack::cntr_objects=" << partitionstack::cntr_objects << endl;
	if (orthogonal::cntr_new)
		cout << "orthogonal::cntr_new=" << orthogonal::cntr_new << endl;
	if (orthogonal::cntr_objects)
		cout << "orthogonal::cntr_objects=" << orthogonal::cntr_objects << endl;
	
	cout << "ACTION:" << endl;
	if (action::cntr_new)
		cout << "action::cntr_new=" << action::cntr_new << endl;
	if (action::cntr_objects)
		cout << "action::cntr_objects=" << action::cntr_objects << endl;
	if (action_on_sets::cntr_new)
		cout << "action_on_sets::cntr_new=" << action_on_sets::cntr_new << endl;
	if (action_on_sets::cntr_objects)
		cout << "action_on_sets::cntr_objects=" << action_on_sets::cntr_objects << endl;
	if (product_action::cntr_new)
		cout << "product_action::cntr_new=" << product_action::cntr_new << endl;
	if (product_action::cntr_objects)
		cout << "product_action::cntr_objects=" << product_action::cntr_objects << endl;
	if (matrix_group::cntr_new)
		cout << "matrix_group::cntr_new=" << matrix_group::cntr_new << endl;
	if (matrix_group::cntr_objects)
		cout << "matrix_group::cntr_objects=" << matrix_group::cntr_objects << endl;
	if (perm_group::cntr_new)
		cout << "perm_group::cntr_new=" << perm_group::cntr_new << endl;
	if (perm_group::cntr_objects)
		cout << "perm_group::cntr_objects=" << perm_group::cntr_objects << endl;
	if (page_storage::cntr_new)
		cout << "page_storage::cntr_new=" << page_storage::cntr_new << endl;
	if (page_storage::cntr_objects)
		cout << "page_storage::cntr_objects=" << page_storage::cntr_objects << endl;
	if (vector_ge::cntr_new)
		cout << "vector_ge::cntr_new=" << vector_ge::cntr_new << endl;
	if (vector_ge::cntr_objects)
		cout << "vector_ge::cntr_objects=" << vector_ge::cntr_objects << endl;
	if (schreier::cntr_new)
		cout << "schreier::cntr_new=" << schreier::cntr_new << endl;
	if (schreier::cntr_objects)
		cout << "schreier::cntr_objects=" << schreier::cntr_objects << endl;
	if (sims::cntr_new)
		cout << "sims::cntr_new=" << sims::cntr_new << endl;
	if (sims::cntr_objects)
		cout << "sims::cntr_objects=" << sims::cntr_objects << endl;
	if (group::cntr_new)
		cout << "group::cntr_new=" << group::cntr_new << endl;
	if (group::cntr_objects)
		cout << "group::cntr_objects=" << group::cntr_objects << endl;
	cout << endl;

	cout << "vector_ge:" << endl;
	dump_memory_chain(vector_ge::allocated_objects);
}


void the_end_quietly(INT t0)
{
	//cout << "the_end_quietly: freeing global data" << endl;
	free_global_data();
	time_check(cout, t0);
	cout << endl;
}

void calc_Kramer_Mesner_matrix_neighboring(generator *gen, 
	INT level, matrix &M, INT verbose_level)
// we assume that we don't use implicit fusion nodes
{
	INT f_v = (verbose_level >= 1);
	INT f1, f2, f3, l1, l2, i, j, k, I, J, len, ii;
	oracle *O;
	
	if (f_v) {
		cout << "calc_Kramer_Mesner_matrix_neighboring level=" << level << endl;
		}
	f1 = gen->first_oracle_node_at_level[level];
	f2 = gen->first_oracle_node_at_level[level + 1];
	f3 = gen->first_oracle_node_at_level[level + 2];
	l1 = f2 - f1;
	l2 = f3 - f2;
	M.m_mn_n(l1, l2);
	for (i = 0; i < l1; i++) {
		I = f1 + i;
		O = &gen->root[I];
		for (k = 0; k < O->nb_extensions; k++) {
			if (O->E[k].type == 1) {
				len = O->E[k].orbit_len;
				J = O->E[k].data;
				j = J - f2;
				M.s_iji(i, j) += len;
				}
			if (O->E[k].type == 2) {
				// fusion node
				len = O->E[k].orbit_len;

				O->store_set(gen, level - 1);
					// stores a set of size level to gen->S
				gen->S[level] = O->E[k].pt;

				for (ii = 0; ii < level + 1; ii++) {
					gen->set[level + 1][ii] = gen->S[ii];
					}
				
				gen->A->element_one(gen->transporter->ith(level + 1), 0);
				
				J = O->apply_fusion_element(gen, 
					level, I /* current_node */, 
					//0 /* my_node */, 0 /* my_extension */, 0 /* my_coset */, 
					k /* current_extension */, level + 1, 
					FALSE /* f_tolerant */, 
					0/*verbose_level - 2*/);
				if (FALSE) {
					cout << "after apply_fusion_element J=" << J << endl;
					}
#if 0
				//cout << "fusion node:" << endl;
				//INT_vec_print(cout, gen->S, level + 1);
				//cout << endl;
				gen->A->element_retrieve(O->E[k].data, gen->Elt1, 0);
	
				gen->A2->map_a_set(gen->S, gen->S0, level + 1, gen->Elt1, 0);
				//INT_vec_print(cout, gen->S0, level + 1);
				//cout << endl;

				INT_vec_heapsort(gen->S0, level + 1); //INT_vec_sort(level + 1, gen->S0);

				//INT_vec_print(cout, gen->S0, level + 1);
				//cout << endl;

				J = gen->find_oracle_node_for_set(level + 1, gen->S0, 0);
#endif
				j = J - f2;
				M.s_iji(i, j) += len;
				}
			}
		}
	if (f_v) {
		cout << "calc_Kramer_Mesner_matrix_neighboring level=" << level << " done" << endl;
		}
}



void Mtk_from_MM(Vector & MM, matrix & Mtk, INT t, INT k, 
	INT f_subspaces, INT q,  INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	matrix M, N;
	INT i;
	
	if (f_v) {
		cout << "Mtk_from_MM(): t = " << t << ", k = " << k << endl;
		}
	if (k == t) {
		matrix &M1 = MM[t - 1].as_matrix();
		INT n = M1.s_n();
		Mtk.m_mn_n(n, n);
		Mtk.one();
		return;
		}
	M = MM[t].as_matrix();
	for (i = t + 2; i <= k; i++) {
		Mtk_via_Mtr_Mrk(t, i - 1, i, f_subspaces, q,  
			M, MM[i - 1].as_matrix(), N, verbose_level - 1);
		N.swap(M);
		}
	M.swap(Mtk);
	if (f_v) {
		cout << "Mtk_from_MM(): t = " << t << ", k = " << k << " done" << endl;
		}
}

void Mtk_via_Mtr_Mrk(INT t, INT r, INT k, INT f_subspaces, INT q,  
	matrix & Mtr, matrix & Mrk, matrix & Mtk, INT verbose_level)
// Computes $M_{tk}$ via a recursion formula:
// $M_{tk} = {{k - t} \choose {k - r}} \cdot M_{t,r} \cdot M_{r,k}$.
{
	INT f_v = (verbose_level >= 1);
	base s, h;
	INT i, j, m, n;
	
	if (f_v) {
		cout << "Mtk_via_Mtr_Mrk(): t = " << t << ", r = " << r << ", k = " << k << endl;
		}
	Mtk.mult(Mtr, Mrk);
		/* Mtk := (k - t) atop (k - r) * M_t,k */
	if (f_subspaces) {
		longinteger_domain D;
		longinteger_object a;
		
		D.q_binomial(a, k - t, r - t, q, 0/* verbose_level*/);
		s.m_i_i(a.as_INT());
		}
	else {
		Binomial(k - t, k - r, s);
		}
	m = Mtk.s_m();
	n = Mtk.s_n();
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			Mtk[i][j].integral_division_exact(s, h);
			Mtk[i][j] = h;
			}
		}
	if (f_v) {
		cout << "Mtk_via_Mtr_Mrk matrix M_{" << t << "," << k 
			<< "} of format " << m << " x " << n << " computed" << endl;
		}
}

void Mtk_sup_to_inf(generator *gen, 
	INT t, INT k, matrix & Mtk_sup, matrix & Mtk_inf, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT Nt, Nk, i, j, a;
	INT *M_sup;
	INT *M_inf;	
	
	if (f_v) {
		cout << "Mtk_sup_to_inf" << endl;
		cout << "M_{" << t << "," << k << "} sup=" << endl;
		cout << Mtk_sup << endl;
		}
	Nt = Mtk_sup.s_m();
	Nk = Mtk_sup.s_n();
	M_sup = NEW_INT(Nt * Nk);
	M_inf = NEW_INT(Nt * Nk);
	for (i = 0; i < Nt; i++) {
		for (j = 0; j < Nk; j++) {
			M_sup[i * Nk + j] = Mtk_sup.s_iji(i, j);
			}
		}
	gen->Asup_to_Ainf(t, k, M_sup, M_inf, verbose_level);
	Mtk_inf.m_mn_n(Nt, Nk);
	for (i = 0; i < Nt; i++) {
		for (j = 0; j < Nk; j++) {
			a = M_inf[i * Nk + j];
			Mtk_inf.m_iji(i, j, a);
			}
		}
	if (f_v) {
		cout << "Mtk_sup_to_inf" << endl;
		cout << "M_{" << t << "," << k << "} inf=" << endl;
		cout << Mtk_inf << endl;
		}
	
}


