// oracle.C
//
// Anton Betten
// December 27, 2004

#include "orbiter.h"

void oracle::init_root_node(generator *gen, INT verbose_level)
// copies gen->SG0 and gen->transversal_length into the oracle structure using store_strong_generators
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "oracle::init_root_node() initializing root node" << endl;
		}

	freeself();

	node = 0;
	prev = -1;
	sv = NULL;


	if (f_v) {
		cout << "storing strong generators" << endl;
		}
	store_strong_generators(gen, gen->Strong_gens /* *gen->SG0, gen->transversal_length */);
		// stores the strong generators into the oracle structure,
		// copies transversal_length into tl
	if (f_v) {
		cout << "init_root_node done" << endl;
		}

}

void oracle::init_extension_node_prepare_H(generator *gen,
	INT prev, INT prev_ex, INT size,
	group &G, longinteger_object &go_G,
	group &H, longinteger_object &go_H,
	INT pt, INT pt_orbit_len,
	INT verbose_level)
// sets up the group H which is the stabilizer of the point pt in G
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT f_vvv = (verbose_level >= 3);
	//INT f_v5 = (verbose_level >= 5);
	//INT f_v10 = (verbose_level >= 10);
	longinteger_object q, r;
	longinteger_domain D;

	H.init(gen->A);

	oracle *Op = &gen->root[prev];


	if (f_v) {
		cout << "oracle::init_extension_node_prepare_H" << endl;
		}


	if (f_vv) {
		gen->print_level_extension_info(size, prev, prev_ex);
		INT_vec_print(cout, gen->S, size);
		cout << "oracle::init_extension_node_prepare_H computing stabilizer of point " << pt
			<< " (of index " << pt_orbit_len << " in a group of order " << go_G;
		if (G.f_has_sims && !go_G.is_one()) {
			cout << " = ";
			G.S->print_group_order_factored(cout);
			}
		cout << ")" << endl;
		cout << "verbose_level=" << verbose_level << endl;
		}

	//cout << "computing point stabilizer" << endl;
	if (f_vv) {
		gen->print_level_extension_info(size, prev, prev_ex);
		INT_vec_print(cout, gen->S, size);
		cout << "oracle::init_extension_node_prepare_H computing stabilizer of point " << pt << " in group of order " << go_G << endl;
		}

	if (gen->f_on_subspaces) {


		action_on_factor_space AF;
		action A_factor_space;
		INT projected_pt;

		if (gen->f_early_test_func) {
			INT i;

			Op->setup_factor_space_action_with_early_test(gen,
				AF, A_factor_space, size - 1,
				verbose_level - 4);
			for (i = 0; i < AF.nb_cosets; i++) {
				if (AF.preimage(i, 0) == pt) {
					if (f_vv) {
						cout << "oracle::init_extension_node_prepare_H: point pt=" << pt << " is coset " << i << endl;
						}
					break;
					}
				}
			if (i == AF.nb_cosets) {
				cout << "oracle::init_extension_node_prepare_H fatal: could not find the coset corresponding to point " << pt << endl;
				exit(1);
				}
			projected_pt = i;
			}
		else {
			// no early_test_func:

			if (f_vvv) {
				gen->print_level_extension_info(size, prev, prev_ex);
				cout << " oracle::init_extension_node_prepare_H  setting up factor space action:" << endl;
				}
			Op->setup_factor_space_action(gen, AF, A_factor_space, size - 1,
				TRUE /*f_compute_tables*/, verbose_level - 4);
			projected_pt = AF.project(pt, verbose_level - 4);
			}

		if (f_vvv) {
			gen->print_level_extension_info(size, prev, prev_ex);
			cout << " oracle::init_extension_node_prepare_H  pt=" << pt << " projected_pt=" << projected_pt << endl;
			}
		if (projected_pt == -1) {
			cout << "oracle::init_extension_node_prepare_H fatal: projected_pt == -1" << endl;
			exit(1);
			}
		if (f_vvv) {
			gen->print_level_extension_info(size, prev, prev_ex);
			cout << " oracle::init_extension_node_prepare_H calling G.point_stabilizer_with_action" << endl;
			}
		G.point_stabilizer_with_action(&A_factor_space, H, projected_pt, verbose_level - 4);
		if (f_vvv) {
			gen->print_level_extension_info(size, prev, prev_ex);
			cout << " oracle::init_extension_node_prepare_H G.point_stabilizer_with_action done" << endl;
			}
		}
	else {
		// action on sets:

		if (Op->sv) {
			if (f_vvv) {
				gen->print_level_extension_info(size, prev, prev_ex);
				cout << " oracle::init_extension_node_prepare_H setting up restricted action from the previous schreier vector:" << endl;
				}
			action AR;


			if (Op->nb_strong_generators) {
				// if G is non-trivial
				AR.induced_action_by_restriction_on_orbit_with_schreier_vector(*gen->A2,
					FALSE /* f_induce_action */, NULL /* old_G */,
					Op->sv, pt, 0 /*verbose_level - 1*/);
				if (f_vvv) {
					gen->print_level_extension_info(size, prev, prev_ex);
					cout << " oracle::init_extension_node_prepare_H calling G.point_stabilizer_with_action:" << endl;
					}
				G.point_stabilizer_with_action(&AR, H, 0 /*pt */, 0 /*verbose_level - 2*/);
				if (f_vvv) {
					gen->print_level_extension_info(size, prev, prev_ex);
					cout << " oracle::init_extension_node_prepare_H G.point_stabilizer_with_action done" << endl;
					}
				}
			else {
				// do nothing, the stabilizer is trivial (since G is trivial)
				vector_ge stab_gens;
				INT *tl;
				INT i;

				stab_gens.init(gen->A);
				stab_gens.allocate(0);
				tl = NEW_INT(gen->A->base_len);
				for (i = 0; i < gen->A->base_len; i++) {
					tl[i] = 1;
					}

				H.init(gen->A);
				H.init_strong_generators(stab_gens, tl);
				FREE_INT(tl);
				}
			}
		else {
			if (f_vvv) {
				gen->print_level_extension_info(size, prev, prev_ex);
				cout << " previous schreier vector not available:" << endl;
				}
			G.point_stabilizer_with_action(gen->A2, H, pt, 0);
			}
		}
	// now H has strong generators only

	if (f_vv) {
		gen->print_level_extension_info(size, prev, prev_ex);
		cout << "oracle::init_extension_node_prepare_H calling schreier_sims for point stabilizer" << endl;
		}
	H.schreier_sims(0);

	if (f_vv) {
		gen->print_level_extension_info(size, prev, prev_ex);
		cout << "oracle::init_extension_node_prepare_H after schreier_sims for point stabilizer" << endl;
		}

	H.group_order(go_H);
	D.integral_division(go_G, go_H, q, r, 0);
	if (f_vv) {
		gen->print_level_extension_info(size, prev, prev_ex);
		INT_vec_print(cout, gen->S, size);
		cout << "oracle::init_extension_node_prepare_H point stabilizer has order ";
		H.print_group_order(cout);
		//cout << endl;
		cout << ", of index = " << q << " in " << go_G << endl;
		//H.S->print(TRUE);
		}
	if (q.as_INT() != pt_orbit_len) {
		gen->print_level_extension_info(size, prev, prev_ex);
		INT_vec_print(cout, gen->S, size);
		cout << "oracle::init_extension_node_prepare_H: fatal: q != pt_orbit_len" << endl;
		cout << "q = " << q.as_INT() << endl;
		cout << "pt_orbit_len = " << pt_orbit_len << endl;
		exit(1);
		}
	if (f_vv) {
		gen->print_level_extension_info(size, prev, prev_ex);
		INT_vec_print(cout, gen->S, size);
		cout << "oracle::init_extension_node_prepare_H point stabilizer is generated by:";
		INT f_print_as_permutation = FALSE;
		if (/*f_v10 &&*/ gen->A2->degree < 100) {
			f_print_as_permutation = TRUE;
			}
		H.print_strong_generators(cout, f_print_as_permutation);
		}

	if (f_v) {
		cout << "oracle::init_extension_node_prepare_H done" << endl;
		}
}

void oracle::init_extension_node_prepare_G(generator *gen,
	INT prev, INT prev_ex, INT size, group &G, longinteger_object &go_G,
	INT verbose_level)
// sets up the group G using the strong generators that are stored
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT f_vvv = (verbose_level >= 3);

	if (f_v) {
		cout << "oracle::init_extension_node_prepare_G" << endl;
		}
	oracle &Op = gen->root[prev];

	G.init(gen->A);
	if (f_vv) {
		gen->print_level_extension_info(size, prev, prev_ex);
		INT_vec_print(cout, gen->S, size);
		cout << "oracle::init_extension_node_prepare_G calling init_strong_generators_by_hdl" << endl;
		INT_vec_print(cout, Op.hdl_strong_generators, Op.nb_strong_generators);
		cout << endl;
		cout << "verbose_level=" << verbose_level << endl;
		}
	G.init_strong_generators_by_hdl(Op.nb_strong_generators, Op.hdl_strong_generators, Op.tl, 0);
	if (f_vv) {
		gen->print_level_extension_info(size, prev, prev_ex);
		INT_vec_print(cout, gen->S, size);
		cout << "oracle::init_extension_node_prepare_G calling schreier_sims for stabilizer with "
			<< Op.nb_strong_generators << " strong generators" << endl;
		}
	//if (f_vvv) {
		gen->print_level_extension_info(size, prev, prev_ex);
		INT_vec_print(cout, gen->S, size);
		cout << "oracle::init_extension_node_prepare_G the strong generators are:" << endl;
		G.print_strong_generators(cout, FALSE /* f_print_as_permutation */);
		//}
	cout<<">>>>G order before";
	G.print_group_order(cout);
	G.schreier_sims(0);
	G.group_order(go_G);
	cout<<">>>>G order after";
	G.print_group_order(cout);
	if (f_vv) {
		gen->print_level_extension_info(size, prev, prev_ex);
		INT_vec_print(cout, gen->S, size);
		cout << "_{" << go_G << "}, previous stabilizer reconstructed" << endl;
		}

	if (f_v) {
		cout << "oracle::init_extension_node_prepare_G done" << endl;
		}

}

INT oracle::get_level(generator *gen)
{
	INT l;

	l = depth_of_node(gen);
	return l;
}

INT oracle::get_node_in_level(generator *gen)
{
	INT l, n;

	l = depth_of_node(gen);
	n = node - gen->first_oracle_node_at_level[l];
	return n;
}

INT oracle::get_nb_of_live_points()
{
	INT n;

	if (sv == NULL) {
		return -1;
		}
	n = sv[0];
	return n;
}

INT oracle::get_nb_of_orbits_under_stabilizer()
{
	INT nb;

	if (sv == NULL) {
		return -1;
		}
	nb = sv_number_of_orbits(sv);
	return nb;
}

void oracle::get_stabilizer_order(generator *gen, longinteger_object &go)
{
	strong_generators *Strong_gens;

	get_stabilizer_generators(gen,
		Strong_gens, 0 /*verbose_level*/);
	Strong_gens->group_order(go);
}

void oracle::get_stabilizer(generator *gen,
	group &G, longinteger_object &go_G,
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	G.init(gen->A);
	G.init_strong_generators_by_hdl(nb_strong_generators, hdl_strong_generators, tl, 0);
	if (f_v) {
		cout << "oracle::get_stabilizer calling schreier_sims for stabilizer with "
			<< nb_strong_generators << " strong generators" << endl;
		}
	G.schreier_sims(verbose_level - 3);
	G.group_order(go_G);
	if (f_v) {
		cout << "oracle::get_stabilizer stabilizer has order " << go_G << endl;
		}
}

void oracle::get_stabilizer_generators(generator *gen,
	strong_generators *&Strong_gens,
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;

	if (f_v) {
		cout << "oracle::get_stabilizer_generators" << endl;
		cout << "oracle::get_stabilizer_generators nb_strong_generators=" << nb_strong_generators << endl;
		}
	Strong_gens = new strong_generators;
	Strong_gens->init_by_hdl(gen->A, hdl_strong_generators, nb_strong_generators, 0);
	//gens.init_by_hdl(gen->A, hdl_strong_generators, nb_strong_generators);
	if (nb_strong_generators == 0) {
		for (i = 0; i < gen->A->base_len; i++) {
			Strong_gens->tl[i] = 1;
			}
		}
	else {
		for (i = 0; i < gen->A->base_len; i++) {
			Strong_gens->tl[i] = oracle::tl[i];
			}
		}
}

oracle::oracle()
{
	null();
}

oracle::~oracle()
{
	freeself();
}

void oracle::null()
{
	nb_strong_generators = 0;
	hdl_strong_generators = NULL;
	tl = NULL;
	nb_extensions = 0;
	E = NULL;
	sv = NULL;
}

void oracle::freeself()
{
	if (hdl_strong_generators) {
#if 0
		cout << "oracle::freeself() deleting hdl_strong_generators: ";
		INT_vec_print(cout, hdl_strong_generators, nb_strong_generators);
		cout << endl;
		cout << "pointer = ";
		print_pointer_hex(cout, hdl_strong_generators);
		cout << endl;
#endif
		FREE_INT(hdl_strong_generators);
		hdl_strong_generators = NULL;
		//cout << "oracle::freeself() deleting hdl_strong_generators done" << endl;
		}
	if (tl) {
		//cout << "oracle::freeself() deleting tl" << endl;
		FREE_INT(tl);
		tl = NULL;
		}
	if (E) {
		//cout << "oracle::freeself() deleting E" << endl;
		delete [] E;
		E = NULL;
		}
	if (sv) {
		//cout << "oracle::freeself() deleting sv" << endl;
		FREE_INT(sv);
		sv = NULL;
		}
	//cout << "oracle::freeself() finished" << endl;
}

void oracle::oracle_depth_breadth_perm_and_inverse(generator *gen, INT max_depth,
	INT &idx, INT hdl, INT cur_depth, INT *perm, INT *perm_inv)
{
	INT i, nxt;

	perm[idx] = hdl;
	perm_inv[hdl] = idx;
	idx++;
	if (cur_depth == max_depth)
		return;
	for (i = 0; i < nb_extensions; i++) {
		if (E[i].type == EXTENSION_TYPE_EXTENSION) {
			nxt = E[i].data;
			if (nxt >= 0) {
				gen->root[nxt].oracle_depth_breadth_perm_and_inverse(gen,
					max_depth, idx, nxt, cur_depth + 1, perm, perm_inv);
				}
			}
		}
}

INT oracle::find_extension_from_point(generator *gen, INT pt, INT verbose_level)
// a -1 means not found
{
	INT i;

	for (i = 0; i < nb_extensions; i++) {
		if (E[i].pt == pt)
			break;
		}
	if (i == nb_extensions) {
		return -1;
		}
	return i;
}

void oracle::print_extensions(ostream &ost)
{
	INT i;

	ost << "Node " << node << ", the extensions are" << endl;
	if (nb_extensions >= 10) {
		ost << "too many to print" << endl;
		return;
		}
	ost << "i : pt : orbit_len : type : to where" << endl;
	for (i = 0; i < nb_extensions; i++) {
		ost << setw(5) << i << " : "
			<< setw(7) << E[i].pt << " : "
			<< setw(5) << E[i].orbit_len << " : ";

		print_extension_type(ost, E[i].type);
		if (E[i].type == EXTENSION_TYPE_FUSION) {
			ost << " -> (" << E[i].data1 << "," << E[i].data2 << ") hdl=" << E[i].data << endl;
			}
		else if (E[i].type == EXTENSION_TYPE_EXTENSION) {
			ost << " -> " << E[i].data << endl;
			}
		else {
			ost << setw(5) << E[i].data << endl;
			}
		if (E[i].type >= NB_EXTENSION_TYPES) {
			ost << "E[i].type >= NB_EXTENSION_TYPES" << endl;
			exit(1);
			}
		}
	cout << "done with node " << node << endl;
}

void oracle::store_strong_generators(generator *gen, strong_generators *Strong_gens /* vector_ge &SG, INT *tl*/)
{
	INT i;

	nb_strong_generators = Strong_gens->gens->len;
	if (nb_strong_generators == 0) {
		hdl_strong_generators = NULL;
		oracle::tl = NULL;
		}
	else {
		hdl_strong_generators = NEW_INT(nb_strong_generators);
		oracle::tl = NEW_INT(gen->A->base_len);
		for (i = 0; i < nb_strong_generators; i++) {
			hdl_strong_generators[i] = gen->A->element_store(Strong_gens->gens->ith(i), FALSE);
			}
		for (i = 0; i < gen->A->base_len; i++) {
			oracle::tl[i] = Strong_gens->tl[i];
			}
		}
}

void oracle::log_current_node_without_group(generator *gen, INT s, ostream &f, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	longinteger_object go;
	INT i;

	if (f_v) {
		cout << "oracle::log_current_node_without_group" << endl;
		}
	store_set_to(gen, s - 1, gen->set0);

	if (f_v) {
		f << "# ***** orbit ***** " << node - gen->first_oracle_node_at_level[s] << " "<< endl;
		}
	f << s << " ";
	for (i = 0; i < s; i++) {
		f << gen->set0[i] << " ";
		}
	f << endl;

	if (f_v) {
		f << "# BEGINCOMMENT" << endl;
		if (gen->f_print_function) {
			(*gen->print_function)(f, s, gen->set0, gen->print_function_data);
			}

		f << "# ENDCOMMENT" << endl;
		}
}

void oracle::log_current_node(generator *gen, INT s, ostream &f, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	longinteger_object go, rk;
	INT i;

	if (f_v) {
		cout << "oracle::log_current_node node=" << node << " s=" << s << endl;
		}
	store_set_to(gen, s - 1, gen->set0);
	if (f_v) {
		cout << "oracle::log_current_node node=" << node << " after store_set_to" << endl;
		}

	if (f_v) {
		f << "# ***** orbit ***** " << node - gen->first_oracle_node_at_level[s] << " "<< endl;
		}
	f << s << " ";
	for (i = 0; i < s; i++) {
		f << gen->set0[i] << " ";
		}

	if (nb_strong_generators == 0) {
		f << " 1" << endl;
		if (f_v) {
			cout << "oracle::log_current_node node=" << node << " done" << endl;
			}
		return;
		}


	if (f_v) {
		cout << "oracle::log_current_node node=" << node << " creating group" << endl;
		}

	group G;

	G.init(gen->A);
	G.init_strong_generators_by_hdl(nb_strong_generators, hdl_strong_generators, tl, FALSE);

#if 0
	if (node == 26) {
		cout << "oracle::log_current_node node=26 with " << nb_strong_generators << " strong generators" << endl;
		strong_generators *Strong_gens;

		Strong_gens = new strong_generators;
		Strong_gens->init_by_hdl(gen->A, hdl_strong_generators, nb_strong_generators, verbose_level);
		Strong_gens->print_generators();
		delete Strong_gens;
		cout << "tl=" << endl;
		INT_vec_print(cout, tl, gen->A->base_len);
		cout << endl;

		//G.schreier_sims(verbose_level);
		}
#endif

	if (f_v) {
		cout << "oracle::log_current_node node=" << node << " before schreier_sims" << endl;
		}
	G.schreier_sims(0);
	if (f_v) {
		cout << "oracle::log_current_node node=" << node << " after schreier_sims" << endl;
		}
	G.group_order(go);
	if (f_v) {
		cout << "oracle::log_current_node node=" << node << " group order = " << go << endl;
		}
	//if (f_v) {
		//cout << "oracle::log_current_node() stabilizer of order " << go << " reconstructed" << endl;
		//}
	if (go.is_one()) {
		go.print_not_scientific(f);
		f << endl;
		//f << go << endl;
		}
	else {
		G.code_ascii(FALSE);
		go.print_not_scientific(f);
		f << " " << G.ascii_coding << endl;
		//f << go << " " << G.ascii_coding << endl;
		}

	if (FALSE /* f_v */) {
		f << "# BEGINCOMMENT" << endl;
		if (gen->f_print_function) {
			(*gen->print_function)(f, s, gen->set0, gen->print_function_data);
			}

		if (!go.is_one()) {
			if (f_v) {
				cout << "oracle::log_current_node node=" << node << " printing generators" << endl;
				}
			G.require_strong_generators();
			f << "tl: ";
			for (i = 0; i < G.A->base_len; i++)
				f << G.tl[i] << " ";
			f << endl;
			f << G.SG->len << " strong generators by rank: " << endl;
			for (i = 0; i < G.SG->len; i++) {
				f << i << " : " << endl;

				G.A->element_print(G.SG->ith(i), f);
				f << endl;
				G.A->element_print_as_permutation(G.SG->ith(i), f);
				f << endl;

#if 0
				G.A->element_rank(rk, G.SG->ith(i), 0);
				f << "\"" << rk << "\", ";
				f << endl;
#endif
				}
			//for (i = 0; i < G.SG->len; i++) {
				//}
			}
		f << "# ENDCOMMENT" << endl;
		}
}

void oracle::log_current_node_after_applying_group_element(generator *gen, INT s, ostream &f, INT hdl, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	longinteger_object go;
	INT i;
	INT *S;
	INT *Elt;
	INT *Elt_inv;
	INT *Elt1;
	INT *Elt2;

	S = NEW_INT(s);
	Elt = NEW_INT(gen->A->elt_size_in_INT);
	Elt_inv = NEW_INT(gen->A->elt_size_in_INT);
	Elt1 = NEW_INT(gen->A->elt_size_in_INT);
	Elt2 = NEW_INT(gen->A->elt_size_in_INT);

	store_set_to(gen, s - 1, gen->set0);
	gen->A->element_retrieve(hdl, Elt, 0);
	//gen->A->element_print(Elt, cout);
	gen->A->element_invert(Elt, Elt_inv, 0);
	for (i = 0; i < s; i++) {
		S[i] = Elt[gen->set0[i]];
		}

	if (f_v) {
		f << "# ***** orbit ***** " << node - gen->first_oracle_node_at_level[s] << " "<< endl;
		}
	f << s << " ";
	for (i = 0; i < s; i++) {
		f << S[i] << " ";
		}
	group G;

	G.init(gen->A);
	G.init_strong_generators_by_hdl(nb_strong_generators, hdl_strong_generators, tl, FALSE);
	G.schreier_sims(0);
	G.group_order(go);
	//if (f_v) {
		//cout << "oracle::log_current_node() stabilizer of order " << go << " reconstructed" << endl;
		//}
	if (go.is_one()) {
		f << go << endl;
		}
	else {
		G.code_ascii(FALSE);
		f << go << " " << G.ascii_coding << endl;
		}

	if (f_v) {
		if (gen->f_print_function) {
			(*gen->print_function)(f, s, S, gen->print_function_data);
			}
		if (!go.is_one()) {
			G.require_strong_generators();
			f << "# ";
			for (i = 0; i < G.A->base_len; i++)
				f << G.tl[i] << " ";
			f << endl;
			for (i = 0; i < G.SG->len; i++) {
				f << "# ";
				//G.A->element_print(G.SG->ith(i), f);
				G.A->element_mult(Elt_inv, G.SG->ith(i), Elt1, FALSE);
				G.A->element_mult(Elt1, Elt, Elt2, FALSE);
				G.A->element_print(Elt2, f);
				//f << endl;
				}
			}
		}
	FREE_INT(S);
	FREE_INT(Elt);
	FREE_INT(Elt_inv);
	FREE_INT(Elt1);
	FREE_INT(Elt2);
}

void oracle::log_current_node_with_candidates(generator *gen, INT lvl, ostream &f, INT verbose_level)
{
	//INT f_v = (verbose_level >= 1);
	INT i;

	store_set_to(gen, lvl - 1, gen->set0);

	f << lvl << " ";
	for (i = 0; i < lvl; i++) {
		f << gen->set0[i] << " ";
		}
	f << -1 << " ";

	INT n;
	INT *subset;
	INT *candidates = NULL;
	INT nb_candidates = 0;
	INT f_subset_is_allocated;

	if (!downstep_get_invariant_subset(
		gen,
		lvl,
		n, subset, f_subset_is_allocated,
		verbose_level)) {
		cout << "oracle::log_current_node_with_candidates downstep_get_invariant_subset returns FALSE" << endl;
		exit(1);
		}
	candidates = NEW_INT(n);

	downstep_apply_early_test(gen, lvl,
		n, subset,
		candidates, nb_candidates,
		verbose_level - 2);
	f << nb_candidates << " ";
	for (i = 0; i < nb_candidates; i++) {
		f << candidates[i] << " ";
		}
	f << -1 << endl;
	if (f_subset_is_allocated) {
		FREE_INT(subset);
		}
	FREE_INT(candidates);
}


INT oracle::depth_of_node(generator *gen)
{
	if (prev == -1) {
		return 0;
		}
	else {
		return gen->root[prev].depth_of_node(gen) + 1;
		}
}

void oracle::store_set(generator *gen, INT i)
// stores a set of size i + 1 to gen->S[]
{
	if (i < 0)
		return;
	gen->S[i] = pt;
	if (i >= 0) {
		if (prev == -1) {
			cout << "store_set prev == -1" << endl;
			exit(1);
			}
		gen->root[prev].store_set(gen, i - 1);
		}
}

void oracle::store_set_with_verbose_level(generator *gen, INT i, INT verbose_level)
// stores a set of size i + 1 to gen->S[]
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "oracle::store_set_with_verbose_level node=" << node << " prev=" << prev << " pt=" << pt << " i=" << i << endl;
		}
	if (i < 0)
		return;
	gen->S[i] = pt;
	if (i >= 0) {
		if (prev == -1) {
			cout << "store_set prev == -1" << endl;
			exit(1);
			}
		gen->root[prev].store_set(gen, i - 1);
		}
}

void oracle::store_set_to(generator *gen, INT i, INT *to)
// stores a set of size i + 1 to 'to'
{
	if (i < 0)
		return;
	to[i] = pt;
	if (i >= 0) {
		if (prev == -1) {
			cout << "store_set_to prev == -1" << endl;
			exit(1);
			}
		gen->root[prev].store_set_to(gen, i - 1, to);
		}
}

void oracle::store_set_to(generator *gen, INT *to)
{
	store_set_to(gen, depth_of_node(gen), to);
}

INT oracle::check_node_and_set_consistency(generator *gen, INT i, INT *set)
{
	if (i < 0)
		return TRUE;
	if (set[i] != pt) {
		cout << "check_node_and_set_consistency() inconsistent" << endl;
		return FALSE;
		}
	if (i >= 0) {
		if (prev == -1) {
			cout << "check_node_and_set_consistency prev == -1" << endl;
			exit(1);
			}
		gen->root[prev].check_node_and_set_consistency(gen, i - 1, set);
		}
	return TRUE;
}

void oracle::print_set_verbose(generator *gen)
{
	INT depth;
	INT *set;

	//cout << "oracle::print_set_verbose" << endl;
	depth = depth_of_node(gen);
	print_set(gen);
	cout << endl;


	set = NEW_INT(depth);
	store_set_to(gen, depth - 1, set /* gen->S0 */);
	if (gen->f_print_function) {
		(*gen->print_function)(cout, depth, set /* gen->S0 */, gen->print_function_data);
		}
	FREE_INT(set);
	//cout << "oracle::print_set_verbose done" << endl;
}

void oracle::print_set(generator *gen)
{
	INT depth, size, i;
	longinteger_object go;
	longinteger_domain D;
	INT *set;

	depth = depth_of_node(gen);
	//cout << "oracle::print_set depth = " << depth << endl;
	size = depth;
	set = NEW_INT(size);
	store_set_to(gen, depth - 1, set /*gen->S0*/);
	INT_set_print(cout, set /*gen->S0*/, size);
	if (nb_strong_generators == 0) {
		cout << "_1";
		}
	else {
		D.multiply_up(go, tl, gen->A->base_len);
		cout << "_{";
		for (i = 0; i < gen->A->base_len; i++) {
			cout << tl[i];
			if (i < gen->A->base_len - 1)
				cout << " * ";
			}
		cout << " = " << go << "}";
		}
	FREE_INT(set);
}

void oracle::print_node(generator *gen)
{
	INT depth;
	INT *set;
	//INT i, depth, node2, len;
	//INT *orbit;

	//orbit = NEW_INT(gen->A->degree);
	depth = depth_of_node(gen);
	cout << "Node " << node << " at depth " << depth << ", prev=" << prev << endl;
	print_set(gen);
	cout << endl;
	//cout << "pt=" << pt << endl;
	cout << "nb_strong_generators=" << nb_strong_generators << endl;
	cout << "nb_extensions=" << nb_extensions << endl;

	set = NEW_INT(depth);
	store_set_to(gen, depth - 1, set /*gen->S0*/);

	if (gen->f_print_function) {
		(*gen->print_function)(cout, depth, set /* gen->S0 */, gen->print_function_data);
		}

	FREE_INT(set);
	print_extensions(gen);

#if 0
	for (i = 0; i < nb_extensions; i++) {
		cout << setw(3) << i << " : " << setw(7) << E[i].pt << " : " << setw(5) << E[i].orbit_len << " : ";
		len = gen->A->compute_orbit_of_point_generators_by_handle(
			nb_strong_generators, hdl_strong_generators, E[i].pt, orbit, 0);
		if (len != E[i].orbit_len) {
			cout << "oracle::print_node len != E[i].orbit_len" << endl;
			cout << "len = " << len << endl;
			cout << "E[i].orbit_len = " << E[i].orbit_len << endl;
			}
		INT_vec_heapsort(orbit, len); // INT_vec_sort(len, orbit);
		if (E[i].type == EXTENSION_TYPE_UNPROCESSED) {
			cout << "unprocessed";
			}
		else if (E[i].type == EXTENSION_TYPE_EXTENSION) {
			cout << "extension to node " << E[i].data;
			}
		else if (E[i].type == EXTENSION_TYPE_FUSION) {
			//cout << "fusion node from ";
			gen->A->element_retrieve(E[i].data, gen->Elt1, FALSE);
			store_set(gen, depth - 1);
			gen->S[depth] = E[i].pt;
			//INT_vec_print(cout, gen->S, depth + 1);
			//cout << " to ";
			gen->A->map_a_set(gen->S, gen->set[0], depth + 1, gen->Elt1, 0);
			//INT_vec_print(cout, gen->set[0], depth + 1);
			INT_vec_heapsort(gen->set[0], depth + 1); // INT_vec_sort(depth + 1, gen->set[0]);
			//cout << " = ";
			//INT_vec_print(cout, gen->set[0], depth + 1);
			node2 = gen->find_oracle_node_for_set(depth + 1, gen->set[0], 0 /* f_tolerant */, 0);
			//cout << node2;
			cout << "fusion to node " << node2;
			}
		else if (E[i].type == EXTENSION_TYPE_PROCESSING) {
			cout << "currently processing";
			}
		cout << " : ";
		INT_vec_print(cout, orbit, len);
		cout << endl;
		}
	FREE_INT(orbit);
#endif
}

void oracle::print_extensions(generator *gen)
{
	//INT i, depth, /*node2,*/ len;
	INT depth;
	INT *orbit;

	depth = depth_of_node(gen);
	cout << "oracle::print_extensions node=" << node << " at depth " << depth << " degree=" << gen->A2->degree << endl;
	print_extensions(cout);
	orbit = NEW_INT(gen->A2->degree);

	if (nb_extensions >= 10) {
		cout << "too many to print (nb_extensions=" << nb_extensions << ")" << endl;
		return;
		}

#if 0
	for (i = 0; i < nb_extensions; i++) {
		cout << setw(3) << i << " : " << setw(7) << E[i].pt << " : " << setw(5) << E[i].orbit_len << " : ";

#if 0
		cout << "before gen->A->compute_orbit_of_point_generators_by_handle nb_strong_generators=" << nb_strong_generators << endl;
#endif

		if (FALSE) {
			len = gen->A2->compute_orbit_of_point_generators_by_handle(
				nb_strong_generators, hdl_strong_generators, E[i].pt, orbit, 0);
			cout << "orbit of length " << len << endl;
			if (len != E[i].orbit_len) {
				cout << "oracle::print_extensions len != E[i].orbit_len" << endl;
				cout << "len = " << len << endl;
				cout << "E[i].orbit_len = " << E[i].orbit_len << endl;
				}
			INT_vec_heapsort(orbit, len);
			}
		if (E[i].type == EXTENSION_TYPE_UNPROCESSED) {
			cout << "unprocessed";
			}
		else if (E[i].type == EXTENSION_TYPE_EXTENSION) {
			cout << "extension to node " << E[i].data;
			}
		else if (E[i].type == EXTENSION_TYPE_FUSION) {
			cout << "fusion node from " << endl;
			store_set_with_verbose_level(gen, depth - 1, 1);
			gen->S[depth] = E[i].pt;
			INT_vec_print(cout, gen->S, depth + 1);

			cout << "fusion handle=" << E[i].data << endl;
			gen->A->element_retrieve(E[i].data, gen->Elt1, FALSE);
			cout << "fusion element:" << endl;
			gen->A2->element_print_quick(gen->Elt1, cout);

			cout << " to " << E[i].data1 << "/" << E[i].data2 << endl;
#if 0
			gen->A2->map_a_set(gen->S, gen->set[0], depth + 1, gen->Elt1, 0);
			INT_vec_print(cout, gen->set[0], depth + 1);
			cout << endl;
			INT_vec_heapsort(gen->set[0], depth + 1);
			cout << " = " << endl;
			INT_vec_print(cout, gen->set[0], depth + 1);
			cout << endl;
			node2 = gen->find_oracle_node_for_set(depth + 1, gen->set[0], 0 /* f_tolerant */, 0);
			cout << "Which is node " << node2 << endl;
			cout << "fusion to node " << node2 << endl;
#endif

			}
		else if (E[i].type == EXTENSION_TYPE_PROCESSING) {
			cout << "currently processing";
			}
		cout << " : ";


		//INT_vec_print(cout, orbit, len);
		cout << endl;
		}
#endif

	FREE_INT(orbit);
}

void oracle::reconstruct_extensions_from_sv(generator *gen, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	INT n, nb, i, j, a, idx;
	INT *pts;
	INT *prev;
	INT *ancestor;
	INT *depth;
	INT *orbit_reps;


	if (f_v) {
		cout << "oracle::reconstruct_extensions_from_sv" << endl;
		}
	n = sv[0];
	nb = sv_number_of_orbits(sv);
	if (f_v) {
		cout << "n=" << n << " nb=" << nb << endl;
		}
	pts = sv + 1;
	prev = pts + n;

	ancestor = NEW_INT(n);
	depth = NEW_INT(n);
	orbit_reps = NEW_INT(nb);
	for (i = 0; i < n; i++) {
		depth[i] = -1;
		ancestor[i] = -1;
		}
	for (i = 0; i < n; i++) {
		schreier_vector_determine_depth_recursion(n, pts, prev, depth, ancestor, i);
		}

	nb_extensions = nb;
	E = new extension[nb];
	for (i = 0; i < nb; i++) {
		E[i].orbit_len = 0;
		E[i].type = EXTENSION_TYPE_UNPROCESSED;
		}
	j = 0;
	for (i = 0; i < n; i++) {
		if (prev[i] == -1) {
			E[j].pt = pts[i];
			orbit_reps[j] = pts[i];
			j++;
			}
		}
	for (i = 0; i < n; i++) {
		a = ancestor[i];
		if (!INT_vec_search(orbit_reps, nb, a, idx)) {
			cout << "oracle::reconstruct_extensions_from_sv  did not find orbit rep" << endl;
			exit(1);
			}
		E[idx].orbit_len++;
		}

	FREE_INT(ancestor);
	FREE_INT(depth);
	FREE_INT(orbit_reps);
}
