// generator.C
//
// Anton Betten
// December 29, 2003

#include "orbiter.h"


INT generator::nb_orbits_at_level(INT level)
{
	INT f, l;

	f = first_oracle_node_at_level[level];
	l = first_oracle_node_at_level[level + 1] - f;
	return l;
}

INT generator::poset_structure_is_contained(INT *set1, INT sz1, INT *set2, INT sz2, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT f_contained;
	INT i, rk1, rk2;

	if (f_v) {
		cout << "poset_structure_is_contained" << endl;
		}
	if (f_vv) {
		cout << "set1: ";
		INT_vec_print(cout, set1, sz1);
		cout << " ; ";
		cout << "set2: ";
		INT_vec_print(cout, set2, sz2);
		cout << endl;
		}
	if (sz1 > sz2) {
		f_contained = FALSE;
		}
	else {
		if (f_on_subspaces) {
			INT *B1, *B2;

			B1 = NEW_INT(sz1 * vector_space_dimension);
			B2 = NEW_INT((sz1 + sz2) * vector_space_dimension);

			for (i = 0; i < sz1; i++) {
				(*unrank_point_func)(B1 + i * vector_space_dimension, set1[i], rank_point_data);
				}
			for (i = 0; i < sz2; i++) {
				(*unrank_point_func)(B2 + i * vector_space_dimension, set2[i], rank_point_data);
				}

			rk1 = F->Gauss_easy(B1, sz1, vector_space_dimension);
			if (rk1 != sz1) {
				cout << "poset_structure_is_contained rk1 != sz1" << endl;
				exit(1);
				}
			
			rk2 = F->Gauss_easy(B2, sz2, vector_space_dimension);
			if (rk2 != sz2) {
				cout << "poset_structure_is_contained rk2 != sz2" << endl;
				exit(1);
				}
			INT_vec_copy(B1, B2 + sz2 * vector_space_dimension, sz1 * vector_space_dimension);
			rk2 = F->Gauss_easy(B2, sz1 + sz2, vector_space_dimension);
			if (rk2 > sz2) {
				f_contained = FALSE;
				}
			else {
				f_contained = TRUE;
				}

			FREE_INT(B1);
			FREE_INT(B2);
			}
		else {
			f_contained = INT_vec_sort_and_test_if_contained(set1, sz1, set2, sz2);
			}
		}
	return f_contained;
}

void generator::print_progress_by_extension(INT size, INT cur, INT prev, INT cur_ex, INT nb_ext_cur, INT nb_fuse_cur)
{
	double progress;
	
	
	progress = level_progress(size);
		
	print_level_info(size + 1, prev);
	cout << " **** Upstep extension " << cur_ex << " / " << root[prev].nb_extensions << " with " 
		<< nb_ext_cur << " new orbits and " 
		<< nb_fuse_cur << " fusion nodes. We now have " 
		<< cur - first_oracle_node_at_level[size + 1] 
		<< " nodes at level " << size + 1;
		cout << ", ";
	print_progress(size + 1, progress);
	print_progress_by_level(size + 1);
}

void generator::print_progress(INT size, INT cur, INT prev, INT nb_ext_cur, INT nb_fuse_cur)
{
	double progress;
	
	
	progress = level_progress(size);
		
	print_level_info(size + 1, prev);
	cout << " **** Upstep finished with " 
		<< nb_ext_cur << " new orbits and " 
		<< nb_fuse_cur << " fusion nodes. We now have " 
		<< cur - first_oracle_node_at_level[size + 1] 
		<< " nodes at level " << size + 1;
		cout << ", ";
	print_progress(size + 1, progress);
	print_progress_by_level(size + 1);
}

void generator::print_progress(INT lvl, double progress)
{
	double progress0;
	INT progress1, progress2;
		
	progress0 = progress * 100.;
	progress2 = (INT) (progress0 * 100.);
	progress1 = progress2 / 100;
	progress2 = progress2 % 100;
	cout << "progress: " << progress1 << "." << setw(2) << progress2 << " % " << endl;
}

void generator::print_progress_by_level(INT lvl)
{
	INT i;
	
	for (i = 0; i < lvl; i++) {
		//remaining = nb_extension_nodes_at_level_total[i] 
		//	- nb_extension_nodes_at_level[i] - nb_fusion_nodes_at_level[i];
		cout << setw(5) << i << " : " << setw(10) << nb_extension_nodes_at_level[i] << " : " 
			<< setw(10) << nb_fusion_nodes_at_level[i] << " : " 
			<< setw(10) << nb_extension_nodes_at_level_total[i] << " : " 
			<< setw(10) << nb_unprocessed_nodes_at_level[i];
		cout << endl;
		}
	//print_statistic_on_callbacks();
}

void generator::print_orbit_numbers(INT depth)
{
	INT nb_nodes, j;
	
	nb_nodes = nb_orbits_at_level(depth);
	cout << "##################################################################################################" << endl;
	print_problem_label();
	cout << "Found " << nb_nodes << " orbits at depth " << depth << endl;
	for (j = 0; j <= depth; j++) {
		cout << j << " : " << nb_orbits_at_level(j) << " orbits" << endl;
		}
	cout << "total: " << first_oracle_node_at_level[depth + 1] << endl;
	//gen->print_statistic_on_callbacks();
	compute_and_print_automorphism_group_orders(depth, cout);
}


void generator::print_statistic_on_callbacks_naked()
{
	cout << A->nb_times_image_of_called - nb_times_image_of_called0 << "/";
	cout << A->nb_times_mult_called - nb_times_mult_called0 << "/";
	cout << A->nb_times_invert_called - nb_times_invert_called0 << "/";
	cout << A->nb_times_retrieve_called - nb_times_retrieve_called0 << "/";
	cout << A->nb_times_store_called - nb_times_store_called0;
}

void generator::print_statistic_on_callbacks()
{
	cout << "# of calls to image_of/mult/invert/retrieve/store: ";
	print_statistic_on_callbacks_naked();
	cout << endl;
}

void generator::get_set_by_level(INT level, INT node, INT *set)
{
	INT n, size;
	
	n = first_oracle_node_at_level[level] + node;
	size = root[n].depth_of_node(this);
	if (size != level) {
		cout << "generator::get_set_by_level size != level" << endl;
		exit(1);
		}
	root[n].store_set_to(this, size - 1, set);
}

void generator::get_set(INT node, INT *set, INT &size)
{
	size = root[node].depth_of_node(this);
	root[node].store_set_to(this, size - 1, set);
}

void generator::print_set_verbose(INT node)
{
	root[node].print_set_verbose(this);
}

void generator::print_set(INT node)
{
	root[node].print_set(this);
}

INT generator::find_oracle_node_for_set(INT len, INT *set, INT f_tolerant, INT verbose_level)
// finds the node that represents s_0,...,s_{len - 1}
{
	INT f_v = (verbose_level >= 1);
	INT ret;
	
	if (f_v) {
		cout << "generator::find_oracle_node_for_set ";
		INT_vec_print(cout, set, len);
		cout << endl;
		}
	if (f_starter) {
		INT i, j, h;
		if (len < starter_size) {
			cout << "generator::find_oracle_node_for_set len < starter_size" << endl;
			cout << "len=" << len << endl;
			exit(1);
			}
		for (i = 0; i < starter_size; i++) {
			for (j = i; j < len; j++) {
				if (set[j] == starter[i]) {
					if (f_v) {
						cout << "found " << i << "-th element of the starter which is " << starter[i] << " at position " << j << endl;
						}
					break;
					}
				}
			if (j == len) {
				cout << "generator::find_oracle_node_for_set did not find " << i << "-th element of the starter" << endl;
				}
			for (h = j; h > i; h--) {
				set[h] = set[h - 1];
				}
			set[i] = starter[i];
			}
		INT from = starter_size;
		INT node = starter_size;
		ret = find_oracle_node_for_set_basic(from, node, len, set, f_tolerant, verbose_level);
		}
	else {
		INT from = 0;
		INT node = 0;
		ret = find_oracle_node_for_set_basic(from, node, len, set, f_tolerant, verbose_level);
		}
	if (ret == -1) {
		if (f_tolerant) {
			if (f_v) {
				cout << "generator::find_oracle_node_for_set ";
				INT_vec_print(cout, set, len);
				cout << " extension not found, we are tolerant, returnning -1" << endl;
				}
			return -1;
			}
		else {
			cout << "generator::find_oracle_node_for_set we should not be here" << endl;
			exit(1);
			}
		}
	return ret;
	
}

INT generator::find_oracle_node_for_set_basic(INT from, INT node, INT len, INT *set, INT f_tolerant, INT verbose_level)
{
	INT i, j, pt;
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 1);

	if (f_vv) {
		cout << "generator::find_oracle_node_for_set_basic looking for set ";
		INT_vec_print(cout, set, len);
		cout << endl;
		cout << "node=" << node << endl;
		cout << "from=" << from << endl;
		cout << "len=" << len << endl;
		cout << "f_tolerant=" << f_tolerant << endl;
		}
	for (i = from; i < len; i++) {
		pt = set[i];
		if (f_vv) {
			cout << "pt=" << pt << endl;
			cout << "calling root[node].find_extension_from_point" << endl;
			}
		j = root[node].find_extension_from_point(this, pt, FALSE);
		if (j == -1) {
			if (f_v) {
				cout << "generator::find_oracle_node_for_set_basic  depth " << i << " no extension for point " << pt << " found" << endl;
				}
			if (f_tolerant) {
				if (f_v) {
					cout << "generator::find_oracle_node_for_set_basic  since we are tolerant, we return -1" << endl;
					}
				return -1;
				}
			else {
				cout << "generator::find_oracle_node_for_set_basic failure in find_extension_from_point" << endl;
				INT_vec_print(cout, set, len);
				cout << endl;
				cout << "node=" << node << endl;
				cout << "from=" << from << endl;
				cout << "i=" << i << endl;
				cout << "pt=" << pt << endl;
				root[node].print_extensions(this);
				exit(1);
				}
			}
		if (root[node].E[j].pt != pt) {
			cout << "generator::find_oracle_node_for_set() root[node].E[j].pt != pt" << endl;
			exit(1);
			}
		if (root[node].E[j].type != EXTENSION_TYPE_EXTENSION && 
			root[node].E[j].type != EXTENSION_TYPE_PROCESSING) {
			cout << "generator::find_oracle_node_for_set() root[node].E[j].type != EXTENSION_TYPE_EXTENSION" << endl;
			cout << "root[node].E[j].type=" << root[node].E[j].type << " = ";
			print_extension_type(cout, root[node].E[j].type);
			cout << endl;
			cout << "generator::find_oracle_node_for_set_basic looking for set ";
			INT_vec_print(cout, set, len);
			cout << endl;
			cout << "node=" << node << endl;
			cout << "from=" << from << endl;
			cout << "i=" << i << endl;
			cout << "node=" << node << endl;
			cout << "f_tolerant=" << f_tolerant << endl;
			cout << "node=" << node << endl;
			cout << "pt=" << pt << endl;
			cout << "j=" << j << endl;
			exit(1);
			}
		node = root[node].E[j].data;
		if (f_v) {
			cout << "depth " << i << " extension " << j << " new node " << node << endl;
			}
		}
	return node;
}

void generator::oracle_depth_breadth_perm_and_inverse(INT max_depth, 
	INT *&perm, INT *&perm_inv, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT idx = 0;
	INT N;

	if (f_v) {
		cout << "generator::oracle_depth_breadth_perm_and_inverse" << endl;
		cout << "max_depth = " << max_depth << endl;
		}

	N = first_oracle_node_at_level[max_depth + 1];
	if (f_v) {
		cout << "N = first_oracle_node_at_level[max_depth + 1] = " << N << endl;
		}
	
	perm = NEW_INT(N);
	perm_inv = NEW_INT(N);
	
	if (f_v) {
		cout << "calling root->oracle_depth_breadth_perm_and_inverse" << endl;
		}
	root->oracle_depth_breadth_perm_and_inverse(this, max_depth, idx, 0, 0, perm, perm_inv);
	
}

INT generator::count_extension_nodes_at_level(INT lvl)
{
	INT prev;

	nb_extension_nodes_at_level_total[lvl] = 0;
	for (prev = first_oracle_node_at_level[lvl]; 
			prev < first_oracle_node_at_level[lvl + 1]; prev++) {
			
		nb_extension_nodes_at_level_total[lvl] += root[prev].nb_extensions;
		
		}
	nb_unprocessed_nodes_at_level[lvl] = nb_extension_nodes_at_level_total[lvl];
	nb_fusion_nodes_at_level[lvl] = 0;
	nb_extension_nodes_at_level[lvl] = 0;
	return nb_extension_nodes_at_level_total[lvl];
}

double generator::level_progress(INT lvl)
{
	return ((double)(nb_fusion_nodes_at_level[lvl] + nb_extension_nodes_at_level[lvl])) / 
			(double) nb_extension_nodes_at_level_total[lvl];
}



void generator::count_automorphism_group_orders(INT lvl, INT &nb_agos, 
	longinteger_object *&agos, INT *&multiplicities, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, l, j, c, h, f_added;
	longinteger_object ago;
	longinteger_object *tmp_agos;
	INT *tmp_multiplicities;
	longinteger_domain D;
	
	l = nb_orbits_at_level(lvl);
	if (f_v) {
		cout << "collecting the automorphism group orders of " << l << " orbits" << endl;
		}
	nb_agos = 0;
	agos = NULL;
	multiplicities = NULL;
	for (i = 0; i < l; i++) {
		get_stabilizer_order(lvl, i, ago);
		f_added = FALSE;
		for (j = 0; j < nb_agos; j++) {
			c = D.compare_unsigned(ago, agos[j]);
			//cout << "comparing " << ago << " with " << agos[j] << " yields " << c << endl;
			if (c >= 0) {
				if (c == 0) {
					multiplicities[j]++;
					}
				else {
					tmp_agos = agos;
					tmp_multiplicities = multiplicities;
					agos = new longinteger_object[nb_agos + 1];
					multiplicities = NEW_INT(nb_agos + 1);
					for (h = 0; h < j; h++) {
						tmp_agos[h].swap_with(agos[h]);
						multiplicities[h] = tmp_multiplicities[h];
						}
					ago.swap_with(agos[j]);
					multiplicities[j] = 1;
					for (h = j; h < nb_agos; h++) {
						tmp_agos[h].swap_with(agos[h + 1]);
						multiplicities[h + 1] = tmp_multiplicities[h];
						}
					nb_agos++;
					if (tmp_agos) {
						delete [] tmp_agos;
						FREE_INT(tmp_multiplicities);
						}
					}
				f_added = TRUE;
				break;
				}
			}
		if (!f_added) {
			// add at the end (including the case that the list is empty)
			tmp_agos = agos;
			tmp_multiplicities = multiplicities;
			agos = new longinteger_object[nb_agos + 1];
			multiplicities = NEW_INT(nb_agos + 1);
			for (h = 0; h < nb_agos; h++) {
				tmp_agos[h].swap_with(agos[h]);
				multiplicities[h] = tmp_multiplicities[h];
				}
			ago.swap_with(agos[nb_agos]);
			multiplicities[nb_agos] = 1;
			nb_agos++;
			if (tmp_agos) {
				delete [] tmp_agos;
				FREE_INT(tmp_multiplicities);
				}
			}
		}
}

void generator::compute_and_print_automorphism_group_orders(INT lvl, ostream &ost)
{

	INT j, nb_agos;
	longinteger_object *agos;
	INT *multiplicities;
	INT N, r, h;
	longinteger_object S, S1, Q;
	longinteger_domain D;
	
	count_automorphism_group_orders(lvl, nb_agos, agos, multiplicities, FALSE);
	S.create(0);
	N = 0;
	for (j = 0; j < nb_agos; j++) {
		N += multiplicities[j];
		for (h = 0; h < multiplicities[j]; h++) {
			D.add(S, agos[j], S1);
			S1.assign_to(S);
			}
		}
	D.integral_division_by_INT(S, N, Q, r);
	

	ost << "(";
	for (j = 0; j < nb_agos; j++) {
		ost << agos[j];
		if (multiplicities[j] == 1) {
			}
		else if (multiplicities[j] >= 10) {
			ost << "^{" << multiplicities[j] << "}";
			}
		else  {
			ost << "^" << multiplicities[j];
			}
		if (j < nb_agos - 1) {
			ost << ", ";
			}
		}
	ost << ") average is " << Q << " + " << r << " / " << N << endl;
	if (nb_agos) {
		delete [] agos;
		FREE_INT(multiplicities);
		}

}

void generator::stabilizer_order(INT node, longinteger_object &go)
{
	if (root[node].nb_strong_generators) {
		go.create_product(A->base_len, root[node].tl);
		}
	else {
		go.create(1);
		}
}

INT generator::check_the_set(INT len, INT *S, INT verbose_level)
// used by lookahead_first_at_level
{
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "generator::check_the_set" << endl;
		}
	if (!f_candidate_check_func) {
		return TRUE;
		}
	
	if (f_v) {
		cout << "generator::check_the_set checking set: ";
		INT_set_print(cout, S, len);
		cout << endl;
		}
	if (!(*candidate_check_func)(len, S, candidate_check_data, verbose_level - 1)) {
		if (f_v) {
			cout << "the set is not accepted" << endl;
			}
		return FALSE;
		}
	if (f_v) {
		cout << "the set is accepted" << endl;
		}
	return TRUE;
}

INT generator::check_the_set_incrementally(INT len, INT *S, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	
	if (!f_candidate_incremental_check_func) {
		return check_the_set(len, S, verbose_level);
		}
	
	if (f_v) {
		cout << "generator::check_the_set_incrementally checking set: ";
		INT_set_print(cout, S, len);
		cout << endl;
		}
	if (!(*candidate_incremental_check_func)(len, S, candidate_incremental_check_data, verbose_level - 1)) {
		if (f_v) {
			cout << "the set is not accepted" << endl;
			}
		return FALSE;
		}
	if (f_v) {
		cout << "the set is accepted" << endl;
		}
	return TRUE;
}



void generator::orbit_length(INT node, INT level, longinteger_object &len)
// uses generator::go for the group order
{
	longinteger_domain D;
	longinteger_object stab_order, quo, rem;

	get_stabilizer_order(level, node, stab_order);
	D.integral_division(go, stab_order, len, rem, 0);
	if (!rem.is_zero()) {
		cout << "generator::orbit_length stabilizer order does not divide group order" << endl;
		exit(1);
		}
}

INT generator::orbit_length_as_INT(INT node, INT level)
{
	longinteger_object len;

	orbit_length(node, level, len);
	return len.as_INT();
	
}


void generator::print_representatives_at_level(INT lvl)
{
	INT i, f, l;
	
	f = first_oracle_node_at_level[lvl];
	l = nb_orbits_at_level(lvl);
	cout << "The " << l << " representatives at level " << lvl << " are:" << endl;
	for (i = 0; i < l; i++) {
		cout << i << " / " << l << " : ";
		root[f + i].print_set(this);
		cout << endl;
		}
}

void generator::print_problem_label()
{
	if (problem_label[0]) {
		cout << problem_label << " ";
		}
}

void generator::print_level_info(INT i, INT prev)
{
	INT t1, dt;
	
	t1 = os_ticks();
	//cout << "generator::print_level_info t0=" << t0 << endl;
	//cout << "generator::print_level_info t1=" << t1 << endl;
	dt = t1 - t0;
	//cout << "generator::print_level_info dt=" << dt << endl;
	
	cout << "Time ";
	time_check_delta(cout, dt);
	print_problem_label();
	cout << " : Level " << i - 1 << " Node " << prev << " = " 
		<< prev - first_oracle_node_at_level[i - 1] 
		<< " / " 
		<< nb_orbits_at_level(i - 1)
		<< " : ";
}

void generator::print_level_extension_info(INT i, 
	INT prev, INT cur_extension)
{
#if 0
	INT t1, dt;
	
	t1 = os_ticks();
	dt = t1 - t0;
	
	cout << "Time ";
	time_check_delta(cout, dt);
	print_problem_label();
#endif
	cout << "Level " << i - 1 << " Node " << prev << " = " 
		<< prev - first_oracle_node_at_level[i - 1] 
		<< " / " 
		<< nb_orbits_at_level(i - 1)
		<< " Extension " << cur_extension 
		<< " / " 
		<< root[prev].nb_extensions 
		<< " : ";
}

void generator::print_level_extension_coset_info(INT i, 
	INT prev, INT cur_extension, INT coset, INT nb_cosets)
{
#if 0
	INT t1, dt;
	
	t1 = os_ticks();
	dt = t1 - t0;
	
	cout << "Time ";
	time_check_delta(cout, dt);
	print_problem_label();
#endif
	cout << "Level " << i - 1 << " Node " << prev << " = " 
		<< prev - first_oracle_node_at_level[i - 1] 
		<< " / " 
		<< nb_orbits_at_level(i - 1)
		<< " Extension " << cur_extension 
		<< " / " 
		<< root[prev].nb_extensions 
		<< " : " 
		<< "Coset " << coset << " / " << nb_cosets << " : ";
}

void generator::recreate_schreier_vectors_up_to_level(INT lvl, INT f_compact, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;
	
	if (f_v) {
		cout << "generator::recreate_schreier_vectors_up_to_level creating Schreier vectors up to level " << lvl << endl;
		}
	for (i = 0; i <= lvl; i++) {
		if (f_v) {
			cout << "generator::recreate_schreier_vectors_up_to_level creating Schreier vectors at level " << i << endl;
			}
		recreate_schreier_vectors_at_level(i, f_compact, verbose_level);
		}
}

void generator::recreate_schreier_vectors_at_level(INT i, INT f_compact, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT f, cur, l, prev, u;
	INT f_recreate_extensions = FALSE;
	INT f_dont_keep_sv = FALSE;

	f = first_oracle_node_at_level[i];
	cur = first_oracle_node_at_level[i + 1];
	l = cur - f;

	if (f_v) {
		cout << "creating Schreier vectors at depth " << i << " for " << l << " orbits" << endl;
		}
	if (f_v) {
		cout << "generator::recreate_schreier_vectors_at_level Testing if a schreier vector file exists" << endl;
		}
	if (test_sv_level_file_binary(i, fname_base)) {

		if (f_v) {
			cout << "generator::recreate_schreier_vectors_at_level Yes, a schreier vector file exists. We will read this file" << endl;
			}

		read_sv_level_file_binary(i, fname_base, FALSE, 0, 0, 
			f_recreate_extensions, f_dont_keep_sv, 
			verbose_level);
		if (f_v) {
			cout << "read Schreier vectors at depth " << i << " from file" << endl;
			}
		return;
		}


	if (f_v) {
		cout << "generator::recreate_schreier_vectors_at_level No, a schreier vector file does not exists. We will create such a file now" << endl;
		}



	for (u = 0; u < l; u++) {
			
		prev = f + u;
			
		if (f_v && !f_vv) {
			cout << ".";
			if (((u + 1) % 50) == 0) {
				cout << "; " << u + 1 << " / " << l << endl;
				}
			if (((u + 1) % 1000) == 0)
				cout << " " << u + 1 << endl;
			}
		else if (f_vv) {
			cout << "generator::recreate_schreier_vectors_at_level " 
				<< i << " node " << u << " / " << l << endl;
			}
			
		root[prev].compute_schreier_vector(this, i, f_compact, verbose_level - 1);
		}
	write_sv_level_file_binary(i, fname_base, FALSE, 0, 0, verbose_level);
	if (f_v) {
		cout << "generator::recreate_schreier_vectors_at_level Written a file with Schreier vectors at depth " << i << endl;
		}
	if (f_v) {
		cout << endl;
		}
		
}

void generator::print_node(INT node)
{
	cout << "generator::print_node node " << node << ":" << endl;
	root[node].print_node(this);
}

void generator::print_tree()
{
	INT i;
	
	cout << "generator::print_tree nb_oracle_nodes_used=" << nb_oracle_nodes_used << endl;
	for (i = 0; i < nb_oracle_nodes_used; i++) {
		print_node(i);
		}
}

void generator::get_table_of_nodes(INT *&Table, INT &nb_rows, INT &nb_cols, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;
	
	if (f_v) {
		cout << "generator::get_table_of_nodes nb_oracle_nodes_used=" << nb_oracle_nodes_used << endl;
		}
	nb_rows = nb_oracle_nodes_used;
	nb_cols = 6;
	Table = NEW_INT(nb_oracle_nodes_used * nb_cols);
		
	for (i = 0; i < nb_oracle_nodes_used; i++) {

		if (f_v) {
			cout << "generator::get_table_of_nodes node " << i << " / " << nb_oracle_nodes_used << endl;
			}

		Table[i * nb_cols + 0] = root[i].get_level(this);
		Table[i * nb_cols + 1] = root[i].get_node_in_level(this);
		Table[i * nb_cols + 2] = root[i].pt;

		longinteger_object go;
			
		root[i].get_stabilizer_order(this, go);
		Table[i * nb_cols + 3] = go.as_INT();
		Table[i * nb_cols + 4] = root[i].get_nb_of_live_points();
		Table[i * nb_cols + 5] = root[i].get_nb_of_orbits_under_stabilizer();
		}
	if (f_v) {
		cout << "generator::get_table_of_nodes done" << endl;
		}
}

INT generator::count_live_points(INT level, INT node_local, INT f_compact, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT node;
	INT *osv;
	INT nb_points;
	INT *pt_list;
	INT *pt_list2;
	INT a, i, j;
	
	node = first_oracle_node_at_level[level] + node_local;
	if (root[node].sv == NULL) {
		root[node].compute_schreier_vector(this, 
			level, f_compact, verbose_level - 2);
		}
	osv = root[node].sv;
	nb_points = osv[0];
	pt_list = osv + 1;
	if (f_v) {
		cout << "pt_list of size " << nb_points << " : ";
		INT_vec_print(cout, pt_list, nb_points);
		cout << endl;
		}

	pt_list2 = NEW_INT(nb_points);
	a = root[node].pt;
	j = 0;
	for (i = 0; i < nb_points; i++) {
		if (pt_list[i] > a) {
			pt_list2[j++] = pt_list[i]; // why is this needed?
			}
		}
	nb_points = j;
	FREE_INT(pt_list2);
	return nb_points;
}

void generator::find_automorphism_group_of_order(INT level, INT order)
{
	INT nb_nodes, node, i, j, elt_order;
	longinteger_object ago;
	INT set[300];
	
	nb_nodes = nb_orbits_at_level(level);
	for (i = 0; i < nb_nodes; i++) {
		node = first_oracle_node_at_level[level] + i;
		if (root[node].nb_strong_generators == 0) {
			ago.create(1);
			}
		else {
			ago.create_product(A->base_len, root[node].tl);
			}
		if (ago.as_INT() == order) {
			cout << "found a node whose automorphism group is order " << order << endl;
			cout << "the node is # " << i << " at level " << level << endl;
			get_set(first_oracle_node_at_level[level] + i, set, level);
			INT_vec_print(cout, set, level);
			cout << endl;
			
			strong_generators *Strong_gens;
			
			get_stabilizer_generators(Strong_gens,  
				level, i, 0  /* verbose_level */);
				
			for (j = 0; j < Strong_gens->gens->len; j++) {
				elt_order = A->element_order(Strong_gens->gens->ith(j));
				cout << "generator " << j << " of order" << elt_order << ":" << endl;
				if (order == elt_order) {
					cout << "CYCLIC" << endl;
					}
				A->element_print(Strong_gens->gens->ith(j), cout);
				A->element_print_as_permutation(Strong_gens->gens->ith(j), cout);
				}
			delete Strong_gens;
			}
		}
}

void generator::get_stabilizer_order(INT level, INT orbit_at_level, longinteger_object &go)
{
	oracle *O;
	INT nd;

	nd = first_oracle_node_at_level[level] + orbit_at_level;
	O = root + nd;

	if (O->nb_strong_generators == 0) {
		go.create(1);
		}
	else {
		longinteger_domain D;

		D.multiply_up(go, O->tl, A->base_len);
		}
}

void generator::get_stabilizer_group(group *&G,  
	INT level, INT orbit_at_level, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	oracle *O;
	INT node;

	if (f_v) {
		cout << "generator::get_stabilizer_group level=" << level << " orbit_at_level=" << orbit_at_level << endl;
		}
	G = new group;
	node = first_oracle_node_at_level[level] + orbit_at_level;
	O = root + node;

	G->init(A);
	if (f_vv) {
		cout << "generator::get_stabilizer_group before G->init_strong_generators_by_hdl" << endl;
		}
	G->init_strong_generators_by_hdl(O->nb_strong_generators, O->hdl_strong_generators, O->tl, FALSE);
	G->schreier_sims(0);
	
	if (f_v) {
		cout << "generator::get_stabilizer_group level=" << level << " orbit_at_level=" << orbit_at_level << " done" << endl;
		}
}

void generator::get_stabilizer_generators(strong_generators *&gens,  
	INT level, INT orbit_at_level, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "generator::get_stabilizer_generators level=" << level << " orbit_at_level=" << orbit_at_level << endl;
		}
	group *G;

	get_stabilizer_group(G, level, orbit_at_level, verbose_level - 1);

	gens = new strong_generators;

	gens->init_from_sims(G->S, 0 /* verbose_level */);
	delete G;
	if (f_v) {
		cout << "generator::get_stabilizer_generators level=" << level << " orbit_at_level=" << orbit_at_level << " done" << endl;
		}

}

void generator::change_extension_type(INT level, INT node, INT cur_ext, INT type, INT verbose_level)
{
	if (type == EXTENSION_TYPE_EXTENSION) {
		// extension node
		if (root[node].E[cur_ext].type != EXTENSION_TYPE_UNPROCESSED && 
			root[node].E[cur_ext].type != EXTENSION_TYPE_PROCESSING) {
			cout << "generator::change_extension_type trying to install extension node, fatal: root[node].E[cur_ext].type != EXTENSION_TYPE_UNPROCESSED && root[node].E[cur_ext].type != EXTENSION_TYPE_PROCESSING" << endl;
			cout << "root[node].ext[cur_ext].type=" << root[node].E[cur_ext].type << endl;
			exit(1);
			}
		nb_extension_nodes_at_level[level]++;
		nb_unprocessed_nodes_at_level[level]--;
		root[node].E[cur_ext].type = EXTENSION_TYPE_EXTENSION;
		}
	else if (type == EXTENSION_TYPE_FUSION) {
		// fusion
		if (root[node].E[cur_ext].type != EXTENSION_TYPE_UNPROCESSED) {
			cout << "generator::change_extension_type trying to install fusion node, fatal: root[node].E[cur_ext].type != EXTENSION_TYPE_UNPROCESSED" << endl;
			cout << "root[node].ext[cur_ext].type=" << root[node].E[cur_ext].type << endl;
			exit(1);
			}
		nb_fusion_nodes_at_level[level]++;
		nb_unprocessed_nodes_at_level[level]--;
		root[node].E[cur_ext].type = EXTENSION_TYPE_FUSION;
		}
}

void generator::orbit_element_unrank(INT depth, INT orbit_idx, INT rank, INT *set, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT *Elt1;
	INT *Elt2;
	INT *the_set;
	oracle *O;


	if (f_v) {
		cout << "generator::orbit_element_unrank depth=" << depth << " orbit_idx=" << orbit_idx << " rank=" << rank << endl;
		}

	Elt1 = NEW_INT(A->elt_size_in_INT);
	Elt2 = NEW_INT(A->elt_size_in_INT);
	the_set = NEW_INT(depth);
	
	O = &root[first_oracle_node_at_level[depth] + orbit_idx];
	coset_unrank(depth, orbit_idx, rank, Elt1, 0/*verbose_level*/);

	A->element_invert(Elt1, Elt2, 0);
	O->store_set_to(this, depth - 1, the_set);
	A2->map_a_set(the_set, set, depth, Elt2, 0/*verbose_level*/);

	FREE_INT(the_set);
	FREE_INT(Elt1);
	FREE_INT(Elt2);
	if (f_v) {
		cout << "generator::orbit_element_unrank ";
		INT_vec_print(cout, set, depth);
		cout << endl;
		}
}

void generator::orbit_element_rank(INT depth, INT &orbit_idx, INT &rank, INT *set, 
	INT f_implicit_fusion, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT *Elt1;
	INT *the_set;
	INT *canonical_set;
	INT i;


	if (f_v) {
		cout << "generator::orbit_element_rank depth=" << depth << " ";
		INT_vec_print(cout, set, depth);
		cout << endl;
		}

	Elt1 = NEW_INT(A->elt_size_in_INT);
	the_set = NEW_INT(depth);
	canonical_set = NEW_INT(depth);
	for (i = 0; i < depth; i++) {
		the_set[i] = set[i];
		}
	
	orbit_idx = trace_set(the_set, depth, depth, 
		canonical_set, Elt1, 
		f_implicit_fusion, verbose_level - 3);

	// now Elt1 is the transporter element that moves 
	// the given set to the orbit representative

	if (f_vv) {
		cout << "generator::orbit_element_rank after trace_set, orbit_idx = " << orbit_idx << endl;
		cout << "transporter:" << endl;
		A->element_print_quick(Elt1, cout);
		cout << "as permutation:" << endl;
		A2->element_print_as_permutation(Elt1, cout);
		}
	if (f_v) {
		cout << "calling coset_rank" << endl;
		}
	rank = coset_rank(depth, orbit_idx, Elt1, verbose_level);
	if (f_v) {
		cout << "after coset_rank, rank=" << rank << endl;
		}
		
	FREE_INT(Elt1);
	FREE_INT(the_set);
	FREE_INT(canonical_set);
	if (f_v) {
		cout << "generator::orbit_element_rank orbit_idx=" << orbit_idx << " rank=" << rank << endl;
		}
}

void generator::coset_unrank(INT depth, INT orbit_idx, INT rank, INT *Elt, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT *the_set;
	group *G1, *G2;
	INT *Elt_gk;
	longinteger_object G_order, U_order;
	oracle *O1, *O2;

	if (f_v) {
		cout << "generator::coset_unrank depth=" << depth << " orbit_idx=" << orbit_idx << endl;
		cout << "action A:" << endl;
		A->print_info();
		cout << "action A2:" << endl;
		A2->print_info();
		}

	O1 = &root[0];
	O2 = &root[first_oracle_node_at_level[depth] + orbit_idx];


	
	G1 = new group;
	G2 = new group;
	the_set = NEW_INT(depth);
	Elt_gk = NEW_INT(A->elt_size_in_INT);
	
	O2->store_set_to(this, depth - 1, the_set);
	
	if (f_v) {
		cout << "the set representing orbit " << orbit_idx 
			<< " at level " << depth << " is ";
		INT_vec_print(cout, the_set, depth);
		cout << endl;
		}
	
	O1->get_stabilizer(this, *G1, G_order, verbose_level - 2);
	O2->get_stabilizer(this, *G2, U_order, verbose_level - 2);


	A->coset_unrank(G1->S, G2->S, rank, Elt, verbose_level);

	delete G1;
	delete G2;
	FREE_INT(the_set);
	FREE_INT(Elt_gk);

}

INT generator::coset_rank(INT depth, INT orbit_idx, INT *Elt, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT rank;
	INT *the_set;
	group *G1, *G2;
	INT *Elt_gk;
	longinteger_object G_order, U_order;
	oracle *O1, *O2;

	if (f_v) {
		cout << "generator::coset_rank depth=" << depth << " orbit_idx=" << orbit_idx << endl;
		cout << "action A:" << endl;
		A->print_info();
		cout << "action A2:" << endl;
		A2->print_info();
		}

	O1 = &root[0];
	O2 = &root[first_oracle_node_at_level[depth] + orbit_idx];


	
	G1 = new group;
	G2 = new group;
	the_set = NEW_INT(depth);
	Elt_gk = NEW_INT(A->elt_size_in_INT);
	
	O2->store_set_to(this, depth - 1, the_set);
	
	if (f_v) {
		cout << "the set representing orbit " << orbit_idx 
			<< " at level " << depth << " is ";
		INT_vec_print(cout, the_set, depth);
		cout << endl;
		}
	
	O1->get_stabilizer(this, *G1, G_order, verbose_level - 2);
	O2->get_stabilizer(this, *G2, U_order, verbose_level - 2);


	rank = A->coset_rank(G1->S, G2->S, Elt, verbose_level);

	delete G1;
	delete G2;
	FREE_INT(the_set);
	FREE_INT(Elt_gk);
	
	return rank;
}

void generator::list_all_orbits_at_level(INT depth)
{
	INT l, i;

	l = nb_orbits_at_level(depth);

	cout << "listing all orbits at depth " << depth << ":" << endl;
	for (i = 0; i < l; i++) {
		list_whole_orbit(depth, i);
		}
}

void generator::list_whole_orbit(INT depth, INT orbit_idx)
{
	INT *set;
	INT rank, len;
	strong_generators *Strong_gens;
	
	set = NEW_INT(depth);
	len = orbit_length_as_INT(orbit_idx, depth);
	cout << "orbit " << orbit_idx << " (=node " << first_oracle_node_at_level[depth] + orbit_idx << ") at depth " << depth << " has length " << len << ":" << endl;

	get_stabilizer_generators(Strong_gens,  
		depth, orbit_idx, 0 /* verbose_level*/);
	cout << "The stabilizer is generated by:" << endl;
	Strong_gens->print_generators();

	INT max_len;
	if (len > 10) {
		max_len = 10;
		}
	else {
		max_len = len;
		}
	for (rank = 0; rank < max_len; rank++) {
		orbit_element_unrank(depth, orbit_idx, rank, set, 0 /* verbose_level */);
		cout << setw(5) << rank << " : ";
		INT_set_print(cout, set, depth);
		cout << endl;
		}	

	FREE_INT(set);
	delete Strong_gens;
}

void generator::print_extensions_at_level(ostream &ost, INT lvl)
{
	INT i, node;
	INT fst, len;
	oracle *O;
	
	ost << "extensions at level " << lvl << ":" << endl;
	fst = first_oracle_node_at_level[lvl];
	len = nb_orbits_at_level(lvl);
	ost << "there are " << len << " nodes at level " << lvl << ":" << endl;
	for (i = 0; i < len; i++) {
		node = fst + i;
		O = root + node;
		ost << "Node " << i << " / " << len << " = " << node << ":" << endl;
		O->print_extensions(ost);
		}
}

void generator::map_to_canonical_k_subset(INT *the_set, INT set_size, INT subset_size, INT subset_rk, 
	INT *reduced_set, INT *transporter, INT &local_idx, INT verbose_level)
// fills reduced_set[set_size - subset_size], transporter and local_idx
// local_idx is the index of the orbit that the subset belongs to 
// (in the list of orbit of subsets of size subset_size)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);

	if (f_v) {
		cout << "generator::map_to_canonical_k_subset" << endl;
		}
	INT *our_set;
	INT *subset;
	INT *canonical_subset;
	INT *Elt1;
	INT f, i, j, k, idx;
	INT reduced_set_size;
	INT f_implicit_fusion = TRUE;
	
	our_set = NEW_INT(set_size);
	subset = NEW_INT(set_size);
	canonical_subset = NEW_INT(set_size);
	Elt1 = NEW_INT(A->elt_size_in_INT);
	reduced_set_size = set_size - subset_size;

	// unrank the k-subset and its complement to our_set[set_size]:
	unrank_k_subset(subset_rk, our_set, set_size, subset_size);
	j = 0;
	k = 0;
	for (i = 0; i < set_size; i++) {
		if (j < subset_size && our_set[j] == i) {
			j++;
			continue;
			}
		our_set[subset_size + k] = i;
		k++;
		}
	for (i = 0; i < set_size; i++) {
		subset[i] = the_set[our_set[i]];
		set[0][i] = subset[i];
		}
	
	A->element_one(generator::transporter->ith(0), FALSE);


	// trace the subset:
	f = first_oracle_node_at_level[subset_size];
	
	if (f_vv) {
		cout << "generator::map_to_canonical_k_subset before trace_set" << endl;
		}
	local_idx = trace_set(subset, set_size, subset_size, 
		canonical_subset, Elt1, 
		f_implicit_fusion, verbose_level - 3);

	idx = local_idx + f;

	if (f_vv) {
		cout << "generator::map_to_canonical_k_subset after trace_set local_idx=" << local_idx << endl;
		}
	if (FALSE) {
		cout << "the transporter is" << endl;
		A->element_print(Elt1, cout);
		cout << endl;
		}
	A->element_move(Elt1, transporter, FALSE);
	for (i = 0; i < reduced_set_size; i++) {
		reduced_set[i] = canonical_subset[subset_size + i];
		}
	if (FALSE) {
		cout << "generator::map_to_canonical_k_subset reduced set = ";
		INT_vec_print(cout, reduced_set, reduced_set_size);
		cout << endl;
		}
	FREE_INT(Elt1);
	FREE_INT(our_set);
	FREE_INT(subset);
	FREE_INT(canonical_subset);
	
	if (f_v) {
		cout << "generator::map_to_canonical_k_subset done" << endl;
		}
}

void generator::get_representative_of_subset_orbit(
	INT *set, INT size, INT local_orbit_no, 
	strong_generators *&Strong_gens, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT fst, node, sz;
	oracle *O;

	if (f_v) {
		cout << "generator::get_representative_of_subset_orbit verbose_level=" << verbose_level << endl;
		}
	fst = first_oracle_node_at_level[size];
	node = fst + local_orbit_no;
	if (f_vv) {
		cout << "generator::get_representative_of_subset_orbit before get_set" << endl;
		}
	get_set(node, set, sz);
	if (sz != size) {
		cout << "get_representative_of_subset_orbit: sz != size" << endl;
		exit(1);
		}
	O = root + node;
	if (f_vv) {
		cout << "generator::get_representative_of_subset_orbit before get_stabilizer_generators" << endl;
		}
	O->get_stabilizer_generators(this, Strong_gens, 0);
	if (f_v) {
		cout << "generator::get_representative_of_subset_orbit done" << endl;
		}
}

void generator::print_fusion_nodes(INT depth)
{
	INT i, f, l, j, h;

	for (i = 0; i <= depth; i++) {
		f = first_oracle_node_at_level[i];
		l = nb_orbits_at_level(i);
		for (j = 0; j < l; j++) {
			oracle *O;

			O = &root[f + j];
			for (h = 0; h < O->nb_extensions; h++) {
				extension *E = O->E + h;

				if (E->type == EXTENSION_TYPE_FUSION) {
					cout << "fusion (" << f + j << "/" << h << ") -> (" << E->data1 << "/" << E->data2 << ")" << endl;
					}
				}
			}
		}
}

void generator::identify(INT *data, INT sz, INT *transporter, INT &orbit_at_level, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	INT f_implicit_fusion = FALSE;
	INT final_node;

	if (f_v) {
		cout << "generator::identify" << endl;
		}
	if (f_v) {
		cout << "generator::identify identifying the set ";
		INT_vec_print(cout, data, sz);
		cout << endl;
		}

	if (f_v) {
		cout << "generator::identify before recognize" << endl;
		}

	recognize(this, data, sz, transporter, f_implicit_fusion, 
		final_node, verbose_level);

	if (f_v) {
		cout << "generator::identify after recognize" << endl;
		}

	longinteger_object go;

	orbit_at_level = final_node - first_oracle_node_at_level[sz];
	get_stabilizer_order(sz, orbit_at_level, go);

	if (f_v) {
		cout << "generator::identify trace returns final_node = " << final_node << " which is isomorphism type " << orbit_at_level << " with ago=" << go << endl;
		}
	if (f_v) {
		cout << "generator::identify transporter:" << endl;
		A->element_print_quick(transporter, cout);
		}

	if (f_v) {
		cout << "generator::identify done" << endl;
		}

}

void generator::test_identify(INT level, INT nb_times, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	INT *transporter;
	INT f_implicit_fusion = FALSE;
	INT final_node;
	INT *Elt;
	INT nb_orbits, cnt, r, r2;
	INT *set1;
	INT *set2;
	sims *S;
	longinteger_object go;

	if (f_v) {
		cout << "generator::test_identify, level = " << level << " nb_times = " << nb_times << endl;
		}

	Elt = NEW_INT(A->elt_size_in_INT);
	transporter = NEW_INT(A->elt_size_in_INT);
	nb_orbits = nb_orbits_at_level(level);
	set1 = NEW_INT(level);
	set2 = NEW_INT(level);

	S = Strong_gens->create_sims(0 /*verbose_level*/);

	S->group_order(go);
	cout << "Group of order " << go << " created" << endl;




	for (cnt = 0; cnt < nb_times; cnt++) {
		r = random_integer(nb_orbits);
		if (f_v) {
			cout << "random orbit " << r << " / " << nb_orbits << endl;
			}
		get_set_by_level(level, r, set1);
		if (f_v) {
			cout << "random orbit " << r << " / " << nb_orbits << " is represented by ";
			INT_vec_print(cout, set1, level);
			cout << endl;
			}
		A->random_element(S, Elt, 0 /* verbose_level */);
		A2->map_a_set_and_reorder(set1, set2, level, Elt, 0 /* verbose_level */);
		cout << "mapped set is ";
		INT_vec_print(cout, set2, level);
		cout << endl;

		recognize(this, set2, level, transporter, f_implicit_fusion, 
			final_node, verbose_level);
		
		r2 = final_node - first_oracle_node_at_level[level];
		if (r2 != r) {
			cout << "recognition fails" << endl;
			exit(1);
			}
		else {
			cout << "recognition is successful" << endl;
			}
		}

	delete S;
	FREE_INT(Elt);
	FREE_INT(transporter);
	FREE_INT(set1);
	FREE_INT(set2);
	if (f_v) {
		cout << "generator::test_identify done" << endl;
		}
}


