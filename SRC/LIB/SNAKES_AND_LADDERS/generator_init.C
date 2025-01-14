// generator_init.C
//
// Anton Betten
// December 29, 2003
//
// moved here from generator.C: July 29, 2014


#include "orbiter.h"

generator::generator()
{
	null();
}

generator::~generator()
{
	freeself();
}

void generator::null()
{
	problem_label[0] = 0;
	
	f_candidate_check_func = FALSE;
	f_candidate_incremental_check_func = FALSE;
	f_print_function = FALSE;
	Elt_memory = NULL;
	A = NULL;
	Strong_gens = NULL;
	//SG0 = NULL;
	//transversal_length = NULL;
	S = NULL;

	tmp_set_apply_fusion = NULL;
	tmp_find_node_for_subspace_by_rank1 = NULL;
	tmp_find_node_for_subspace_by_rank2 = NULL;
	tmp_find_node_for_subspace_by_rank3 = NULL;


	nb_times_trace = 0;
	nb_times_trace_was_saved = 0;
	
	sz = 0;
	transporter = NULL;
	set = NULL;
	
	
	nb_oracle_nodes_used = 0;
	root = NULL;
	first_oracle_node_at_level = NULL;
	set0 = NULL;
	set1 = NULL;
	set3 = NULL;
	nb_extension_nodes_at_level_total = NULL;
	nb_extension_nodes_at_level = NULL;
	nb_fusion_nodes_at_level = NULL;
	nb_unprocessed_nodes_at_level = NULL;


	f_prefix = FALSE;
	prefix[0] = 0;
	
	f_max_depth = FALSE;
	
	f_extend = FALSE;
	f_recover = FALSE;
		
	f_w = FALSE;
	f_W = FALSE;
	f_t = FALSE;
	f_T = FALSE;
	f_log = FALSE;
	f_Log = FALSE;
	f_print_only = FALSE;
	f_find_group_order = FALSE;
	
	f_has_invariant_subset_for_root_node = FALSE;
	invariant_subset_for_root_node = NULL;
	
	verbose_level = 0;
	verbose_level_group_theory = 0;
	
	fname_base[0] = 0;
	
	xmax = 1000000;
	ymax = 1000000;
	radius = 300;
	
	f_starter = FALSE;
	f_downstep_split = FALSE;
	f_upstep_split = FALSE;
	f_downstep_collate = FALSE;
	f_upstep_collate = FALSE;
	split_mod = 0;
	split_case = 0;
	
	f_do_group_extension_in_upstep = TRUE;

	f_allowed_to_show_group_elements = FALSE;	
	downstep_orbits_print_max_orbits = 25;
	downstep_orbits_print_max_points_per_orbit = 50;

	f_on_subspaces = FALSE;
	rank_point_func = NULL;
	unrank_point_func = NULL;

	f_early_test_func = FALSE;
	early_test_func = NULL;
	f_its_OK_to_not_have_an_early_test_func = FALSE;

	depth = 0;

	CFI = NULL;
	
	t0 = os_ticks();
}

void generator::freeself()
{
	INT i;
	INT f_v = FALSE;
	
	if (f_v) {
		cout << "generator::freeself" << endl;
		}
	if (Elt_memory) {
		FREE_INT(Elt_memory);
		}

	// do not free Strong_gens
	

	if (f_v) {
		cout << "generator::freeself deleting S" << endl;
		}
	if (S) {
		FREE_INT(S);
		}
	if (tmp_set_apply_fusion) {
		FREE_INT(tmp_set_apply_fusion);
		}
	if (tmp_find_node_for_subspace_by_rank1) {
		FREE_INT(tmp_find_node_for_subspace_by_rank1);
		}
	if (tmp_find_node_for_subspace_by_rank2) {
		FREE_INT(tmp_find_node_for_subspace_by_rank2);
		}
	if (tmp_find_node_for_subspace_by_rank3) {
		FREE_INT(tmp_find_node_for_subspace_by_rank3);
		}

	if (f_v) {
		cout << "generator::freeself deleting transporter and set[]" << endl;
		}
	if (transporter) {
		delete transporter;
		for (i = 0; i <= sz; i++) {
			FREE_INT(set[i]);			
			}
		FREE_PINT(set);
		}
	if (f_v) {
		cout << "generator::freeself  before exit_oracle" << endl;
		}
	exit_oracle();
	if (f_v) {
		cout << "generator::freeself  after exit_oracle" << endl;
		}


	if (f_v) {
		cout << "generator::freeself done" << endl;
		}
	null();
}

void generator::usage()
{
	cout << "generator options:" << endl;
	cout << "-v <n>" << endl;
	cout << "  verbose level n" << endl;
	cout << "-gv <n>" << endl;
	cout << "  verbose level n for all group theory related things" << endl;
	cout << "-w" << endl;
	cout << "  write output in level files (only last level)" << endl;
	cout << "-W" << endl;
	cout << "  write output in level files (all levels)" << endl;
	cout << "-depth <n>" << endl;
	cout << "  compute up to depth n" << endl;
	cout << "-prefix <s>" << endl;
	cout << "  use s as prefix for all output files" << endl;
	cout << "-extend <f> <t> <r> <m> <s>" << endl;
	cout << "  extend all partial solutions congruent" << endl;
	cout << "  to <r> mod <m> from file <s> from level <f> to level <t>" << endl;
	cout << "-x xmax" << endl;
	cout << "   specifies horizontal size (default 1000)" << endl;
	cout << "-y ymax" << endl;
	cout << "   specifies vertical size (default 1000)" << endl;
	cout << "-rad r" << endl;
	cout << "   specifies radius (default 300)" << endl;
	cout << "-t" << endl;
	cout << "   draw tree at last level only" << endl;
	cout << "-T" << endl;
	cout << "   draw tree each level" << endl;
	cout << "-log" << endl;
	cout << "   log nodes at the end only" << endl;
	cout << "-Log" << endl;
	cout << "   log nodes at every step" << endl;
	cout << "-r <fname>" << endl;
	cout << "   recover from data file <fname>" << endl;
	cout << "-printonly <fname>" << endl;
	cout << "   print only (no computation), to be used with -r" << endl;
	cout << "-findgroup <order>" << endl;
	cout << "   find group of order <order>" << endl;
	cout << "-downstep_split <mod> <case>" << endl;
	cout << "   Process downstep for all cases congruent <case> modulo <mod>" << endl;
	cout << "-upstep_split <mod> <case>" << endl;
	cout << "   Process upstep for all cases congruent <case> modulo <mod>" << endl;
	cout << "-downstep_collate <mod>" << endl;
	cout << "   Collates the data from the parallel downstep runs modulo <mod>" << endl;
	cout << "-upstep_collate <mod>" << endl;
	cout << "   Collates the data from the parallel downstep runs modulo <mod>" << endl;
}

void generator::read_arguments(int argc, const char **argv, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;
	
	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-v") == 0) {
			i++;
			generator::verbose_level = atoi(argv[i]);
			if (f_v) {
				cout << "-v " << generator::verbose_level << endl;
				}
			}
		else if (strcmp(argv[i], "-gv") == 0) {
			i++;
			verbose_level_group_theory = atoi(argv[i]);
			if (f_v) {
				cout << "-gv " << verbose_level_group_theory << endl;
				}
			}
		else if (strcmp(argv[i], "-prefix") == 0) {
			i++;
			f_prefix = TRUE;
			strcpy(prefix, argv[i]);
			if (f_v) {
				cout << "-prefix " << prefix << endl;
				}
			}
		else if (strcmp(argv[i], "-w") == 0) {
			f_w = TRUE;
			if (f_v) {
				cout << "-w" << endl;
				}
			}
		else if (strcmp(argv[i], "-W") == 0) {
			f_W = TRUE;
			if (f_v) {
				cout << "-W" << endl;
				}
			}
		else if (strcmp(argv[i], "-t") == 0) {
			f_t = TRUE;
			if (f_v) {
				cout << "-t" << endl;
				}
			}
		else if (strcmp(argv[i], "-T") == 0) {
			f_T = TRUE;
			if (f_v) {
				cout << "-T" << endl;
				}
			}
		else if (strcmp(argv[i], "-log") == 0) {
			f_log = TRUE;
			if (f_v) {
				cout << "-log" << endl;
				}
			}
		else if (strcmp(argv[i], "-Log") == 0) {
			f_Log = TRUE;
			if (f_v) {
				cout << "-Log" << endl;
				}
			}
		else if (strcmp(argv[i], "-x") == 0) {
			xmax = atoi(argv[i + 1]);
			i++;
			if (f_v) {
				cout << "-x " << xmax << endl;
				}
			}
		else if (strcmp(argv[i], "-y") == 0) {
			ymax = atoi(argv[i + 1]);
			i++;
			if (f_v) {
				cout << "-y " << ymax << endl;
				}
			}
		else if (strcmp(argv[i], "-rad") == 0) {
			radius = atoi(argv[i + 1]);
			i++;
			if (f_v) {
				cout << "-rad " << radius << endl;
				}
			}
		else if (strcmp(argv[i], "-depth") == 0) {
			f_max_depth = TRUE;
			max_depth = atoi(argv[++i]);
			if (f_v) {
				cout << "-depth " << max_depth << endl;
				}
			}
		else if (strcmp(argv[i], "-extend") == 0) {
			f_extend = TRUE;
			extend_from = atoi(argv[++i]);
			extend_to = atoi(argv[++i]);
			extend_r = atoi(argv[++i]);
			extend_m = atoi(argv[++i]);
			strcpy(extend_fname, argv[++i]);
			if (f_v) {
				cout << "-extend from level " << extend_from 
					<< " to level " << extend_to 
					<< " cases congruent " << extend_r
					<< " mod " << extend_m
					<< " from file " << extend_fname << endl;
				}
			}
		else if (strcmp(argv[i], "-recover") == 0) {
			f_recover = TRUE;
			recover_fname = argv[++i];
			if (f_v) {
				cout << "-recover " << recover_fname << endl; 
				}
			}
		else if (strcmp(argv[i], "-printonly") == 0) {
			f_print_only = TRUE;
			if (f_v) {
				cout << "-printonly" << endl; 
				}
			}
		else if (strcmp(argv[i], "-findgroup") == 0) {
			f_find_group_order = TRUE;
			find_group_order = atoi(argv[++i]);
			if (f_v) {
				cout << "-findgroup " << find_group_order << endl;
				}
			}
		else if (strcmp(argv[i], "-downstep_split") == 0) {
			f_downstep_split = TRUE;
			split_mod = atoi(argv[++i]);
			split_case = atoi(argv[++i]);
			if (f_v) {
				cout << "-downstep_split " << split_mod << " " << split_case << endl; 
				}
			}
		else if (strcmp(argv[i], "-upstep_split") == 0) {
			f_upstep_split = TRUE;
			split_mod = atoi(argv[++i]);
			split_case = atoi(argv[++i]);
			if (f_v) {
				cout << "-upstep_split " << split_mod << " " << split_case << endl; 
				}
			}
		else if (strcmp(argv[i], "-downstep_collate") == 0) {
			f_downstep_collate = TRUE;
			split_mod = atoi(argv[++i]);
			split_case = 0;
			if (f_v) {
				cout << "-downstep_collate " << split_mod << endl; 
				}
			}
		else if (strcmp(argv[i], "-upstep_collate") == 0) {
			f_upstep_collate = TRUE;
			split_mod = atoi(argv[++i]);
			split_case = 0;
			if (f_v) {
				cout << "-upstep_collate " << split_mod << endl; 
				}
			}
		}
}

void generator::init(action *A, action *A2, 
	strong_generators *gens, 
	INT sz, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	//INT f_v6 = (verbose_level >= 6);
	INT i;
	
	if (f_v) {
		cout << "generator::init" << endl;
		cout << "generator::init sz = " << sz << endl;
		cout << "generator::init A->degree=" << A->degree << endl;
		cout << "generator::init A2->degree=" << A2->degree << endl;
		cout << "generator::init sz = " << sz << endl;
		}
	t0 = os_ticks();

	progress_epsilon = 0.005;

	
	if (f_vv) {
		cout << "generator::init action A:" << endl;
		A->print_info();
		cout << "generator::init action A2:" << endl;
		A2->print_info();
		}
	gens->group_order(go);

	generator::A = A;
	generator::A2 = A2;

	if (f_v) {
		cout << "generator::init group order is ";
		cout << go << endl;
		}
	
	generator::sz = sz;
	
	if (f_v) {
		cout << "generator::init sz = " << sz << endl;
		}
	
	Strong_gens = gens;



	if (f_vv) {
		cout << "generator::init allocating S of size " << sz << endl;
		}
	S = NEW_INT(sz);
	for (i = 0; i < sz; i++) {
		S[i] = i;
		}

	tmp_set_apply_fusion = NEW_INT(sz);

	if (f_vv) {
		cout << "generator::init allocating Elt_memory" << endl;
		}


	Elt_memory = NEW_INT(5 * A->elt_size_in_INT);
	Elt1 = Elt_memory + 0 * A->elt_size_in_INT;
	Elt2 = Elt_memory + 1 * A->elt_size_in_INT;
	Elt3 = Elt_memory + 2 * A->elt_size_in_INT;
	Elt4 = Elt_memory + 3 * A->elt_size_in_INT;
	Elt5 = Elt_memory + 4 * A->elt_size_in_INT;
	
	transporter = new vector_ge;
	transporter->init(A);
	transporter->allocate(sz + 1);
	A->element_one(transporter->ith(0), FALSE);
	
	set = NEW_PINT(sz + 1);
	for (i = 0; i <= sz; i++) {
		set[i] = NEW_INT(sz);
		}
		
	nb_oracle_nodes_used = 0;
	nb_oracle_nodes_allocated = 0;

	nb_times_image_of_called0 = A->nb_times_image_of_called;
	nb_times_mult_called0 = A->nb_times_mult_called;
	nb_times_invert_called0 = A->nb_times_invert_called;
	nb_times_retrieve_called0 = A->nb_times_retrieve_called;
	nb_times_store_called0 = A->nb_times_store_called;

	if (f_v) {
		cout << "generator::init done" << endl;
		}
}



void generator::initialize_with_starter(action *A_base, action *A_use, 
	strong_generators *gens, 
	INT depth, 
	BYTE *prefix, 
	INT starter_size, 
	INT *starter, 
	strong_generators *Starter_Strong_gens, 
	INT *starter_live_points, 
	INT starter_nb_live_points, 
	void *starter_canonize_data, 
	INT (*starter_canonize)(INT *Set, INT len, INT *Elt, void *data, INT verbose_level), 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);

	if (f_v) {
		cout << "generator::initialize_with_starter" << endl;
		}

	generator::f_prefix = TRUE;
	strcpy(generator::fname_base, prefix);

	generator::depth = depth;
	downstep_orbits_print_max_orbits = 50;
	downstep_orbits_print_max_points_per_orbit = INT_MAX;
	

	// !!!
	//f_allowed_to_show_group_elements = TRUE;

	if (f_vv) {
		cout << "generator::initialize_with_starter calling gen->init" << endl;
		}
	init(A_base, A_use, 
		gens, 
		//*gens, tl, 
		depth, verbose_level - 2);
	

	if (f_vv) {
		cout << "generator::initialize_with_starter calling init_starter" << endl;
		}
	init_starter(starter_size, 
		starter, 
		Starter_Strong_gens, 
		starter_live_points, 
		starter_nb_live_points, 
		starter_canonize_data, 
		starter_canonize, 
		verbose_level - 2);

	INT nb_oracle_nodes = 1000;
	
	if (f_vv) {
		cout << "generator::initialize_with_starter calling gen->init_oracle" << endl;
		}
	init_oracle(nb_oracle_nodes, verbose_level - 1);
	if (f_vv) {
		cout << "generator::initialize_with_starter calling gen->init_root_node" << endl;
		}
	init_root_node(verbose_level);

	if (f_v) {
		cout << "generator::initialize_with_starter done" << endl;
		}
}

void generator::initialize(action *A_base, action *A_use, 
	strong_generators *gens, 
	INT depth, 
	BYTE *prefix, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);

	if (f_v) {
		cout << "generator::initialize" << endl;
		cout << "generator::initialize depth = " << depth << endl;
		}

	generator::f_prefix = TRUE;
	strcpy(generator::fname_base, prefix);

	generator::depth = depth;
	downstep_orbits_print_max_orbits = 50;
	downstep_orbits_print_max_points_per_orbit = INT_MAX;
	

	// !!!
	//f_allowed_to_show_group_elements = TRUE;

	if (f_vv) {
		cout << "generator::initialize calling gen->init" << endl;
		}
	init(A_base, A_use, 
		gens, 
		depth, verbose_level - 2);
	
	INT nb_oracle_nodes = 1000;
	
	if (f_vv) {
		cout << "generator::initialize calling gen->init_oracle" << endl;
		}
	init_oracle(nb_oracle_nodes, verbose_level - 1);
	if (f_vv) {
		cout << "generator::initialize calling gen->init_root_node" << endl;
		}
	init_root_node(verbose_level - 1);

	if (f_v) {
		cout << "generator::initialize done" << endl;
		}
}

void generator::init_root_node_invariant_subset(
	INT *invariant_subset, INT invariant_subset_size, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "generator::init_root_node_invariant_subset" << endl;
		}
	f_has_invariant_subset_for_root_node = TRUE;
	invariant_subset_for_root_node = invariant_subset;
	invariant_subset_for_root_node_size = invariant_subset_size;
	if (f_v) {
		cout << "generator::init_root_node_invariant_subset installed invariant subset of size " << invariant_subset_size << endl;
		}
}


void generator::init_root_node(INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	
	if (f_v) {
		cout << "generator::init_root_node" << endl;
		}
	if (f_starter) {
		INT i;
		
		first_oracle_node_at_level[0] = 0;
		root[0].freeself();
		root[0].node = 0;
		root[0].prev = -1;
		root[0].nb_strong_generators = 0;
		root[0].sv = NULL;
		for (i = 0; i < starter_size; i++) {
			
			nb_extension_nodes_at_level_total[i] = 0;
			nb_extension_nodes_at_level[i] = 0;
			nb_fusion_nodes_at_level[i] = 0;
			nb_unprocessed_nodes_at_level[i] = 0;
			
			if (f_vv) {
				cout << "generator::init_root_node initializing node at level " << i << endl;
				}
			first_oracle_node_at_level[i + 1] = first_oracle_node_at_level[i] + 1;
			root[i].E = new extension[1];
			root[i].nb_extensions = 1;
			root[i].E[0].type = EXTENSION_TYPE_EXTENSION;
			root[i].E[0].data = i + 1;
			root[i + 1].freeself();
			root[i + 1].node = i + 1;
			root[i + 1].prev = i;
			root[i + 1].pt = starter[i];
			root[i + 1].nb_strong_generators = 0;
			root[i + 1].sv = NULL;
			}
		if (f_vv) {
			cout << "generator::init_root_node storing strong generators" << endl;
			}
		root[starter_size].store_strong_generators(this, starter_strong_gens);
		first_oracle_node_at_level[starter_size + 1] = starter_size + 1;
		if (f_vv) {
			cout << "i : first_oracle_node_at_level[i]" << endl;
			for (i = 0; i <= starter_size + 1; i++) {
				cout << i << " : " << first_oracle_node_at_level[i] << endl;
				}
			}
		}
	else {
		root[0].init_root_node(this, verbose_level - 1);
		}
	if (f_v) {
		cout << "generator::init_root_node done" << endl;
		}
}

void generator::init_oracle(INT nb_oracle_nodes, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;

	if (f_v) {
		cout << "generator::init_oracle" << endl;
		}
	root = new oracle[nb_oracle_nodes];
	for (i = 0; i < nb_oracle_nodes; i++) {
		root[i].node = i;
		}
	nb_oracle_nodes_allocated = nb_oracle_nodes;
	nb_oracle_nodes_used = 0;
	oracle_nodes_increment = nb_oracle_nodes;
	oracle_nodes_increment_last = nb_oracle_nodes;
	first_oracle_node_at_level = NEW_INT(sz + 2);
	first_oracle_node_at_level[0] = 0;
	first_oracle_node_at_level[1] = 1;
	set0 = NEW_INT(sz + 1);
	set1 = NEW_INT(sz + 1);
	set3 = NEW_INT(sz + 1);
	nb_extension_nodes_at_level_total = NEW_INT(sz + 1);
	nb_extension_nodes_at_level = NEW_INT(sz + 1);
	nb_fusion_nodes_at_level = NEW_INT(sz + 1);
	nb_unprocessed_nodes_at_level = NEW_INT(sz + 1);
	for (i = 0; i < sz + 1; i++) {
		nb_extension_nodes_at_level_total[i] = 0;
		nb_extension_nodes_at_level[i] = 0;
		nb_fusion_nodes_at_level[i] = 0;
		nb_unprocessed_nodes_at_level[i] = 0;
		}
	if (f_v) {
		cout << "generator::init_oracle done" << endl;
		}
}


void generator::exit_oracle()
{
	if (root) {
		delete [] root;
		root = NULL;
		}
	if (set0) {
		FREE_INT(set0);
		set0 = NULL;
		}
	if (set1) {
		FREE_INT(set1);
		set1 = NULL;
		}
	if (set3) {
		FREE_INT(set3);
		set3 = NULL;
		}
	if (first_oracle_node_at_level) {
		FREE_INT(first_oracle_node_at_level);
		first_oracle_node_at_level = NULL;
		}

	if (nb_extension_nodes_at_level_total) {
		FREE_INT(nb_extension_nodes_at_level_total);
		nb_extension_nodes_at_level_total = NULL;
		}
	if (nb_extension_nodes_at_level) {
		FREE_INT(nb_extension_nodes_at_level);
		nb_extension_nodes_at_level = NULL;
		}
	if (nb_fusion_nodes_at_level) {
		FREE_INT(nb_fusion_nodes_at_level);
		nb_fusion_nodes_at_level = NULL;
		}
	if (nb_unprocessed_nodes_at_level) {
		FREE_INT(nb_unprocessed_nodes_at_level);
		nb_unprocessed_nodes_at_level = NULL;
		}
}

void generator::reallocate()
{
	INT increment_new;
	INT verbose_level = 0;
	
	increment_new = oracle_nodes_increment + oracle_nodes_increment_last;
	reallocate_to(nb_oracle_nodes_allocated + oracle_nodes_increment, verbose_level - 1);
	oracle_nodes_increment_last = oracle_nodes_increment;
	oracle_nodes_increment = increment_new;
	
}

void generator::reallocate_to(INT new_number_of_nodes, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;
	oracle *new_root;
	
	if (f_v) {
		cout << "generator::reallocate_to" << endl;
		}
	if (new_number_of_nodes < nb_oracle_nodes_allocated) {
		cout << "generator::reallocate_to new_number_of_nodes < nb_oracle_nodes_allocated" << endl;
		exit(1);
		}
	if (f_v) {
		cout << "generator::reallocate_to from " << nb_oracle_nodes_allocated << " to " << new_number_of_nodes << endl;
		}
	new_root = new oracle[new_number_of_nodes];
	for (i = 0; i < nb_oracle_nodes_allocated; i++) {
		new_root[i] = root[i];
		root[i].null();
		}
	delete [] root;
	root = new_root;
	nb_oracle_nodes_allocated = new_number_of_nodes;
	if (f_v) {
		cout << "generator::reallocate_to done" << endl;
		}
}


void generator::init_check_func(
	INT (*candidate_check_func)(INT len, INT *S, void *data, INT verbose_level), 
	void *candidate_check_data)
{
	f_candidate_check_func = TRUE;
	generator::candidate_check_func = candidate_check_func;
	generator::candidate_check_data = candidate_check_data;
}

void generator::init_incremental_check_func(
	INT (*candidate_incremental_check_func)(INT len, INT *S, void *data, INT verbose_level), 
	void *candidate_incremental_check_data)
{
	f_candidate_incremental_check_func = TRUE;
	generator::candidate_incremental_check_func = candidate_incremental_check_func;
	generator::candidate_incremental_check_data = candidate_incremental_check_data;
}

void generator::init_starter(INT starter_size, 
	INT *starter, 
	strong_generators *starter_strong_gens, 
	INT *starter_live_points, 
	INT starter_nb_live_points, 
	void *starter_canonize_data, 
	INT (*starter_canonize)(INT *Set, INT len, INT *Elt, void *data, INT verbose_level), 
	INT verbose_level)
// Does not initialize the first starter nodes. This is done in init_root_node 
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "generator::init_starter starter: ";
		INT_vec_print(cout, starter, starter_size);
		cout << " with " << starter_strong_gens->gens->len << " strong generators and " 
			<< starter_nb_live_points << " live points" << endl;
		}
	f_starter = TRUE;
	generator::starter_size = starter_size;
	generator::starter = starter;
	generator::starter_strong_gens = starter_strong_gens; 
	generator::starter_live_points = starter_live_points;
	generator::starter_nb_live_points = starter_nb_live_points;
	generator::starter_canonize_data = starter_canonize_data;
	generator::starter_canonize = starter_canonize;
	starter_canonize_Elt = NEW_INT(A->elt_size_in_INT);
}

void generator::init_vector_space_action(INT vector_space_dimension, 
	finite_field *F, 
	INT (*rank_point_func)(INT *v, void *data), 
	void (*unrank_point_func)(INT *v, INT rk, void *data),
	void *data,  
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "generator::init_vector_space_action vector_space_dimension=" << vector_space_dimension << endl;
		}
	f_on_subspaces = TRUE;
	generator::vector_space_dimension = vector_space_dimension;
	generator::F = F;
	generator::rank_point_func = rank_point_func;
	generator::unrank_point_func = unrank_point_func;
	generator::rank_point_data = data;

	tmp_find_node_for_subspace_by_rank1 = NEW_INT(vector_space_dimension);
	tmp_find_node_for_subspace_by_rank2 = NEW_INT(sz * vector_space_dimension);
	tmp_find_node_for_subspace_by_rank3 = NEW_INT(vector_space_dimension);
}

void generator::init_early_test_func(
	void (*early_test_func)(INT *S, INT len, 
		INT *candidates, INT nb_candidates, 
		INT *good_candidates, INT &nb_good_candidates, 
		void *data, INT verbose_level), 
	void *data,  
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "generator::init_early_test_func" << endl;
		}
	f_early_test_func = TRUE;
	generator::early_test_func = early_test_func;
	early_test_func_data = data;
}


