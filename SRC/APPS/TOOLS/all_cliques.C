// all_cliques.C
// 
// Anton Betten
// January 28, 2015
//
// 
//
//

#include "orbiter.h"
#include "discreta.h"


// global data:

INT t0; // the system time when the program started

void print_orbits_at_level(generator *gen, INT level, INT verbose_level);
void early_test_function_cliques(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	void *data, INT verbose_level);
void early_test_function_cocliques(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	void *data, INT verbose_level);

int main(int argc, char **argv)
{
	INT i, j;
	t0 = os_ticks();
	INT verbose_level = 0;
	INT f_file = FALSE;	
	const BYTE *fname = NULL;
	INT f_all_cliques = FALSE;
	INT f_all_cocliques = FALSE;
	INT f_draw_poset = FALSE;
	INT f_embedded = FALSE;
	INT nb_print_level = 0;
	INT print_level[1000];

	
	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-v") == 0) {
			verbose_level = atoi(argv[++i]);
			cout << "-v " << verbose_level << endl;
			}
		else if (strcmp(argv[i], "-file") == 0) {
			f_file = TRUE;
			fname = argv[++i];
			cout << "-file " << fname << endl;
			}
		else if (strcmp(argv[i], "-all_cliques") == 0) {
			f_all_cliques = TRUE;
			cout << "-all_cliques " << endl;
			}
		else if (strcmp(argv[i], "-all_cocliques") == 0) {
			f_all_cocliques = TRUE;
			cout << "-all_cocliques " << endl;
			}
		else if (strcmp(argv[i], "-draw_poset") == 0) {
			f_draw_poset = TRUE;
			cout << "-draw_poset " << endl;
			}
		else if (strcmp(argv[i], "-embedded") == 0) {
			f_embedded = TRUE;
			cout << "-embedded " << endl;
			}
		else if (strcmp(argv[i], "-print_level") == 0) {
			print_level[nb_print_level++] = atoi(argv[++i]);
			cout << "-embedded " << endl;
			}

		}

	if (!f_file) {
		cout << "Please specify the file name using -file <fname>" << endl;
		exit(1);
		}
	colored_graph *CG;

	CG = new colored_graph;

	CG->load(fname, verbose_level);




	INT *Adj;
	action *Aut;
	longinteger_object ago;

	cout << "computing automorphism group of the graph:" << endl;
	//Aut = create_automorphism_group_of_colored_graph_object(CG, verbose_level);


	Adj = NEW_INT(CG->nb_points * CG->nb_points);
	INT_vec_zero(Adj, CG->nb_points * CG->nb_points);
	for (i = 0; i < CG->nb_points; i++) {
		for (j = i + 1; j < CG->nb_points; j++) {
			if (CG->is_adjacent(i, j)) {
				Adj[i * CG->nb_points + j] = 1;
				}
			}
		}
	Aut = create_automorphism_group_of_graph(Adj, CG->nb_points, verbose_level);

	Aut->group_order(ago);	
	cout << "ago=" << ago << endl;

	action *Aut_on_points;
	INT *points;
	
	Aut_on_points = new action;
	points = NEW_INT(CG->nb_points);
	for (i = 0; i < CG->nb_points; i++) {
		points[i] = i;
		}

	Aut_on_points->induced_action_by_restriction(*Aut, 
		TRUE /* f_induce_action */, Aut->Sims, 
		CG->nb_points /* nb_points */, points, verbose_level);
	
	Aut_on_points->group_order(ago);	
	cout << "ago on points = " << ago << endl;

	{
	schreier S;
	strong_generators SG;
	
	Aut_on_points->compute_strong_generators_from_sims(verbose_level);
	SG.init_from_sims(Aut_on_points->Sims, verbose_level);
	Aut_on_points->compute_all_point_orbits(S, 
		*SG.gens, verbose_level);

		/*all_point_orbits(S, verbose_level);*/
	cout << "has " << S.nb_orbits << " orbits on points" << endl;
	}
	
	BYTE prefix[1000];
	generator *gen;
	INT nb_orbits, depth;
	


	if (f_all_cliques) {

		strcpy(prefix, fname);
		replace_extension_with(prefix, "_cliques");


		compute_orbits_on_subsets(gen, 
			CG->nb_points /* target_depth */,
			prefix, 
			FALSE /* f_W */, FALSE /* f_w */,
			Aut_on_points, Aut_on_points, 
			Aut_on_points->Strong_gens, 
			early_test_function_cliques,
			CG, 
			NULL, 
			NULL, 
			verbose_level);
		}
	else {

		strcpy(prefix, fname);
		replace_extension_with(prefix, "_cocliques");

		compute_orbits_on_subsets(gen, 
			CG->nb_points /* target_depth */,
			prefix, 
			FALSE /* f_W */, FALSE /* f_w */,
			Aut_on_points, Aut_on_points, 
			Aut_on_points->Strong_gens, 
			early_test_function_cocliques,
			CG, 
			NULL, 
			NULL, 
			verbose_level);
		}

	for (depth = 0; depth < CG->nb_points; depth++) {
		nb_orbits = gen->nb_orbits_at_level(depth);
		if (nb_orbits == 0) {
			depth--;
			break;
			}
		}

	if (f_all_cliques) {
		cout << "the largest cliques have size " << depth << endl;
		for (i = 0; i <= depth; i++) {
			nb_orbits = gen->nb_orbits_at_level(i);
			cout << setw(3) << i << " : " << setw(3) << nb_orbits << endl;
			}
		}
	else if (f_all_cocliques) {
		cout << "the largest cocliques have size " << depth << endl;
		for (i = 0; i <= depth; i++) {
			nb_orbits = gen->nb_orbits_at_level(i);
			cout << setw(3) << i << " : " << setw(3) << nb_orbits << endl;
			}
		}

	if (f_draw_poset) {
		gen->draw_poset(gen->fname_base, depth, 0 /* data1 */, f_embedded, verbose_level);
		}

	print_orbits_at_level(gen, depth, verbose_level);
	if (nb_print_level) {
		for (i = 0; i < nb_print_level; i++) {
			print_orbits_at_level(gen, print_level[i], verbose_level);
			}
		}

	FREE_INT(Adj);
	FREE_INT(points);
	delete Aut_on_points;
	delete Aut;
	delete CG;

	the_end(t0);
	//the_end_quietly(t0);

}

void print_orbits_at_level(generator *gen, INT level, INT verbose_level)
{
	INT *set;
	longinteger_object go, ol, ago;
	longinteger_domain D;
	INT i, nb_orbits;

	set = NEW_INT(level);
	nb_orbits = gen->nb_orbits_at_level(level);


	gen->A->group_order(ago);
	cout << "group order " << ago << endl;
	cout << "The " << nb_orbits << " orbits at level " << level << " are:" << endl;
	cout << "orbit : representative : stabilizer order : orbit length" << endl;
	for (i = 0; i < nb_orbits; i++) {
		
		gen->get_set_by_level(level, i, set);

		strong_generators *gens;
		gen->get_stabilizer_generators(gens,  
			level, i, 0 /*verbose_level*/);
		gens->group_order(go);
		D.integral_division_exact(ago, go, ol);

		
		cout << "Orbit " << i << " is the set ";
		INT_vec_print(cout, set, level);
		cout << " : " << go << " : " << ol << endl;
		//cout << endl;

		
		}

	FREE_INT(set);
}

void early_test_function_cliques(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	void *data, INT verbose_level)
{
	colored_graph *CG = (colored_graph *) data;
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "early_test_function for set ";
		print_set(cout, len, S);
		cout << endl;
		}

	CG->early_test_func_for_clique_search(S, len, 
		candidates, nb_candidates, 
		good_candidates, nb_good_candidates, 
		verbose_level - 2);


	if (f_v) {
		cout << "early_test_function done" << endl;
		}
}

void early_test_function_cocliques(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	void *data, INT verbose_level)
{
	colored_graph *CG = (colored_graph *) data;
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "early_test_function for set ";
		print_set(cout, len, S);
		cout << endl;
		}

	CG->early_test_func_for_coclique_search(S, len, 
		candidates, nb_candidates, 
		good_candidates, nb_good_candidates, 
		verbose_level - 2);


	if (f_v) {
		cout << "early_test_function done" << endl;
		}
}


