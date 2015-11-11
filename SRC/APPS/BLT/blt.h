// blt.h
// 

#include "discreta.h"

typedef  int * pint;

typedef class blt_generator blt_generator;
typedef class orbit_stitching orbit_stitching;




// global data and global functions:

extern INT t0; // the system time when the program started



// ##################################################################################################
// blt_generator.C
// ##################################################################################################



class blt_generator {

public:
	const BYTE *override_poly;
	const BYTE *starter_directory_name;
	const BYTE *starter_prefix;
	
	generator *gen;
	action *A;
		
	orthogonal *O;
	INT f_orthogonal_allocated;
	
	INT f_BLT;
	INT f_ovoid;
	
	INT f_semilinear; // from the command line
	INT epsilon; // the type of the quadric (0, 1 or -1)
	INT n; // algebraic dimension
	INT q; // field order
	finite_field *F;
	
	INT target_size;

	INT nb_sol; // number of solutions so far

	INT f_override_schreier_depth;
	INT override_schreier_depth;

	INT f_override_n;
	INT override_n;
	
	INT f_override_epsilon;
	INT override_epsilon;

	
	void read_arguments(int argc, const char **argv);
	blt_generator();
	~blt_generator();
	void init_basic(INT q, 
		const BYTE *starter_directory_name, 
		const BYTE *starter_prefix, 
		int argc, const char **argv, 
		INT verbose_level);
	void init_group(INT verbose_level);
	void init_orthogonal(INT verbose_level);
	void init_orthogonal_hash(INT verbose_level);
	void init2(INT verbose_level);
	INT pair_test(INT a, INT x, INT y, INT verbose_level);
		// We assume that a is an element of a set S of size at least two such that 
		// S \cup \{ x \} is BLT and 
		// S \cup \{ y \} is BLT.
		// In order to test of S \cup \{ x, y \} is BLT, we only need to test 
		// the triple \{ x,y,a\}
	INT check_conditions(INT len, INT *S, INT verbose_level);
	INT collinearity_test(INT *line, INT len, INT verbose_level);
	void print(ostream &ost, INT *S, INT len);

	void create_system(INT level, 
		INT orbit_at_level, 
		INT level_of_candidates_file, 
		const BYTE *output_prefix, INT width, INT verbose_level);
	void create_graphs(INT starter_depth, 
		INT orbit_at_level_r, INT orbit_at_level_m, 
		INT level_of_candidates_file, 
		const BYTE *output_prefix, 
		INT f_lexorder_test, INT f_eliminate_graphs_if_possible, 
		INT verbose_level);
	INT create_graph(INT level, 
		INT orbit_at_level, INT level_of_candidates_file, 
		const BYTE *output_prefix, 
		INT f_lexorder_test, INT f_eliminate_graphs_if_possible, 
		INT &nb_vertices, BYTE *graph_fname_base, 
		colored_graph *&CG,  
		INT verbose_level);

	void compute_colors(INT orbit_at_level, 
		INT *starter, INT starter_sz, 
		INT special_line, 
		INT *candidates, INT nb_candidates, 
		INT *&point_color, INT &nb_colors, 
		INT verbose_level);
	void compute_adjacency_list_fast(INT first_point_of_starter, 
		INT *points, INT nb_points, INT *point_color, 
		UBYTE *&bitvector_adjacency, INT &bitvector_length_in_bits, INT &bitvector_length, 
		INT verbose_level);
	void write_problem_to_file_wassermann(
		const BYTE *output_prefix, 
		INT *Incma, 
		INT nb_free_pts, INT nb_candidates, 
		INT starter_level, INT starter_case, 
		INT width);
	void early_test_func(INT *S, INT len, 
		INT *candidates, INT nb_candidates, 
		INT *good_candidates, INT &nb_good_candidates, 
		INT verbose_level);
	INT check_function_incremental(INT len, INT *S, INT verbose_level);
	void Law_71(INT verbose_level);

	void report(isomorph &Iso, INT verbose_level);
	void subset_orbits(isomorph &Iso, INT verbose_level);
};

void print_set(ostream &ost, INT len, INT *S, void *data);
INT check_conditions(INT len, INT *S, void *data, INT verbose_level);
void early_test_func_callback(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	void *data, INT verbose_level);
INT check_function_incremental_callback(INT len, INT *S, void *data, INT verbose_level);
void callback_report(isomorph *Iso, void *data, INT verbose_level);
void callback_subset_orbits(isomorph *Iso, void *data, INT verbose_level);




// ##################################################################################################
// orbit_stitching.C
// ##################################################################################################


class orbit_stitching {
public:
	BYTE label[1000];
	
	action *A;
	matrix_group *M;
	orthogonal *O;
	finite_field *F;
	
	INT print_interval;
	
	vector_ge *gens;
	schreier *S;

	INT *Elt1, *Elt2, *Elt3, *Elt4, *Elt5, *Elt6;
	
	INT *good_orbits;
	INT nb_good_orbits;
	INT *orbit_length;
	INT *orbit_length_sorted;
	INT *orbit_length_sorting_perm;
	INT *orbit_length_sorting_perm_inv;
	INT orbit_length_nb_types;
	INT *orbit_length_type_first;
	INT *orbit_length_type_len;
	
	INT max_orbit_length;
	INT *short_orbits;
	INT nb_short_orbits;
	INT *long_orbits;
	INT nb_long_orbits;
	
	INT *orbit;
	
	INT *set;
	INT *test_set;
	
	INT adjacency_length; // = nb_good_orbits choose 2
	INT *adjacency; // [adjacency_length]

	INT starter_size; // size after preselecting orbits
	
	INT *search_set;
	INT search_set_size;
	
	INT search_depth;
	
	clique_finder *CF;
	INT search_steps;
	INT nb_sol;
	
	INT case_no;
	ofstream *fp_out;
	ofstream *fp_summary;

	orbit_stitching();
	~orbit_stitching();
	void null();
	void free();
	void init(action *A, vector_ge *gens, 
		INT print_interval, INT case_no,
		INT verbose_level);
	void compute_orbits(INT verbose_level);
	INT get_orbit(INT idx, INT *orbit);
	void compute_adjacency_list(INT verbose_level);
	void compute_adjacency_list_short_orbits(INT verbose_level);
	void print_adjacency_list_short_orbits();
	void preselect(INT nb_preselected, INT *preselected, INT verbose_level);
	void search(INT f_solution_file, const BYTE *solution_file_name, INT verbose_level);
	void clique_found(INT *current_clique, INT verbose_level);
	void add_point(INT pt, 
		INT current_clique_size, INT *current_clique, 
		INT verbose_level);
	void delete_point(INT pt, 
		INT current_clique_size, INT *current_clique, 
		INT verbose_level);
	INT is_viable(
		INT pt, 
		INT current_clique_size, INT *current_clique, 
		INT verbose_level);

};

void orbit_stitching_call_back_clique_found(clique_finder *CF, INT verbose_level);
void orbit_stitching_call_back_add_point(clique_finder *CF, 
	INT current_clique_size, INT *current_clique, 
	INT pt, INT verbose_level);
void orbit_stitching_call_back_delete_point(clique_finder *CF, 
	INT current_clique_size, INT *current_clique, 
	INT pt, INT verbose_level);
INT orbit_stitching_call_back_is_viable(clique_finder *CF, 
	INT current_clique_size, INT *current_clique, 
	INT pt, INT verbose_level);
	


