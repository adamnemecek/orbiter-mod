// hyperoval.h
//
// Anton Betten
// December 6, 2004


typedef class arc_generator arc_generator;

extern INT t0; // the system time when the program started

void print_usage();


#define MAX_FILES 1000


class arc_generator {

public:

	INT q;
	finite_field *F;

	INT nb_points_total;
	INT target_size;
	INT block_size;

	BYTE prefix_starter[1000]; // directory prefix
	BYTE base_fname[1000]; // no directory prefix, created in init()

	const BYTE *prefix_lift;
	INT f_starter_directory_name;
	const BYTE *starter_directory_name;
	INT f_solution_directory_name;
	const BYTE *solution_directory_name;

	INT f_compute_starter;
	INT f_starter_size;
	INT starter_size;
	INT f_build_db;
	INT f_draw_poset;
	INT f_lift;
	INT f_read_solutions;
	INT f_compute_orbits;
	INT f_classify;
	INT f_decomposition_matrix;
	INT f_report;

	INT f_recognize;
	INT *recognize_set;
	INT recognize_set_sz;

	const BYTE *solution_file[MAX_FILES];
	INT nb_solution_files;

	INT f_draw_system;
	const BYTE *fname_system;
	INT f_write_tree;
	const BYTE *fname_tree;


	INT f_lex;
	INT f_split;
	INT split_r;
	INT split_m;

	INT f_solve;
	INT f_save;
	INT f_read_instead;

	INT f_no_arc_testing;


	action *A;
	
	grassmann *Grass;
	action_on_grassmannian *AG;
	action *A_on_lines;



	
	//matrix_group *M;
	
	projective_space *P2;
	
	rank_checker rc;
		
	INT verbose_level;
	INT f_write_graph;
	INT f_maxdepth;
	INT maxdepth;
	


	void main(int argc, const char **argv);

	arc_generator();
	~arc_generator();
	void null();
	void freeself();
	void prepare_prefix_iso(BYTE *prefix_iso);
	void prepare_prefix_starter();
	void prepare_fname_base(BYTE *fname_base);
	generator *prepare_generator(int argc, const char **argv, INT verbose_level);
	void compute_starter(int argc, const char **argv, 
		INT f_recognize, INT *set, INT sz, 
		INT verbose_level);
	void build_db(int argc, const char **argv, INT verbose_level);
	void init(int argc, const char **argv);
	void read_arguments(int argc, const char **argv);


	void early_test_func(INT *S, INT len, 
		INT *candidates, INT nb_candidates, 
		INT *good_candidates, INT &nb_good_candidates, 
		INT verbose_level);
	INT check_arc(INT *S, INT len, INT verbose_level);
	void print_set_in_affine_plane(ostream &ost, INT len, INT *S);
	void point_unrank(INT *v, INT rk);
	INT point_rank(INT *v);
	void lifting_prepare_function_new(exact_cover *E, INT starter_case, 
		INT *candidates, INT nb_candidates, strong_generators *Strong_gens, 
		diophant *&Dio, INT *&col_labels, 
		INT &f_ruled_out, 
		INT verbose_level);
	// compute the incidence matrix of tangent lines versus candidate points
	// extended by external lines versus candidate points
#if 0
	void lifting_prepare_function(exact_cover *E, INT starter_case, 
		INT *candidates, INT nb_candidates, strong_generators *Strong_gens, 
		INT *&Inc, INT &nb_free_points, INT &nb_live_blocks2, INT *&live_blocks2, INT &nb_needed, INT &f_ruled_out, 
		INT verbose_level);
	// compute the incidence matrix of tangent lines versus candidate points
#endif
	INT arc_test(INT *S, INT len, INT verbose_level);
	void report(isomorph &Iso, INT verbose_level);
	void report_decompositions(isomorph &Iso, ofstream &f, INT orbit, 
		INT *data, INT verbose_level);
	void report_stabilizer(isomorph &Iso, ofstream &f, INT orbit, INT verbose_level);
};


INT callback_arc_test(exact_cover *EC, INT *S, INT len, void *data, INT verbose_level);
INT check_arc(INT len, INT *S, void *data, INT verbose_level);
INT placebo_test_function(INT len, INT *S, void *data, INT verbose_level);
void arc_generator_early_test_function(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	void *data, INT verbose_level);
void placebo_early_test_function(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	void *data, INT verbose_level);
void arc_generator_lifting_prepare_function_new(exact_cover *EC, INT starter_case, 
	INT *candidates, INT nb_candidates, strong_generators *Strong_gens, 
	diophant *&Dio, INT *&col_labels, 
	INT &f_ruled_out, 
	INT verbose_level);
#if 0
void arc_generator_lifting_prepare_function(exact_cover *E, INT starter_case, 
	INT *candidates, INT nb_candidates, strong_generators *Strong_gens, 
	INT *&Inc, INT &nb_free_points, INT &nb_live_blocks2, INT *&live_blocks2, INT &nb_needed, 
	INT &f_has_RHS, INT *&RHS, INT &f_ruled_out, 
	INT verbose_level);
void arc_generator_lifting_cleanup_function(exact_cover *E, INT starter_case, 
	INT *Inc, INT nb_free_points, INT nb_live_blocks2, INT *live_blocks2, 
	INT f_has_RHS, INT *RHS, 
	INT verbose_level);
#endif
void print_arc(ostream &ost, INT len, INT *S, void *data);
void print_point(ostream &ost, INT pt, void *data);
void callback_arc_report(isomorph *Iso, void *data, INT verbose_level);



