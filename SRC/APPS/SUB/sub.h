// sub.h
//
// Jayant Apte
// November 16, 2015
//

typedef class sub_generator sub_generator;

extern INT t0; // the system time when the program started

void print_usage();
int main(int argc, const char **argv);

// ##################################################################################################
// sub_generator:
// ##################################################################################################

class sub_generator {

public:

	INT n;
    INT nk;
    INT* klist;
	INT q;
	INT d;

	//INT nmk; // n - k


	finite_field *F; // F_q

	action *A; // PGL(n - k, q)
	subspaces *subs;
	action_on_subspaces *Asub;
	action *A_ind;
	generator *gen;


	INT f_irreducibility_test;
	INT f_semilinear;

	INT schreier_depth; // = 1000;
	INT f_use_invariant_subset_if_available; // = TRUE;
	INT f_debug; // = FALSE;
	INT f_lex; // = FALSE;

	INT f_draw_poset;
	INT f_print_data_structure;

	//rank_checker rc;

	void read_arguments(int argc, const char **argv);
	sub_generator();
	~sub_generator();
	void null();
	void freeself();
	void init(int argc, const char **argv);
	void print(ostream &ost, INT len, INT *S);
	void main();
};

void check_klist_early_test_func(INT *S, INT len,
    INT *candidates, INT nb_candidates,
    INT *good_candidates, INT &nb_good_candidates,
    void *data, INT verbose_level);
INT check_klist(INT len, INT *S, void *data, INT verbose_level);
void print_code(ostream &ost, INT len, INT *S, void *data);
