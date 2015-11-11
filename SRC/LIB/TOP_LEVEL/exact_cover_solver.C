// exact_cover_solver.C
// 
// Anton Betten
//
// started:    4/22/13
// 
//
//

#include "orbiter.h"

exact_cover_solver::exact_cover_solver()
{
	null();
}

exact_cover_solver::~exact_cover_solver()
{
	freeself();
}

void exact_cover_solver::null()
{
	Solutions = NULL;
}

void exact_cover_solver::freeself()
{
	if (Solutions) {
		FREE_INT(Solutions);
		}
	null();
}

void exact_cover_solver::init(INT f_prefix, const BYTE *prefix, 
	BYTE *fname_mask, BYTE *fname_solutions_mask, 
	INT *Inc, INT nb_rows, INT nb_cols, INT *col_labels,
	INT f_has_RHS, INT *RHS, 
	INT case_number, INT number_of_cases, 
	INT nb_needed, INT f_solve, INT f_DLX, INT f_read_from_file, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	

	if (f_v) {
		cout << "exact_cover_solver::init" << endl;
		cout << "f_solve=" << f_solve << endl;
		cout << "f_DLX=" << f_DLX << endl;
		cout << "f_read_from_file=" << f_read_from_file << endl;
		cout << "f_has_RHS=" << f_has_RHS << endl;
		cout << "nb_rows=" << nb_rows << endl;
		cout << "nb_cols=" << nb_cols << endl;
		cout << "nb_needed=" << nb_needed << endl;
		}
	exact_cover_solver::f_prefix = f_prefix;
	exact_cover_solver::prefix = prefix;
	exact_cover_solver::fname_mask = fname_mask;
	exact_cover_solver::fname_solutions_mask = fname_solutions_mask;
	exact_cover_solver::Inc = Inc;
	exact_cover_solver::col_labels = col_labels;
	exact_cover_solver::f_has_RHS = f_has_RHS;
	exact_cover_solver::RHS = RHS;
	exact_cover_solver::nb_rows = nb_rows;
	exact_cover_solver::nb_cols = nb_cols;
	exact_cover_solver::nb_needed = nb_needed;
	exact_cover_solver::f_solve = f_solve;
	exact_cover_solver::f_DLX = f_DLX;
	exact_cover_solver::f_read_from_file = f_read_from_file;
	exact_cover_solver::case_number = case_number;
	exact_cover_solver::number_of_cases = number_of_cases;
	nb_backtrack = 0;
	dt = 0;



	if (f_v) {
		cout << "exact_cover_solver::init() done" << endl;
		}
}

void exact_cover_solver::solve(INT f_draw_system, const BYTE *fname_system, 
	INT f_write_tree, const BYTE *fname_tree, 
	INT verbose_level)
// the solutions  that come back have col_label applied already
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT i, j, a, b;
	
	nb_sol = 0;
	dt = 0;


	if (f_v) {
		cout << "exact_cover_solver::solve() verbose_level=" << verbose_level << endl;
		}

	if (f_read_from_file) {
		read_solution_file(verbose_level - 2);
		}
	else {
		if (f_solve) {
			solve_diophant(Inc, nb_rows, nb_cols, nb_needed, 
				f_has_RHS, RHS, 
				Solutions, nb_sol, nb_backtrack, dt, 
				f_DLX, 
				f_draw_system, fname_system, 
				f_write_tree, fname_tree,  
				verbose_level - 2);
				// TOP_LEVEL/extra.C
				
			if (f_vv) {
				cout << "exact_cover_solver::solve " << case_number << " / " << number_of_cases << " after solve, nb_sol=" << nb_sol << endl;
				}
			}
		else {
			write_problem_to_file(verbose_level - 2);
			if (f_vv) {
				cout << "exact_cover_solver::handle_case " << case_number << " / " << number_of_cases << " written problem to file" << endl;
				}
			}
		}
	if (nb_sol) {
		for (i = 0; i < nb_sol; i++) {
			if (f_vv) {
				cout << "exact_cover_solver::solve " << case_number << " / " << number_of_cases << " solution " << i << " / " << nb_sol << " : ";
				for (j = 0; j < nb_needed; j++) {
					a = Solutions[i * nb_needed + j];
					cout << " " << a;
					}
				cout << endl;
				}
			for (j = 0; j < nb_needed; j++) {
				a = Solutions[i * nb_needed + j];
				b = col_labels[a];
				Solutions[i * nb_needed + j] = b;
				}
			}
		}
}

void exact_cover_solver::write_problem_to_file(INT verbose_level)
{
	BYTE fname[1000];
	
	if (f_prefix) {
		strcpy(fname, prefix);
		}
	else {
		fname[0] = 0;
		}
	sprintf(fname + strlen(fname), fname_mask, case_number);

	write_exact_cover_problem_to_file(Inc, nb_rows, nb_cols, fname);
		// GALOIS/util.C
}

void exact_cover_solver::read_solution_file(INT verbose_level)
{
	BYTE fname[1000];
	INT sol_length;
	
	if (f_prefix) {
		strcpy(fname, prefix);
		}
	else {
		fname[0] = 0;
		}
	sprintf(fname + strlen(fname), fname_solutions_mask, case_number);

	::read_solution_file(fname, 
		Inc, nb_rows, nb_cols, 
		Solutions, sol_length, nb_sol, 
		verbose_level - 1);
}


