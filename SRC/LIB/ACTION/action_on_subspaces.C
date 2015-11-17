// action_on_subspaces.C
//
// Anton Betten
// July 20, 2009

#include "galois.h"
#include "action.h"

INT action_on_subspaces::cntr_new = 0;
INT action_on_subspaces::cntr_objects = 0;
INT action_on_subspaces::f_debug_memory = FALSE;

void *action_on_subspaces::operator new(size_t bytes)
{
	cntr_new++;
	cntr_objects++;
	if (f_debug_memory) {
		cout << "action_on_subspaces::operator new bytes=" << bytes
			<< " cntr_new=" << cntr_new
			<< " cntr_objects=" << cntr_objects
			<< endl;
		}
	return malloc(bytes);
}

void *action_on_subspaces::operator new[](size_t bytes)
{
	INT n;

	n = bytes / sizeof(action_on_subspaces);
	cntr_new++;
	cntr_objects += n;
	if (f_debug_memory) {
		cout << "action_on_subspaces::operator new[] n=" << n
			<< " bytes=" << bytes
			<< " cntr_new=" << cntr_new
			<< " cntr_objects=" << cntr_objects
			<< endl;
		}
	return malloc(bytes);
}

void action_on_subspaces::operator delete(void *ptr, size_t bytes)
{
	if (f_debug_memory) {
		cout << "action_on_subspaces::operator delete bytes=" << bytes
			<< " cntr_new=" << cntr_new
			<< " cntr_objects=" << cntr_objects
			<< endl;
		}
	cntr_new--;
	cntr_objects--;
	return ::free(ptr);
}

void action_on_subspaces::operator delete[](void *ptr, size_t bytes)
{
	INT n;

	n = bytes / sizeof(action_on_subspaces);
	if (f_debug_memory) {
		cout << "action_on_subspaces::operator delete[] n=" << n
			<< " cntr_new=" << cntr_new
			<< " cntr_objects=" << cntr_objects
			<< endl;
		}
	cntr_new--;
	cntr_objects -= n;
	return ::free(ptr);
}

action_on_subspaces::action_on_subspaces()
{
	null();
}

action_on_subspaces::~action_on_subspaces()
{
	free();
}

void action_on_subspaces::null()
{
	M = NULL;
	M1 = NULL;
	M2 = NULL;
	S = NULL;
	subspace_basis = NULL;
	subspace_basis2 = NULL;
}

void action_on_subspaces::free()
{
	INT f_v = TRUE;

	if (M1) {
		if (f_v) {
			cout << "action_on_subspaces::free before free M1" << endl;
			}
		FREE_INT(M1);
		}
	if (M2) {
		if (f_v) {
			cout << "action_on_subspaces::free before free M2" << endl;
			}
		FREE_INT(M2);
		}
#if 0
	if (S) {
		delete S;
		}
#endif
	if (subspace_basis) {
		if (f_v) {
			cout << "action_on_subspaces::free before free subspace_basis" << endl;
			}
		FREE_INT(subspace_basis);
		}
	if (subspace_basis2) {
		if (f_v) {
			cout << "action_on_subspaces::free before free subspace_basis2" << endl;
			}
		FREE_INT(subspace_basis2);
		}
	null();
}

void action_on_subspaces::init(action &A, subspaces *S, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	action_on_subspaces::S = S;
	n = S->n;
	q = S->q;
	F = S->F;
	if (f_v) {
		cout << "action_on_subspaces::init" << endl;
		cout << "n=" << n << endl;
		cout << "q=" << q << endl;
		}
	low_level_point_size = n * n;


	M1 = NEW_INT(n * n); // note: make sure this works
	M2 = NEW_INT(n * n);
	if (A.type_G == matrix_group_t) {
		M = A.G.matrix_grp;
		}
	else
	{
		cout<< "invalid group for action on subspaces" <<endl;
		exit(1);
	}
	//else {
		//action *sub = A.subaction;
		//M = sub->G.matrix_grp;
		//}

	if (f_v) {
		cout << "action_on_subspaces::init done" << endl;
		}
}



INT action_on_subspaces::compute_image_INT(action *A, INT *Elt,
	INT i, INT verbose_level)
{
    INT k;
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT h, j;

	if (f_v) {
		cout << "action_on_subspaces::compute_image_INT_ordinary i = " << i << endl;
		cout << "using action " << A->label << endl;
		}

	S->unrank_INT(i, verbose_level - 1);
    k=S->unrank_k(i,verbose_level-1);
	if (f_vv) {
		cout << "action_on_subspaces::compute_image_INT_ordinary after S->unrank_INT" << endl;
		print_integer_matrix_width(cout, S->G[k]->M, k, S->n, S->n, M->GFq->log10_of_q);
		}
	for (h = 0; h < k; h++) {
		A->element_image_of_low_level(S->G[k]->M + h * n, M1 + h * n, Elt, verbose_level - 1);
		}
    INT_vec_copy(M1, S->G[k]->M, k * n);
	j = S->rank_INT(k,verbose_level - 1);
	if (f_v) {
		cout << "action_on_subspaces::compute_image_INT_ordinary image of " << i << " is " << j << endl;
		}
	return j;
}
