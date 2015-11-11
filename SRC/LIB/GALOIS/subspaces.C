// grassmann.C
//
// Jayant Apte
// 11/10/2015
//
//
//
//
//

#include "galois.h"

subspaces::subspaces()
{
	F = NULL;
	base_cols = NULL;
	coset = NULL;
	M = NULL;
	G = NULL;
}

subspaces::~subspaces()
{
	int i;
	//cout << "subspaces::~subspaces 1" << endl;
	if (base_cols) {
		FREE_INT(base_cols);
		}
	//cout << "subspaces::~subspaces 2" << endl;
	if (coset) {
		FREE_INT(coset);
		}
	//cout << "subspaces::~subspaces 3" << endl;
	if (M) {
		FREE_INT(M);
		}
	//cout << "subspaces::~subspaces 4" << endl;
	for (i = 0; i <= n; i++)
	{
	if (G[i]) {
		delete G[i];
		}
	}
	//cout << "subspaces::~subspaces 5" << endl;
	free(G);
}

void subspaces::init(INT n, finite_field *F, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	subspaces::n = n;
	subspaces::F = F;
	q = F->q;


	if (f_v) {
		cout << "subspaces::init n=" << n << " q=" << q << endl;
		}


	base_cols = NEW_INT(n);
	coset = NEW_INT(n);
	G = (grassmann*)malloc((n+1)*sizeof(grassmann));
	for (i = 0; i <= n; i++){
		G[i]->init(n, i, F, verbose_level);
		}
}

INT subspaces::nb_of_subspaces(INT verbose_level)
{
	INT nb,i;
	nb = 0;
	for (i = 0;i <= n; i++) {
		nb = nb + generalized_binomial(n, i, q);
		}
	return nb;
}

INT subspaces::nb_points_covered(INT verbose_level)
{
	INT nb,i;
	nb = 0;
	for (i = 0; i <= n; i++){
		nb = nb + generalized_binomial(i, 1, q);
		}
	return nb;
}


INT subspaces::unrank_k(INT rk, INT verbose_level)
{
	INT nb, nblast, i, nbx;
	nbx=0;
	for (i=0;i<=n;i++)
	{
		nb =  generalized_binomial(n,i,q);
		nbx = nbx + nb;
		if (rk < nbx)
		break;
	}
	return i;
}

// ToDo: void subspaces::points_covered(INT *the_points, INT verbose_level)


void subspaces::unrank_INT_here(INT *Mtx, INT rk, INT verbose_level)
{
	int k=unrank_k(rk,verbose_level);
	G[k]->unrank_INT(rk, verbose_level);
	INT_vec_copy(G[k]->M, Mtx, k * n);
}

INT subspaces::rank_INT_here(INT *Mtx, INT k, INT verbose_level)
{
	INT_vec_copy(Mtx, M, k * n);
	return rank_INT(k, verbose_level);
}

void subspaces::unrank_INT(INT rk, INT verbose_level)
{
	INT k = unrank_k(rk,verbose_level);
	G[k]->unrank_INT(rk, verbose_level);
}

INT subspaces::rank_INT(INT k, INT verbose_level)
{
	return G[k]->rank_INT(verbose_level);
}

void subspaces::print(INT k)
{
	print_integer_matrix_width(cout, G[k]->M, k, n, n, F->log10_of_q + 1);
}

// ToDo: INT subspaces::dimension_of_join(INT rk1, INT rk2, INT verbose_level)
