// grassmann.C
// 
// Anton Betten
// 7/5/2009
//
//
// 
//
//

#include "galois.h"

grassmann::grassmann()
{
	F = NULL;
	base_cols = NULL;
	coset = NULL;
	M = NULL;
	G = NULL;
}

grassmann::~grassmann()
{
	//cout << "grassmann::~grassmann 1" << endl;
	if (base_cols) {
		FREE_INT(base_cols);
		}
	//cout << "grassmann::~grassmann 2" << endl;
	if (coset) {
		FREE_INT(coset);
		}
	//cout << "grassmann::~grassmann 3" << endl;
	if (M) {
		FREE_INT(M);
		}
	//cout << "grassmann::~grassmann 4" << endl;
	if (G) {
		delete G;
		}
}

void grassmann::init(INT n, INT k, finite_field *F, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	
	grassmann::n = n;
	grassmann::k = k;
	grassmann::F = F;
	q = F->q;
	

	if (f_v) {
		cout << "grassmann::init n=" << n << " k=" << k << " q=" << q << endl;
		}


	base_cols = NEW_INT(n);
	coset = NEW_INT(n);
	M = NEW_INT(k * n);
	if (k > 1) {
		G = new grassmann;
		G->init(n - 1, k - 1, F, verbose_level);
		}
	else {
		G = NULL;
		}
}

INT grassmann::nb_of_subspaces(INT verbose_level)
{
	INT nb;

	nb = generalized_binomial(n, k, q);
	return nb;
}

INT grassmann::nb_points_covered(INT verbose_level)
{
	INT nb;

	nb = generalized_binomial(k, 1, q);
	return nb;
}

void grassmann::points_covered(INT *the_points, INT verbose_level)
{
	INT *v;
	INT *w;
	INT i, nb, a;

	v = NEW_INT(k);
	w = NEW_INT(n);
	nb = nb_points_covered(0 /* verbose_level*/);
	for (i = 0; i < nb; i++) {
		PG_element_unrank_modified(*F, v, 1, k, i);
		F->mult_vector_from_the_left(v, M, w, k, n);
		PG_element_rank_modified(*F, w, 1, n, a);
		the_points[i] = a;
		}
	FREE_INT(v);
	FREE_INT(w);
}

void grassmann::unrank_INT_here(INT *Mtx, INT rk, INT verbose_level)
{
	unrank_INT(rk, verbose_level);
	INT_vec_copy(M, Mtx, k * n);
}

INT grassmann::rank_INT_here(INT *Mtx, INT verbose_level)
{
	INT_vec_copy(Mtx, M, k * n);
	return rank_INT(verbose_level);
}

void grassmann::unrank_INT(INT rk, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT r, h, a, A, nb_free_cols, Q, b, c, i, j;
	
	if (f_v) {
		cout << "unrank_INT " << rk << endl;
		}
	if (k == 0) {
		return;
		}
	// null the first row:
	for (j = 0; j < n; j++) {
		M[j] = 0;
		}
	
	
	// find out the value of h:
	r = rk;
	if (f_v) {
		cout << "r=" << r << endl;
		}
	h = 0;
	while (h < n) {
		a = generalized_binomial(n - h - 1, k - 1, q);
		if (f_v) {
			cout << "[" << n - h - 1 << " choose " << k - 1 << "]_" << q << " = " << a << endl;
			}
		nb_free_cols = n - h - 1 - (k - 1);
		Q = i_power_j(q, nb_free_cols);
		if (f_v) {
			cout << "Q=" << Q << endl;
			}
		A = a * Q;
		if (f_v) {
			cout << "A=" << A << endl;
			}
		if (r < A) {
			break;
			}
		r -= A;
		if (f_v) {
			cout << "r=" << r << endl;
			}
		h++;
		}
	if (h == n) {
		cout << "grassmann::unrank_INT h == n" << endl;
		cout << "h=" << h << endl;
		cout << "r=" << r << endl;
		exit(1);
		}
	
	// now h has been determined
	if (f_v) {
		cout << "grassmann::unrank_INT " << rk << " h=" << h << " nb_free_cols=" << nb_free_cols << endl;
		}
	base_cols[0] = h;
	M[h] = 1;
	
	// find out the coset number b and the rank c of the subspace:	
	b = r / a;
	c = r % a;
	if (f_v) {
		cout << "r=" << r << " coset " << b << " subspace rank " << c << endl;
		}
	
	// unrank the coset:
	if (nb_free_cols) {
		AG_element_unrank(q, coset, 1, nb_free_cols, b);
		}
	if (f_v) {
		cout << "grassmann::unrank_INT coset " << b << " = ";
		INT_vec_print(cout, coset, nb_free_cols);
		cout << endl;
		}
	
	//unrank the subspace (if there is one)
	if (k > 1) {
		G->n = n - h - 1;
		G->unrank_INT(c, verbose_level - 1);
		for (j = 0; j < k - 1; j++) {
			base_cols[j + 1] = G->base_cols[j] + h + 1;
			}
		}
	if (f_v) {
		cout << "grassmann::unrank_INT calling INT_vec_complement n=" << n << " k=" << k << " : ";
		INT_vec_print(cout, base_cols, k);
		cout << endl;
		}
	INT_vec_complement(base_cols, n, k);
	
	// fill in the coset:
	if (k == 1) {
		for (j = 0; j < nb_free_cols; j++) {
			M[h + 1 + j] = coset[j];
			}
		}
	else {
		for (j = 0; j < nb_free_cols; j++) {
			M[h + 1 + G->base_cols[G->k + j]] = coset[j];
			}
		}


	// copy the subspace (rows i=1,..,k-1):
	if (k > 1) {
		for (i = 0; i < G->k; i++) {
		
			// zero beginning:
			for (j = 0; j <= h; j++) {
				M[(1 + i) * n + j] = 0;
				}
			
			// the non-trivial part of the row:
			for (j = 0; j < G->n; j++) {
				M[(1 + i) * n + h + 1 + j] = G->M[i * G->n + j];
				}
			}
		}
	if (f_v) {
		cout << "unrank " << rk << ", we found the matrix" << endl;
		print_integer_matrix_width(cout, M, k, n, n, F->log10_of_q + 1);
		cout << "grassmann::unrank_INT base_cols = ";
		INT_vec_print(cout, base_cols, k);
		cout << endl;
		cout << "grassmann::unrank_INT complement = ";
		INT_vec_print(cout, base_cols + k, n - k);
		cout << endl;
		}
	if (f_v) {
		cout << "unrank_INT " << rk << " finished" << endl;
		}
}

INT grassmann::rank_INT(INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT k1, r, h, a, A, nb_free_cols, Q, b, c, i, j;
	
	r = 0;
	if (f_v) {
		cout << "rank_INT " << endl;
		print_integer_matrix_width(cout, M, k, n, n, F->log10_of_q + 1);
		}
	if (k == 0) {
		return 0;
		}
	k1 = F->Gauss_INT(M, FALSE /*f_special */, TRUE /* f_complete */, base_cols, 
		FALSE /* f_P */, NULL, k, n, n, 0 /* verbose_level */);
	
	if (f_v) {
		cout << "after Gauss:" << endl;
		print_integer_matrix_width(cout, M, k, n, n, F->log10_of_q + 1);
		}
	if (k1 != k) {
		cout << "error, does not have full rank" << endl;
		exit(1);
		}
	if (f_v) {
		cout << "base_cols: ";
		INT_vec_print(cout, base_cols, k);
		cout << endl;
		}
	

	if (f_v) {
		cout << "calling INT_vec_complement n=" << n << " k=" << k << " : ";
		INT_vec_print(cout, base_cols, k);
		cout << endl;
		}
	INT_vec_complement(base_cols, n, k);
	if (f_v) {
		cout << "complement : ";
		INT_vec_print(cout, base_cols + k, n - k);
		cout << endl;
		}

	for (h = 0; h < base_cols[0]; h++) {
		nb_free_cols = n - h - 1 - (k - 1);
		Q = i_power_j(q, nb_free_cols);
		a = generalized_binomial(n - h - 1, k - 1, q);
		A = a * Q;
		r += A;
		}
	nb_free_cols = n - h - 1 - (k - 1);
	a = generalized_binomial(n - h - 1, k - 1, q);
	
	// now h has been determined
	if (f_v) {
		cout << "rank h=" << h << " nb_free_cols=" << nb_free_cols << " r=" << r << endl;
		}

	// copy the subspace (rows i=1,..,k-1):
	if (k > 1) {
		G->n = n - h - 1;
		for (i = 0; i < G->k; i++) {
		
			// the non-trivial part of the row:
			for (j = 0; j < G->n; j++) {
				G->M[i * G->n + j] = M[(1 + i) * n + h + 1 + j];
				}
			}
		}


	// rank the subspace (if there is one)
	if (k > 1) {
		c = G->rank_INT(verbose_level - 1);
		}
	else {
		c = 0;
		}

	// get in the coset:
	if (k == 1) {
		for (j = 0; j < nb_free_cols; j++) {
			coset[j] = M[h + 1 + j];
			}
		}
	else {
		for (j = 0; j < nb_free_cols; j++) {
			coset[j] = M[h + 1 + G->base_cols[G->k + j]];
			}
		}
	// rank the coset:
	if (nb_free_cols) {
		AG_element_rank(q, coset, 1, nb_free_cols, b);
		}
	else {
		b = 0;
		}
	if (f_v) {
		cout << "coset " << b << " = ";
		INT_vec_print(cout, coset, nb_free_cols);
		cout << endl;
		}

		
	// compose the rank from the coset number b and the rank c of the subspace:
	r += b * a + c;	
	if (f_v) {
		cout << "r=" << r << " coset " << b << " subspace rank " << c << endl;
		}
	return r;
}

void grassmann::unrank_longinteger_here(INT *Mtx, longinteger_object &rk, INT verbose_level)
{
	unrank_longinteger(rk, verbose_level);
	INT_vec_copy(M, Mtx, k * n);
}

void grassmann::rank_longinteger_here(INT *Mtx, longinteger_object &rk, INT verbose_level)
{
	INT_vec_copy(Mtx, M, k * n);
	rank_longinteger(rk, verbose_level);
}

void grassmann::unrank_longinteger(longinteger_object &rk, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	longinteger_object r, r1, a, A, mA, Q, b, c;
	longinteger_domain D;
	INT i, j, h, nb_free_cols;
	
	if (f_v) {
		cout << "unrank_longinteger " << rk << endl;
		}
	if (k == 0) {
		return;
		}
	// null the first row:
	for (j = 0; j < n; j++) {
		M[j] = 0;
		}
	
	
	// find out the value of h:
	rk.assign_to(r);
	h = 0;
	while (h < n) {
		D.q_binomial(a, n - h - 1, k - 1, q, 0);
		if (f_v) {
			cout << "[" << n - h - 1 << " choose " << k - 1 << "]_" << q << " = " << a << endl;
			}
		nb_free_cols = n - h - 1 - (k - 1);
		Q.create_i_power_j(q, nb_free_cols);
		D.mult(a, Q, A);
		//A = a * Q;
		if (D.compare(r, A) < 0) {
			break;
			}
		A.assign_to(mA);
		mA.negate();
		D.add(r, mA, r1);
		r.swap_with(r1);
		//r -= A;
		h++;
		}
	if (h == n) {
		cout << "grassmann::unrank_longinteger h == n" << endl;
		exit(1);
		}
	
	// now h has been determined
	if (f_v) {
		cout << "grassmann::unrank_longinteger " << rk << " h=" << h << " nb_free_cols=" << nb_free_cols << endl;
		}
	base_cols[0] = h;
	M[h] = 1;
	
	// find out the coset number b and the rank c of the subspace:	
	D.integral_division(r, a, b, c, 0);
	//b = r / a;
	//c = r % a;
	if (f_v) {
		cout << "r=" << r << " coset " << b << " subspace rank " << c << endl;
		}
	
	// unrank the coset:
	if (nb_free_cols) {
		AG_element_unrank_longinteger(q, coset, 1, nb_free_cols, b);
		}
	if (f_v) {
		cout << "grassmann::unrank_longinteger coset " << b << " = ";
		INT_vec_print(cout, coset, nb_free_cols);
		cout << endl;
		}
	
	//unrank the subspace (if there is one)
	if (k > 1) {
		G->n = n - h - 1;
		G->unrank_longinteger(c, verbose_level - 1);
		for (j = 0; j < k - 1; j++) {
			base_cols[j + 1] = G->base_cols[j] + h + 1;
			}
		}
	if (f_v) {
		cout << "grassmann::unrank_longinteger calling INT_vec_complement n=" << n << " k=" << k << " : ";
		INT_vec_print(cout, base_cols, k);
		cout << endl;
		}
	INT_vec_complement(base_cols, n, k);
	
	// fill in the coset:
	if (k == 1) {
		for (j = 0; j < nb_free_cols; j++) {
			M[h + 1 + j] = coset[j];
			}
		}
	else {
		for (j = 0; j < nb_free_cols; j++) {
			M[h + 1 + G->base_cols[G->k + j]] = coset[j];
			}
		}


	// copy the subspace (rows i=1,..,k-1):
	if (k > 1) {
		for (i = 0; i < G->k; i++) {
		
			// zero beginning:
			for (j = 0; j <= h; j++) {
				M[(1 + i) * n + j] = 0;
				}
			
			// the non-trivial part of the row:
			for (j = 0; j < G->n; j++) {
				M[(1 + i) * n + h + 1 + j] = G->M[i * G->n + j];
				}
			}
		}
	if (f_v) {
		cout << "unrank_longinteger " << rk << ", we found the matrix" << endl;
		print_integer_matrix_width(cout, M, k, n, n, F->log10_of_q + 1);
		cout << "grassmann::unrank_longinteger base_cols = ";
		INT_vec_print(cout, base_cols, k);
		cout << endl;
		cout << "grassmann::unrank_longinteger complement = ";
		INT_vec_print(cout, base_cols + k, n - k);
		cout << endl;
		}
	if (f_v) {
		cout << "unrank_longinteger " << rk << " finished" << endl;
		}
}

void grassmann::rank_longinteger(longinteger_object &r, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	longinteger_object r1, a, A, Q, b, c, tmp1, tmp2;
	longinteger_domain D;
	INT k1, nb_free_cols, h, i, j;
	
	r.create(0);
	if (f_v) {
		cout << "grassmann::rank_longinteger " << endl;
		print_integer_matrix_width(cout, M, k, n, n, F->log10_of_q + 1);
		}
	if (k == 0) {
		return;
		}
	k1 = F->Gauss_INT(M, FALSE /*f_special */, TRUE /* f_complete */, base_cols, 
		FALSE /* f_P */, NULL, k, n, n, 0 /* verbose_level */);
	
	if (f_v) {
		cout << "after Gauss:" << endl;
		print_integer_matrix_width(cout, M, k, n, n, F->log10_of_q + 1);
		}
	if (k1 != k) {
		cout << "error, does not have full rank" << endl;
		exit(1);
		}
	if (f_v) {
		cout << "base_cols: ";
		INT_vec_print(cout, base_cols, k);
		cout << endl;
		}
	

	if (f_v) {
		cout << "calling INT_vec_complement for ";
		INT_vec_print(cout, base_cols, k);
		cout << endl;
		}
	INT_vec_complement(base_cols, n, k);
	if (f_v) {
		cout << "yields: ";
		INT_vec_print(cout, base_cols + k, n - k);
		cout << endl;
		}

	for (h = 0; h < base_cols[0]; h++) {
		nb_free_cols = n - h - 1 - (k - 1);
		Q.create_i_power_j(q, nb_free_cols);
		if (f_v) {
			cout << "create_i_power_j q=" << q << " nb_free_cols=" << nb_free_cols << " yields " << Q << endl;
			}
		D.q_binomial(a, n - h - 1, k - 1, q, 0);
		if (f_v) {
			cout << "q_binomial [" << n - h - 1 << "," << k - 1 << "]_" << q << " = " << a << endl;
			}
		D.mult(a, Q, A);
		D.add(r, A, r1);
		r.swap_with(r1);
		}
	nb_free_cols = n - h - 1 - (k - 1);
	D.q_binomial(a, n - h - 1, k - 1, q, 0);
	if (f_v) {
		cout << "q_binomial [" << n - h - 1 << "," << k - 1 << "]_" << q << " = " << a << endl;
		}
	
	// now h has been determined
	if (f_v) {
		cout << "grassmann::rank_longinteger h=" << h << " nb_free_cols=" << nb_free_cols << " r=" << r << endl;
		}

	// copy the subspace (rows i=1,..,k-1):
	if (k > 1) {
		G->n = n - h - 1;
		for (i = 0; i < G->k; i++) {
		
			// the non-trivial part of the row:
			for (j = 0; j < G->n; j++) {
				G->M[i * G->n + j] = M[(1 + i) * n + h + 1 + j];
				}
			}
		}


	// rank the subspace (if there is one)
	if (k > 1) {
		G->rank_longinteger(c, verbose_level);
		}
	else {
		c.create(0);
		}
	if (f_v) {
		cout << "rank of subspace by induction is " << c << endl;
		}

	// get in the coset:
	if (k == 1) {
		for (j = 0; j < nb_free_cols; j++) {
			coset[j] = M[h + 1 + j];
			}
		}
	else {
		for (j = 0; j < nb_free_cols; j++) {
			coset[j] = M[h + 1 + G->base_cols[G->k + j]];
			}
		}
	// rank the coset:
	if (nb_free_cols) {
		AG_element_rank_longinteger(q, coset, 1, nb_free_cols, b);
		if (f_v) {
			cout << "AG_element_rank_longinteger for coset ";
			INT_vec_print(cout, coset, nb_free_cols);
			cout << " yields " << b << endl;
			}
		}
	else {
		b.create(0);
		}
	if (f_v) {
		cout << "coset " << b << " = ";
		INT_vec_print(cout, coset, nb_free_cols);
		cout << endl;
		}

		
	// compose the rank from the coset number b and the rank c of the subspace:
	if (f_v) {
		cout << "computing r:=" << r << " + " << b << " * " << a << " + " << c << endl;
		}
	D.mult(b, a, tmp1);
	if (f_v) {
		cout << b << " * " << a << " = " << tmp1 << endl;
		}
	D.add(tmp1, c, tmp2);
	if (f_v) {
		cout << tmp1 << " + " << c << " = " << tmp2 << endl;
		}
	D.add(r, tmp2, r1);
	r.swap_with(r1);
	//r += b * a + c;	
	if (f_v) {
		cout << "r=" << r << " coset " << b << " subspace rank " << c << endl;
		}
}

void grassmann::print()
{
	print_integer_matrix_width(cout, M, k, n, n, F->log10_of_q + 1);
}

INT grassmann::dimension_of_join(INT rk1, INT rk2, INT verbose_level)
{
	INT *A;
	INT i, r;

	A = NEW_INT(2 * k * n);
	unrank_INT(rk1, 0);
	for (i = 0; i < k * n; i++) {
		A[i] = M[i];
		}
	unrank_INT(rk2, 0);
	for (i = 0; i < k * n; i++) {
		A[k * n + i] = M[i];
		}
	r = F->rank_of_rectangular_matrix(A, 2 * k, n, 0 /*verbose_level*/);
	return r;
}

void grassmann::unrank_INT_here_and_extend_basis(INT *Mtx, INT rk, INT verbose_level)
// Mtx must be n x n
{
	INT f_v = (verbose_level >= 1);
	INT r, i;
	INT *base_cols;
	INT *embedding;

	if (f_v) {
		cout << "grassmann::unrank_INT_here_and_extend_basis" << endl;
		}
	unrank_INT(rk, verbose_level);
	INT_vec_copy(M, Mtx, k * n);
	base_cols = NEW_INT(n);
	embedding = base_cols + k;
	r = F->base_cols_and_embedding(k, n, Mtx, base_cols, embedding, 0/*verbose_level*/);
	if (r != k) {
		cout << "r != k" << endl;
		exit(1);
		}
	INT_vec_zero(Mtx + k * n, (n - k) * n);
	for (i = 0; i < n - k; i++) {
		Mtx[(k + i) * n + embedding[i]] = 1;
		}
	
	FREE_INT(base_cols);
	if (f_v) {
		cout << "grassmann::unrank_INT_here_and_extend_basis done" << endl;
		}
}





