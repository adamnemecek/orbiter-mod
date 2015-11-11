// galois_global.C
//
// Anton Betten
//
// started: Oct 16, 2013
//
// unipoly stuff:
// started:  November 16, 2002
// moved here: Oct 16, 2013




#include "galois.h"

void test_unipoly()
{
	finite_field GFp;
	INT p = 2;
	unipoly_object m, a, b, c;
	unipoly_object elts[4];
	INT i, j;
	INT verbose_level = 0;
	
	GFp.init(p, verbose_level);
	unipoly_domain FX(&GFp);
	
	FX.create_object_by_rank(m, 7);
	FX.create_object_by_rank(a, 5);
	FX.create_object_by_rank(b, 55);
	FX.print_object(a, cout); cout << endl;
	FX.print_object(b, cout); cout << endl;

	unipoly_domain Fq(&GFp, m);
	Fq.create_object_by_rank(c, 2);
	for (i = 0; i < 4; i++) {
		Fq.create_object_by_rank(elts[i], i);
		cout << "elt_" << i << " = ";
		Fq.print_object(elts[i], cout); cout << endl;
		}
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			Fq.print_object(elts[i], cout);
			cout << " * ";
			Fq.print_object(elts[j], cout);
			cout << " = ";
			Fq.mult(elts[i], elts[j], c);
			Fq.print_object(c, cout); cout << endl;
			
			FX.mult(elts[i], elts[j], a);
			FX.print_object(a, cout); cout << endl;
			}
		}
	
}

void test_unipoly2()
{
	finite_field Fq;
	INT q = 4, p = 2, i;
	INT verbose_level = 0;
	
	Fq.init(q, verbose_level);
	unipoly_domain FX(&Fq);
	
	unipoly_object a;
	
	FX.create_object_by_rank(a, 0);
	for (i = 1; i < q; i++) {
		FX.minimum_polynomial(a, i, p, TRUE);
		//cout << "minpoly_" << i << " = ";
		//FX.print_object(a, cout); cout << endl;
		}
	
}

BYTE *search_for_primitive_polynomial_of_given_degree(INT p, INT degree, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	finite_field Fp;
	
	Fp.init(p, 0 /*verbose_level*/);
	unipoly_domain FX(&Fp);
	
	unipoly_object m;
	longinteger_object rk;
	
	FX.create_object_by_rank(m, 0);
	
	if (f_v) {
		cout << "search_for_primitive_polynomial_of_given_degree p=" << p << " degree=" << degree << endl;
		}
	FX.get_a_primitive_polynomial(m, degree, verbose_level - 1);
	FX.rank_longinteger(m, rk);
	
	BYTE *s;
	INT i, j;
	if (f_v) {
		cout << "found a polynomial. It's rank is " << rk << endl;
		}
	
	s = NEW_BYTE(rk.len() + 1);
	for (i = rk.len() - 1, j = 0; i >= 0; i--, j++) {
		s[j] = '0' + rk.rep()[i];
		}
	s[rk.len()] = 0;
	
	if (f_v) {
		cout << "created string " << s << endl;
		}
	
	return s;
}


void search_for_primitive_polynomials(INT p_min, INT p_max, INT n_min, INT n_max, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT d, p;
	

	longinteger_f_print_scientific = FALSE;

	
	if (f_v) {
		cout << "search_for_primitive_polynomials p_min=" << p_min << " p_max=" << p_max << " n_min=" << n_min << " n_max=" << n_max << endl;
		}
	for (p = p_min; p <= p_max; p++) {
		if (!is_prime(p)) {
			continue;
			}
		if (f_v) {
			cout << "considering the prime " << p << endl;
			}

			{
			finite_field Fp;
			Fp.init(p, 0 /*verbose_level*/);
			unipoly_domain FX(&Fp);
	
			unipoly_object m;
			longinteger_object rk;
	
			FX.create_object_by_rank(m, 0);
	
			for (d = n_min; d <= n_max; d++) {
				if (f_v) {
					cout << "d=" << d << endl;
					}
				FX.get_a_primitive_polynomial(m, d, verbose_level - 1);
				FX.rank_longinteger(m, rk);
				//cout << d << " : " << rk << " : ";
				cout << "\"" << rk << "\", // ";
				FX.print_object(m, cout);
				cout << endl;
				}
			FX.delete_object(m);
			}
		}
}


void make_linear_irreducible_polynomials(INT q, INT &nb, INT *&table, INT verbose_level)
{
	INT i;
	
	finite_field F;
	F.init(q, 0 /*verbose_level*/);
#if 0
	if (f_no_eigenvalue_one) {
		nb = q - 2;
		table = NEW_INT(nb * 2);
		for (i = 0; i < nb; i++) {
			table[i * 2 + 0] = F.negate(i + 2);
			table[i * 2 + 1] = 1;
			}
		}
	else {
#endif
		nb = q - 1;
		table = NEW_INT(nb * 2);
		for (i = 0; i < nb; i++) {
			table[i * 2 + 0] = F.negate(i + 1);
			table[i * 2 + 1] = 1;
			}
#if 0
		}
#endif
}



void gl_random_matrix(INT k, INT q, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 1);
	INT *M;
	INT *M2;
	finite_field F;
	unipoly_object char_poly;

	if (f_v) {
		cout << "gl_random_matrix" << endl;
		}
	F.init(q, 0 /*verbose_level*/);
	M = NEW_INT(k * k);
	M2 = NEW_INT(k * k);
	
	F.random_invertible_matrix(M, k, verbose_level - 2);

	cout << "Random invertible matrix:" << endl;
	INT_matrix_print(M, k, k);

	
	{
	unipoly_domain U(&F);



	U.create_object_by_rank(char_poly, 0);
		
	U.characteristic_polynomial(M, k, char_poly, verbose_level - 2);

	cout << "The characteristic polynomial is ";
	U.print_object(char_poly, cout);
	cout << endl;

	U.substitute_matrix_in_polynomial(char_poly, M, M2, k, verbose_level);
	cout << "After substitution, the matrix is " << endl;
	INT_matrix_print(M2, k, k);

	U.delete_object(char_poly);

	}
	FREE_INT(M);
	FREE_INT(M2);

}

void save_colored_graph(const BYTE *fname, INT nb_vertices, INT nb_colors, 
	INT *vertex_labels, INT *vertex_colors, 
	INT *data, INT data_sz, 
	UBYTE *bitvector_adjacency, INT bitvector_length,
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	FILE *fp;
	INT i;

	if (f_v) {
		cout << "save_colored_graph" << endl;
		}

	
	fp = fopen(fname, "wb");

	fwrite_INT4(fp, nb_vertices);
	fwrite_INT4(fp, nb_colors);
	fwrite_INT4(fp, data_sz);
	for (i = 0; i < data_sz; i++) {
		fwrite_INT4(fp, data[i]);
		}
	for (i = 0; i < nb_vertices; i++) {
		if (vertex_labels) {
			fwrite_INT4(fp, vertex_labels[i]);
			}
		else {
			fwrite_INT4(fp, 0);
			}
		fwrite_INT4(fp, vertex_colors[i]);
		}
	fwrite_UBYTEs(fp, bitvector_adjacency, bitvector_length);
	fclose(fp);


	if (f_v) {
		cout << "save_colored_graph done" << endl;
		}
}


void load_colored_graph(const BYTE *fname, INT &nb_vertices, INT &nb_colors, 
	INT *&vertex_labels, INT *&vertex_colors, 
	INT *&user_data, INT &user_data_size, 
	UBYTE *&bitvector_adjacency, INT &bitvector_length,
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	FILE *fp;
	INT i, L;

	if (f_v) {
		cout << "load_colored_graph" << endl;
		}
	fp = fopen(fname, "rb");

	nb_vertices = fread_INT4(fp);
	nb_colors = fread_INT4(fp);


	L = (nb_vertices * (nb_vertices - 1)) >> 1;

	bitvector_length = (L + 7) >> 3;

	user_data_size = fread_INT4(fp);
	user_data = NEW_INT(user_data_size);
	
	for (i = 0; i < user_data_size; i++) {
		user_data[i] = fread_INT4(fp);
		}

	vertex_labels = NEW_INT(nb_vertices);
	vertex_colors = NEW_INT(nb_vertices);
	
	for (i = 0; i < nb_vertices; i++) {
		vertex_labels[i] = fread_INT4(fp);
		vertex_colors[i] = fread_INT4(fp);
		if (vertex_colors[i] >= nb_colors) {
			cout << "colored_graph::load" << endl;
			cout << "vertex_colors[i] >= nb_colors" << endl;
			cout << "vertex_colors[i]=" << vertex_colors[i] << endl;
			cout << "i=" << i << endl;
			cout << "nb_colors=" << nb_colors << endl;
			exit(1);
			}
		}

	bitvector_adjacency = NEW_UBYTE(bitvector_length);
	fread_UBYTEs(fp, bitvector_adjacency, bitvector_length);


	fclose(fp);
	if (f_v) {
		cout << "load_colored_graph done" << endl;
		}
}

INT is_diagonal_matrix(INT *A, INT n)
{
	INT i, j;
	
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i == j) {
				continue;
				}
			else {
				if (A[i * n + j]) {
					return FALSE;
					}
				}
			}
		}
	return TRUE;
}

INT is_association_scheme(INT *color_graph, INT n, INT *&Pijk, 
	INT *&colors, INT &nb_colors, INT verbose_level)
// color_graph[n * n]
// added Dec 22, 2010. Originally in BLT_ANALYZE/analyze_plane_invariant.C
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	//INT f_vvv = (verbose_level >= 3);
	INT N;
	INT *M1;
	INT k, i, j;
	INT ret = FALSE;
	
	if (f_v) {
		cout << "is_association_scheme" << endl;
		}
	N = (n * (n - 1)) / 2;
	M1 = NEW_INT(N);
	k = 0;
	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			M1[k++] = color_graph[i * n + j];
			}
		}
	if (k != N) {
		cout << "N=" << N << endl;
		cout << "k=" << k << endl;
		exit(1);
		}

	classify Cl;

	Cl.init(M1, N, FALSE, 0);
	nb_colors = Cl.nb_types + 1;
	colors = NEW_INT(nb_colors);
	colors[0] = color_graph[0];
	for (i = 0; i < Cl.nb_types; i++) {
		colors[i + 1] = Cl.data_sorted[Cl.type_first[i]];
		}

	if (f_vv) {
		cout << "colors (the 0-th color is the diagonal color): ";
		INT_vec_print(cout, colors, nb_colors);
		cout << endl;
		}

	INT C = nb_colors;
	INT *M = color_graph;
	INT pijk, pijk1, u, v, w, u0, v0;
	
	Pijk = NEW_INT(C * C * C);
	INT_vec_zero(Pijk, C * C * C);
	for (k = 0; k < C; k++) {
		for (i = 0; i < C; i++) {
			for (j = 0; j < C; j++) {
				pijk = -1;
				for (u = 0; u < n; u++) {
					for (v = 0; v < n; v++) {
						//if (v == u) continue;
						if (M[u * n + v] != colors[k])
							continue;
						// now: edge (u,v) is colored k
						pijk1 = 0;
						for (w = 0; w < n; w++) {
							//if (w == u)continue;
							//if (w == v)continue;
							if (M[u * n + w] != colors[i])
								continue;
							if (M[v * n + w] != colors[j])
								continue;
							//cout << "i=" << i << " j=" << j << " k=" << k << " u=" << u << " v=" << v << " w=" << w << " increasing pijk" << endl;
							pijk1++;
							} // next w
						//cout << "i=" << i << " j=" << j << " k=" << k << " u=" << u << " v=" << v << " pijk1=" << pijk1 << endl;
						if (pijk == -1) {
							pijk = pijk1;
							u0 = u;
							v0 = v;
							//cout << "u=" << u << " v=" << v << " p_{" << i << "," << j << "," << k << "}=" << Pijk[i * C * C + j * C + k] << endl;
							}
						else {
							if (pijk1 != pijk) {
								//FREE_INT(Pijk);
								//FREE_INT(colors);

								cout << "not an association scheme" << endl;
								cout << "k=" << k << endl;
								cout << "i=" << i << endl;
								cout << "j=" << j << endl;
								cout << "u0=" << u0 << endl;
								cout << "v0=" << v0 << endl;
								cout << "pijk=" << pijk << endl;
								cout << "u=" << u << endl;
								cout << "v=" << v << endl;
								cout << "pijk1=" << pijk1 << endl;
								//exit(1);

								goto done;
								}
							}
						} // next v
					} // next u
				Pijk[i * C * C + j * C + k] = pijk;
				} // next j
			} // next i
		} // next k

	ret = TRUE;

	if (f_v) {
		cout << "it is an association scheme" << endl;


		if (f_v) {
			print_Pijk(Pijk, C);
			}

		if (C == 3 && colors[1] == 0 && colors[2] == 1) {
			INT k, lambda, mu;

			k = Pijk[2 * C * C + 2 * C + 0]; // p220;
			lambda = Pijk[2 * C * C + 2 * C + 2]; // p222;
			mu = Pijk[2 * C * C + 2 * C + 1]; // p221;
			cout << "it is an srg(" << n << "," << k << "," << lambda << "," << mu << ")" << endl;
			}


		}
	

done:
	FREE_INT(M1);
	return ret;
}

void print_Pijk(INT *Pijk, INT nb_colors)
{
	INT i, j, k;
	INT C = nb_colors;
	
	for (k = 0; k < C; k++) {
		INT *Mtx;

		Mtx = NEW_INT(C * C);
		for (i = 0; i < C; i++) {
			for (j = 0; j < C; j++) {
				Mtx[i * C + j] = Pijk[i * C * C + j * C + k];
				}
			}
		cout << "P^{(" << k << ")}=(p_{i,j," << k << "})_{i,j}:" << endl;
		print_integer_matrix_width(cout, Mtx, C, C, C, 3);
		FREE_INT(Mtx);
		}
}




void write_colored_graph(ofstream &ost, BYTE *label, 
	INT point_offset, 
	INT nb_points, 
	INT f_has_adjacency_matrix, INT *Adj, 
	INT f_has_adjacency_list, INT *adj_list, 
	INT f_has_bitvector, UBYTE *bitvector_adjacency, 
	INT f_has_is_adjacent_callback, 
	INT (*is_adjacent_callback)(INT i, INT j, void *data), 
	void *is_adjacent_callback_data, 
	INT f_colors, INT nb_colors, INT *point_color, 
	INT f_point_labels, INT *point_label)
{
	INT i, j, d, w, aij;
	
	cout << "write_graph " << label 
		<< " with " << nb_points 
		<< " points, point_offset=" <<  point_offset
		<< endl;
	w = INT_log10(nb_points);
	cout << "w=" << w << endl;
	ost << "<GRAPH label=\"" << label << "\" num_pts=\"" << nb_points 
		<< "\" f_has_colors=\"" <<  f_colors
		<< "\" num_colors=\"" << nb_colors 
		<< "\" point_offset=\"" <<  point_offset
		<< "\" f_point_labels=\"" <<  f_point_labels
		<< "\">" 
		<< endl;
	for (i = 0; i < nb_points; i++) {
		d = 0;
		for (j = 0; j < nb_points; j++) {
			if (j == i) {
				continue;
				}
			if (f_has_adjacency_matrix) {
				aij = Adj[i * nb_points + j];
				}
			else if (f_has_adjacency_list) {
				INT h;
				if (i < j) {
					h = ij2k(i, j, nb_points);
					}
				else {
					h = ij2k(j, i, nb_points);
					}
				aij = adj_list[h];
				}
			else if (f_has_bitvector) {
				INT h;
				if (i < j) {
					h = ij2k(i, j, nb_points);
					}
				else {
					h = ij2k(j, i, nb_points);
					}
				aij = bitvector_s_i(bitvector_adjacency, h);
				}
			else if (f_has_is_adjacent_callback) {
				aij = (*is_adjacent_callback)(i, j, is_adjacent_callback_data);
				}
			else {
				cout << "write_colored_graph cannot determine adjacency" << endl;
				}

			if (aij) {
				d++;
				}
			}
		ost << setw(w) << i + point_offset << " " << setw(w) << d << " ";
		for (j = 0; j < nb_points; j++) {
			if (j == i) {
				continue;
				}
			if (f_has_adjacency_matrix) {
				aij = Adj[i * nb_points + j];
				}
			else if (f_has_adjacency_list) {
				INT h;
				if (i < j) {
					h = ij2k(i, j, nb_points);
					}
				else {
					h = ij2k(j, i, nb_points);
					}
				aij = adj_list[h];
				}
			else if (f_has_bitvector) {
				INT h;
				if (i < j) {
					h = ij2k(i, j, nb_points);
					}
				else {
					h = ij2k(j, i, nb_points);
					}
				aij = bitvector_s_i(bitvector_adjacency, h);
				}
			else if (f_has_is_adjacent_callback) {
				aij = (*is_adjacent_callback)(i, j, is_adjacent_callback_data);
				}
			else {
				cout << "write_colored_graph cannot determine adjacency" << endl;
				}
			if (aij) {
				ost << setw(w) << j + point_offset << " ";
				}
			}
		ost << endl;
	
		
		}
	
	if (f_colors) {
		ost << endl;
		for (j = 0; j < nb_colors; j++) {
			d = 0;
			for (i = 0; i < nb_points; i++) {
				if (point_color[i] == j)
					d++;
				}
			ost << setw(w) << j + point_offset << " " << setw(w) << d << " ";
			for (i = 0; i < nb_points; i++) {
				if (point_color[i] == j)
					ost << setw(w) << i + point_offset << " ";
				}
			ost << endl;
			}
		}
	
	if (f_point_labels) {
		ost << endl;
		for (i = 0; i < nb_points; i++) {
			ost << setw(w) << i + point_offset << " " 
				<< setw(6) << point_label[i] << endl;
			}
		}

	ost << "</GRAPH>" << endl;
	
}


int str2int(string &str)
{
	int i, res, l;
	
	l = str.length();
	res = 0;
	for (i = 0; i < l; i++) {
		res = (res * 10) + (str[i] - 48);
		}
	return res;
}




