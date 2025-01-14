// sequences.C
// 
// Anton Betten
// August 6, 2015

#include "orbiter.h"

INT distance(INT n, INT *seq1, INT *seq2);


int main(int argc, char **argv)
{
	INT verbose_level = 0;
	INT i, j, N, s;
	INT *Adj;
	INT *seq1;
	INT *seq2;
	INT *seq3;
	INT f_n = FALSE;
	INT n; // sequence length
	INT f_q = FALSE;
	INT q; // alphabet size
	INT f_d = FALSE;
	INT d; // minimum distance
	INT f_set = FALSE;
	INT sz = 0;
	INT set[1000];

	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-v") == 0) {
			verbose_level = atoi(argv[++i]);
			cout << "-v " << verbose_level << endl;
			}
		else if (strcmp(argv[i], "-n") == 0) {
			f_n = TRUE;
			n = atoi(argv[++i]);
			cout << "-n " << n << endl;
			}
		else if (strcmp(argv[i], "-q") == 0) {
			f_q = TRUE;
			q = atoi(argv[++i]);
			cout << "-q " << q << endl;
			}
		else if (strcmp(argv[i], "-d") == 0) {
			f_d = TRUE;
			d = atoi(argv[++i]);
			cout << "-d " << d << endl;
			}
		else if (strcmp(argv[i], "-set") == 0) {
			f_set = TRUE;
			sz = atoi(argv[++i]);
			for (j = 0; j < sz; j++) {
				set[j] = atoi(argv[++i]);
				}
			cout << "-set " << sz;
			INT_vec_print(cout, set, sz);
			cout << endl;
			}
		}

	if (!f_n) {
		cout << "Please use option -n <n>" << endl;
		exit(1);
		}
	if (!f_q) {
		cout << "Please use option -q <q>" << endl;
		exit(1);
		}
	if (!f_d) {
		cout << "Please use option -d <d>" << endl;
		exit(1);
		}
#if 0
	if (!f_set) {
		cout << "Please use option -set <sz> <s_1> ... <s_sz>" << endl;
		exit(1);
		}
#endif

	N = i_power_j(q, n);
	

	Adj = NEW_INT(N * N);
	INT_vec_zero(Adj, N * N);

	seq1 = NEW_INT(n);
	seq2 = NEW_INT(n);
	seq3 = NEW_INT(n);
	
	for (i = 0; i < N; i++) {
		AG_element_unrank(q, seq1, 1, n, i);
		for (j = i + 1; j < N; j++) {
			AG_element_unrank(q, seq2, 1, n, j);

			s = distance(n, seq1, seq2);
			if (s >= d) {
				Adj[i * N + j] = 1;
				Adj[j * N + 1] = 1;
				}
			}
		}

	colored_graph *CG;
	BYTE fname[1000];

	CG = new colored_graph;
	CG->init_adjacency_no_colors(N, Adj, verbose_level);

	sprintf(fname, "Sequences_%ld_%ld_%ld.colored_graph", n, q, d);

	CG->save(fname, verbose_level);

	delete CG;



	if (f_set) {
		cout << "The set of size " << sz << " is ";
		INT_vec_print(cout, set, sz);
		cout << endl;

		INT *D, *Adj2;
		INT a, b;
		
		D = NEW_INT(sz * sz);
		Adj2 = NEW_INT(sz * sz);
		INT_vec_zero(D, sz * sz);
		INT_vec_zero(Adj2, sz * sz);
		
		for (i = 0; i < sz; i++) {
			a = set[i];
			AG_element_unrank(q, seq1, 1, n, a);
			for (j = 0; j < sz; j++) {
				if (j == i) {
					continue;
					}
				b = set[j];
				AG_element_unrank(q, seq2, 1, n, b);

				s = distance(n, seq1, seq2);
				D[i * sz + j] = s;
				}
			}
		cout << "The distance matrix is:" << endl;
		INT_matrix_print(D, sz, sz);
		cout << endl;

		for (i = 0; i < sz; i++) {
			for (j = 0; j < sz; j++) {
				if (j == i) {
					continue;
					}
				if (D[i * sz + j] > d) {
					Adj2[i * sz + j] = 1;
					}
				}
			}

		cout << "The adjacency matrix is:" << endl;
		INT_matrix_print(Adj2, sz, sz);
		cout << endl;


		colored_graph *CG;
		BYTE fname[1000];

		CG = new colored_graph;
		CG->init_adjacency_no_colors(sz, Adj2, verbose_level);

		sprintf(fname, "Sequences_%ld_%ld_%ld_set_%ld.colored_graph", n, q, d, sz);

		CG->save(fname, verbose_level);

		delete CG;


		}




	FREE_INT(Adj);
	FREE_INT(seq1);
	FREE_INT(seq2);
	FREE_INT(seq3);
}

INT distance(INT n, INT *seq1, INT *seq2)
{
	INT h, s;
	
	s = 0;
	for (h = 0; h < n; h++) {
		if (seq1[h] != seq2[h]) {
			s++;
			}
		}
	return s;
}

