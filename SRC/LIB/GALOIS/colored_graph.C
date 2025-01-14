// colored_graph.C
//
// Anton Betten
//
// started:  October 28, 2012




#include "galois.h"

colored_graph::colored_graph()
{
	null();
}

colored_graph::~colored_graph()
{
	freeself();
}

void colored_graph::null()
{
	user_data = NULL;
	points = NULL;
	point_color = NULL;
	bitvector_adjacency = NULL;
	f_has_list_of_edges = FALSE;
	nb_edges = 0;
	list_of_edges = NULL;
}

void colored_graph::freeself()
{
	if (user_data) {
		FREE_INT(user_data);
		}
	if (points) {
		FREE_INT(points);
		}
	if (point_color) {
		FREE_INT(point_color);
		}
	if (bitvector_adjacency) {
		FREE_UBYTE(bitvector_adjacency);
		}
	if (list_of_edges) {
		FREE_INT(list_of_edges);
		}
	null();
}

void colored_graph::compute_edges(INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, j, nb, a;

	if (f_v) {
		cout << "colored_graph::compute_edges" << endl;
		}
	if (f_has_list_of_edges) {
		cout << "colored_graph::compute_edges f_has_list_of_edges" << endl;
		exit(1);
		}
	nb = 0;
	for (i = 0; i < nb_points; i++) {
		for (j = i + 1; j < nb_points; j++) {
			if (is_adjacent(i, j)) {
				nb++;
				}
			}
		}
	list_of_edges = NEW_INT(nb);
	nb_edges = 0;
	for (i = 0; i < nb_points; i++) {
		for (j = i + 1; j < nb_points; j++) {
			if (is_adjacent(i, j)) {
				a = ij2k(i, j, nb_points);
				list_of_edges[nb_edges++] = a;
				}
			}
		}
	if (nb_edges != nb) {
		cout << "colored_graph::compute_edges nb_edges != nb" << endl;
		exit(1);
		}

	f_has_list_of_edges = TRUE;
	if (f_v) {
		cout << "colored_graph::compute_edges done" << endl;
		}
}


INT colored_graph::is_adjacent(INT i, INT j)
{
	if (i == j) {
		return FALSE;
		}
	if (i > j) {
		return is_adjacent(j, i);
		}
	INT k;
	
	k = ij2k(i, j, nb_points);
	return bitvector_s_i(bitvector_adjacency, k);
}

void colored_graph::set_adjacency(INT i, INT j, INT a)
{
	INT k;
	k = ij2k(i, j, nb_points);
	bitvector_m_ii(bitvector_adjacency, k, a);
}

void colored_graph::print()
{
	INT i, j, aij;
	
	cout << "colored graph with " << nb_points << " points and " << nb_colors << " colors" << endl;

	cout << "i : points[i] : point_color[i]" << endl;
	for (i = 0; i < nb_points; i++) {
		cout << i << " : " << points[i] << " : " << point_color[i] << endl;
		}
	
	classify C;

	C.init(point_color, nb_points, FALSE, 0);
	cout << "point colors: ";
	C.print_naked(TRUE);
	cout << endl;

	
	cout << "Adjacency:" << endl;
	for (i = 0; i < nb_points; i++) {
		for (j = 0; j < nb_points; j++) {
			aij = is_adjacent(i, j);
			cout << aij;
			}
		cout << endl;
		}

	INT *A;
	INT I, J, f1, l1, f2, l2, ii, jj, idx1, idx2;

	A = NEW_INT(nb_points * nb_points);
	for (I = 0; I < C.nb_types; I++) {
		f1 = C.type_first[I];
		l1 = C.type_len[I];
		for (J = 0; J < C.nb_types; J++) {
			f2 = C.type_first[J];
			l2 = C.type_len[J];
			cout << "block (" << I << "," << J << ")" << endl;
			for (i = 0; i < l1; i++) {
				ii = f1 + i;
				idx1 = C.sorting_perm_inv[ii];
				for (j = 0; j < l2; j++) {
					jj = f2 + j;
					idx2 = C.sorting_perm_inv[jj];
					aij = is_adjacent(idx1, idx2);
					cout << aij;
					}
				cout << endl;
				}
			cout << endl;
			}
		}
	FREE_INT(A);	
}

void colored_graph::init(INT nb_points, INT nb_colors, 
	INT *colors, UBYTE *bitvec, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;
	
	if (f_v) {
		cout << "colored_graph::init" << endl;
		cout << "nb_points=" << nb_points << endl;
		cout << "nb_colors=" << nb_colors << endl;
		}
	colored_graph::nb_points = nb_points;
	colored_graph::nb_colors = nb_colors;
	
	L = (nb_points * (nb_points - 1)) >> 1;

	bitvector_length = (L + 7) >> 3;

	user_data_size = 0;
	
	points = NEW_INT(nb_points);
	for (i = 0; i < nb_points; i++) {
		points[i] = i;
		}
	point_color = NEW_INT(nb_points);
	INT_vec_copy(colors, point_color, nb_points);
	
	bitvector_adjacency = bitvec;

	if (f_v) {
		cout << "colored_graph::init" << endl;
		}

}

void colored_graph::init_adjacency(INT nb_points, INT nb_colors, 
	INT *colors, INT *Adj, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, j, k;
	INT bitvector_length;
	UBYTE *bitvec;


	if (f_v) {
		cout << "colored_graph::init_adjacency" << endl;
		cout << "nb_points=" << nb_points << endl;
		cout << "nb_colors=" << nb_colors << endl;
		}
	L = (nb_points * (nb_points - 1)) >> 1;

	bitvector_length = (L + 7) >> 3;
	bitvec = NEW_UBYTE(bitvector_length);
	for (i = 0; i < bitvector_length; i++) {
		bitvec[i] = 0;
		}
	for (i = 0; i < nb_points; i++) {
		for (j = i + 1; j < nb_points; j++) {
			if (Adj[i * nb_points + j]) {
				k = ij2k(i, j, nb_points);
				bitvector_m_ii(bitvec, k, 1);
				}
			}
		}
	init(nb_points, nb_colors, 
		colors, bitvec, verbose_level);

	// do not free bitvec here

	if (f_v) {
		cout << "colored_graph::init_adjacency" << endl;
		}

}

void colored_graph::init_adjacency_no_colors(INT nb_points, INT *Adj, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT *vertex_colors;

	if (f_v) {
		cout << "colored_graph::init_adjacency_no_colors" << endl;
		cout << "nb_points=" << nb_points << endl;
		}
	vertex_colors = NEW_INT(nb_points);
	INT_vec_zero(vertex_colors, nb_points);

	init_adjacency(nb_points, 1 /* nb_colors */, 
		vertex_colors, Adj, verbose_level);

	FREE_INT(vertex_colors);
	if (f_v) {
		cout << "colored_graph::init_adjacency_no_colors done" << endl;
		}
}

void colored_graph::init_user_data(INT *data, INT data_size, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;
	
	if (f_v) {
		cout << "colored_graph::init_user_data" << endl;
		}
	user_data_size = data_size;
	user_data = NEW_INT(data_size);
	for (i = 0; i < data_size; i++) {
		user_data[i] = data[i];
		}
	if (f_v) {
		cout << "colored_graph::init_user_data done" << endl;
		}
}

void colored_graph::save(const BYTE *fname, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "colored_graph::save" << endl;
		}

	save_colored_graph(fname, nb_points, nb_colors, 
		points, point_color, 
		user_data, user_data_size, 
		bitvector_adjacency, bitvector_length,
		verbose_level - 1);
		// GALOIS/galois_global.C
	
	if (f_v) {
		cout << "colored_graph::save done" << endl;
		}
}

void colored_graph::load(const BYTE *fname, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	FILE *fp;
	BYTE ext[1000];
	INT i;
	
	if (file_size(fname) <= 0) {
		cout << "colored_graph::load file is empty or does not exist" << endl;
		exit(1);
		}
	
	if (f_v) {
		cout << "colored_graph::load Reading file " << fname << " of size " << file_size(fname) << endl;
		}


	get_extension_if_present(fname, ext);
	strcpy(fname_base, fname);
	chop_off_extension_if_present(fname_base, ext);
	if (f_v) {
		cout << "fname_base=" << fname_base << endl;
		}

	fp = fopen(fname, "rb");

	nb_points = fread_INT4(fp);
	nb_colors = fread_INT4(fp);


	L = (nb_points * (nb_points - 1)) >> 1;

	bitvector_length = (L + 7) >> 3;

	user_data_size = fread_INT4(fp);
	user_data = NEW_INT(user_data_size);
	
	for (i = 0; i < user_data_size; i++) {
		user_data[i] = fread_INT4(fp);
		}

	points = NEW_INT(nb_points);
	point_color = NEW_INT(nb_points);


	if (f_v) {
		cout << "colored_graph::load the graph has " << nb_points << " points and " << nb_colors << " colors" << endl;
		}
	
	for (i = 0; i < nb_points; i++) {
		points[i] = fread_INT4(fp);
		point_color[i] = fread_INT4(fp);
		if (point_color[i] >= nb_colors) {
			cout << "colored_graph::load" << endl;
			cout << "point_color[i] >= nb_colors" << endl;
			cout << "point_color[i]=" << point_color[i] << endl;
			cout << "i=" << i << endl;
			cout << "nb_colors=" << nb_colors << endl;
			exit(1);
			}
		}

	bitvector_adjacency = NEW_UBYTE(bitvector_length);
	fread_UBYTEs(fp, bitvector_adjacency, bitvector_length);


	fclose(fp);
	if (f_v) {
		cout << "colored_graph::load Read file " << fname << " of size " << file_size(fname) << endl;
		}
}

void colored_graph::all_rainbow_cliques(ofstream *fp, INT f_output_solution_raw, 
	INT f_maxdepth, INT maxdepth, 
	INT f_tree, INT f_decision_nodes_only, const BYTE *fname_tree,  
	INT print_interval, 
	INT &search_steps, INT &decision_steps, INT &nb_sol, INT &dt, 
	INT verbose_level)
{
	rainbow_cliques *R;

	R = new rainbow_cliques;
	R->search(this, fp, f_output_solution_raw, 
		f_maxdepth, maxdepth, 
		f_tree, f_decision_nodes_only, fname_tree,  
		print_interval, 
		search_steps, decision_steps, nb_sol, dt, 
		verbose_level);
}

void colored_graph::all_rainbow_cliques_with_additional_test_function(ofstream *fp, INT f_output_solution_raw, 
	INT f_maxdepth, INT maxdepth, 
	INT f_restrictions, INT *restrictions, 
	INT f_tree, INT f_decision_nodes_only, const BYTE *fname_tree,  
	INT print_interval, 
	INT f_has_additional_test_function,
	void (*call_back_additional_test_function)(rainbow_cliques *R, void *user_data, 
		INT current_clique_size, INT *current_clique, 
		INT nb_pts, INT &reduced_nb_pts, 
		INT *pt_list, INT *pt_list_inv, 
		INT verbose_level), 
	INT f_has_print_current_choice_function,
	void (*call_back_print_current_choice)(clique_finder *CF, 
		INT depth, void *user_data, INT verbose_level), 
	void *user_data, 
	INT &search_steps, INT &decision_steps, INT &nb_sol, INT &dt, 
	INT verbose_level)
{
	rainbow_cliques *R;

	R = new rainbow_cliques;
	R->search_with_additional_test_function(this, fp, f_output_solution_raw, 
		f_maxdepth, maxdepth, 
		f_restrictions, restrictions, 
		f_tree, f_decision_nodes_only, fname_tree,  
		print_interval, 
		f_has_additional_test_function,
		call_back_additional_test_function, 
		f_has_print_current_choice_function,
		call_back_print_current_choice, 
		user_data, 
		search_steps, decision_steps, nb_sol, dt, 
		verbose_level);
}

void colored_graph::draw_on_circle(char *fname, 
	INT xmax_in, INT ymax_in, INT xmax_out, INT ymax_out,
	INT f_labels, INT f_embedded, 
	double tikz_global_scale, double tikz_global_line_width)
{
	BYTE fname_full[1000];
	
	sprintf(fname_full, "%s.mp", fname);
	{
	mp_graphics G;
	G.setup(fname, 0, 0, 
		xmax_in /* ONE_MILLION */, ymax_in /* ONE_MILLION */, 
		xmax_out, ymax_out, 
		f_embedded, 
		tikz_global_scale, tikz_global_line_width);
	

	//G.header();
	//G.begin_figure(1000 /* factor_1000 */);
	
	draw_on_circle_2(G, f_labels);


	G.finish(cout, TRUE);
	}
	cout << "written file " << fname_full << " of size " << file_size(fname_full) << endl;
	
}

void colored_graph::draw_on_circle_2(mp_graphics &G, INT f_labels)
{
	INT n = nb_points;
	INT i, j;
	INT *Px, *Py;
	INT *Px1, *Py1;
	double phi = 360. / (double) n;
	INT rad1 = 500000;
	INT rad2 = 5000;
	//BYTE str[1000];
	
	Px = NEW_INT(n);
	Py = NEW_INT(n);
	Px1 = NEW_INT(n);
	Py1 = NEW_INT(n);
	
	for (i = 0; i < n; i++) {
		on_circle_int(Px, Py, i, ((INT)(90. + (double)i * phi)) % 360, rad1);
		//cout << "i=" << i << " Px=" << Px[i] << " Py=" << Py[i] << endl;
		}

	if (f_labels) {
		INT rad_big;

		rad_big = (INT)((double)rad1 * 1.1);
		cout << "rad_big=" << rad_big << endl;
		for (i = 0; i < n; i++) {
			on_circle_int(Px1, Py1, i, ((INT)(90. + (double)i * phi)) % 360, rad_big);
			//cout << "i=" << i << " Px=" << Px[i] << " Py=" << Py[i] << endl;
			}
		}
	for (i = 0; i < n; i++) {

		if (i) {
			//continue;
			}
		
		for (j = i + 1; j < n; j++) {
			if (is_adjacent(i, j)) {
				G.polygon2(Px, Py, i, j);
				}
			}
		}
	for (i = 0; i < n; i++) {
		G.sf_interior(100);
		G.sf_color(0);
		G.circle(Px[i], Py[i], rad2);

		G.sf_interior(0);
		G.sf_color(0);
		G.circle(Px[i], Py[i], rad2);
		}
	if (f_labels) {
		BYTE str[1000];
		for (i = 0; i < n; i++) {
			sprintf(str, "%ld", i);
			G.aligned_text(Px1[i], Py1[i], "", str);
			}
		}
	
	FREE_INT(Px);
	FREE_INT(Py);
	FREE_INT(Px1);
	FREE_INT(Py1);
}



void colored_graph::draw(const BYTE *fname, 
	INT xmax_in, INT ymax_in, INT xmax_out, INT ymax_out,
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_dots = FALSE;
	UBYTE *D = NULL;
	INT len, i, j, k;
	INT nb_vertices;
	
	if (f_v) {
		cout << "colored_graph::draw" << endl;
		}


	nb_vertices = nb_points;

	len = (nb_vertices * nb_vertices + 7) >> 3;
	if (f_v) {
		cout << "colored_graph::draw len = " << len << endl;
		}
	D = NEW_UBYTE(len);
	for (i = 0; i < len; i++) {
		D[i] = 0;
		}
	for (i = 0; i < nb_vertices; i++) {
		for (j = i + 1; j < nb_vertices; j++) {
			k = ij2k(i, j, nb_vertices);
			if (bitvector_s_i(bitvector_adjacency, k)) {
				bitvector_m_ii(D, i * nb_vertices + j, 1);
				bitvector_m_ii(D, j * nb_vertices + i, 1);
				}
			}
		}
	
	draw_bitmatrix(fname, f_dots, 
		FALSE, 0, NULL, 0, NULL, 
		TRUE /* f_bitmatrix */, D, NULL, 
		nb_vertices, nb_vertices, 
		xmax_in, ymax_in, xmax_out, ymax_out, 
		FALSE, NULL);
	

	FREE_UBYTE(D);
	
	if (f_v) {
		cout << "colored_graph::draw done" << endl;
		}
}

void colored_graph::draw_partitioned(const BYTE *fname, 
	INT xmax_in, INT ymax_in, INT xmax_out, INT ymax_out,
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_dots = FALSE;
	UBYTE *D = NULL;
	//INT xmax_out = 1000000;
	//INT ymax_out = 1000000;
	INT len, i, j, k, ii, jj;
	INT nb_vertices;
	
	if (f_v) {
		cout << "colored_graph::draw_partitioned" << endl;
		}


	nb_vertices = nb_points;

	len = (nb_vertices * nb_vertices + 7) >> 3;
	if (f_v) {
		cout << "colored_graph::draw_partitioned len = " << len << endl;
		}
	D = NEW_UBYTE(len);
	for (i = 0; i < len; i++) {
		D[i] = 0;
		}

	classify C;

	C.init(point_color, nb_vertices, FALSE, 0);
	if (f_v) {
		cout << "colored_graph::draw_partitioned we found " << C.nb_types << " classes" << endl;
		}
	
	
	for (i = 0; i < nb_vertices; i++) {
		ii = C.sorting_perm_inv[i];
		for (j = i + 1; j < nb_vertices; j++) {
			jj = C.sorting_perm_inv[j];
			k = ij2k(ii, jj, nb_vertices);
			if (bitvector_s_i(bitvector_adjacency, k)) {
				bitvector_m_ii(D, i * nb_vertices + j, 1);
				bitvector_m_ii(D, j * nb_vertices + i, 1);
				}
			}
		}
	
	INT *part;

	part = NEW_INT(C.nb_types + 1);
	for (i = 0; i < C.nb_types; i++) {
		part[i] = C.type_first[i];
		}
	part[C.nb_types] = nb_vertices;

	draw_bitmatrix(fname, f_dots, 
		TRUE, C.nb_types, part, C.nb_types, part, 
		TRUE /* f_bitmatrix */, D, NULL, 
		nb_vertices, nb_vertices, 
		xmax_in, ymax_in, xmax_out, ymax_out, 
		TRUE /*f_has_labels*/, C.sorting_perm_inv /*labels*/);
		// GALOIS/draw.C

	FREE_UBYTE(D);
	FREE_INT(part);
	
	if (f_v) {
		cout << "colored_graph::draw_partitioned done" << endl;
		}
}

colored_graph *colored_graph::compute_neighborhood_subgraph(INT pt, 
	fancy_set *&vertex_subset, fancy_set *&color_subset, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	colored_graph *S;
	INT *color_in_graph;
	INT *color_in_subgraph;
	INT i, j, l, len, ii, jj, c, idx;
	INT nb_points_subgraph;
	UBYTE *bitvec;

	if (f_v) {
		cout << "colored_graph::compute_neighborhood_subgraph of point " << pt << endl;
		}
	if (f_v) {
		cout << "The graph has " << nb_points << " vertices and " << nb_colors << " colors" << endl;
		}
	S = new colored_graph;
	vertex_subset = new fancy_set;
	color_subset = new fancy_set;
	color_in_graph = NEW_INT(nb_points);
	color_in_subgraph = NEW_INT(nb_points);

	vertex_subset->init(nb_points, 0 /* verbose_level */);
	color_subset->init(nb_colors, 0 /* verbose_level */);
	
	for (i = 0; i < nb_points; i++) {
		if (i == pt) {
			continue;
			}
		if (is_adjacent(i, pt)) {
			c = point_color[i];
			color_in_graph[vertex_subset->k] = c;
			vertex_subset->add_element(i);
			color_subset->add_element(c);
			}
		}


	nb_points_subgraph = vertex_subset->k;

	color_subset->sort();

	if (f_v) {
		cout << "The subgraph has " << nb_points_subgraph << " vertices and " << color_subset->k << " colors" << endl;
		}

	for (i = 0; i < nb_points_subgraph; i++) {
		c = color_in_graph[i];
		if (!INT_vec_search(color_subset->set, color_subset->k, c, idx)) {
			cout << "error, did not find color" << endl;
			exit(1);
			}
		color_in_subgraph[i] = idx;
		}
	
	l = (nb_points_subgraph * (nb_points_subgraph - 1)) >> 1;
	len = (l + 7) >> 3;
	bitvec = NEW_UBYTE(len);
	for (i = 0; i < len; i++) {
		bitvec[i] = 0;
		}
	S->init(nb_points_subgraph, color_subset->k, color_in_subgraph, bitvec, verbose_level);
	for (i = 0; i < nb_points_subgraph; i++) {
		ii = vertex_subset->set[i];
		for (j = i + 1; j < nb_points_subgraph; j++) {
			jj = vertex_subset->set[j];
			if (is_adjacent(ii, jj)) {
				S->set_adjacency(i, j, 1);
				S->set_adjacency(j, i, 1);
				}
			}
		}
	FREE_INT(color_in_graph);
	FREE_INT(color_in_subgraph);
	if (f_v) {
		cout << "colored_graph::compute_neighborhood_subgraph done" << endl;
		}
	return S;
}

colored_graph *colored_graph::compute_neighborhood_subgraph_with_additional_test_function(INT pt, 
	fancy_set *&vertex_subset, fancy_set *&color_subset, 
	INT (*test_function)(colored_graph *CG, INT test_point, INT pt, void *test_function_data, INT verbose_level),
	void *test_function_data, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	colored_graph *S;
	INT *color_in_graph;
	INT *color_in_subgraph;
	INT i, j, l, len, ii, jj, c, idx;
	INT nb_points_subgraph;
	UBYTE *bitvec;

	if (f_v) {
		cout << "colored_graph::compute_neighborhood_subgraph_with_additional_test_function of point " << pt << endl;
		}
	if (f_v) {
		cout << "The graph has " << nb_points << " vertices and " << nb_colors << " colors" << endl;
		}
	S = new colored_graph;
	vertex_subset = new fancy_set;
	color_subset = new fancy_set;
	color_in_graph = NEW_INT(nb_points);
	color_in_subgraph = NEW_INT(nb_points);

	vertex_subset->init(nb_points, 0 /* verbose_level */);
	color_subset->init(nb_colors, 0 /* verbose_level */);
	
	for (i = 0; i < nb_points; i++) {
		if (i == pt) {
			continue;
			}
		if (is_adjacent(i, pt)) {

			if ((*test_function)(this, i, pt, test_function_data, 0 /*verbose_level*/)) {
				c = point_color[i];
				color_in_graph[vertex_subset->k] = c;
				vertex_subset->add_element(i);
				color_subset->add_element(c);
				}
			}
		}


	nb_points_subgraph = vertex_subset->k;

	color_subset->sort();

	if (f_v) {
		cout << "The subgraph has " << nb_points_subgraph << " vertices and " << color_subset->k << " colors" << endl;
		}

	for (i = 0; i < nb_points_subgraph; i++) {
		c = color_in_graph[i];
		if (!INT_vec_search(color_subset->set, color_subset->k, c, idx)) {
			cout << "error, did not find color" << endl;
			exit(1);
			}
		color_in_subgraph[i] = idx;
		}
	
	l = (nb_points_subgraph * (nb_points_subgraph - 1)) >> 1;
	len = (l + 7) >> 3;
	bitvec = NEW_UBYTE(len);
	for (i = 0; i < len; i++) {
		bitvec[i] = 0;
		}
	S->init(nb_points_subgraph, color_subset->k, color_in_subgraph, bitvec, verbose_level);
	for (i = 0; i < nb_points_subgraph; i++) {
		ii = vertex_subset->set[i];
		for (j = i + 1; j < nb_points_subgraph; j++) {
			jj = vertex_subset->set[j];
			if (is_adjacent(ii, jj)) {
				S->set_adjacency(i, j, 1);
				S->set_adjacency(j, i, 1);
				}
			}
		}
	FREE_INT(color_in_graph);
	FREE_INT(color_in_subgraph);
	if (f_v) {
		cout << "colored_graph::compute_neighborhood_subgraph_with_additional_test_function done" << endl;
		}
	return S;
}


void colored_graph::export_to_magma(const BYTE *fname, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, j;
	INT *neighbors;
	INT nb_neighbors;

	if (f_v) {
		cout << "colored_graph::export_to_magma" << endl;
		}
	{
		ofstream fp(fname);

		neighbors = NEW_INT(nb_points);
		fp << "Gamma := Graph< " << nb_points << " | [" << endl;
		for (i = 0; i < nb_points; i++) {


			nb_neighbors = 0;
			for (j = 0; j < nb_points; j++) {
				if (j == i) {
					continue;
					}
				if (is_adjacent(i, j)) {
					neighbors[nb_neighbors++] = j;
					}
				}

			fp << "{";
			for (j = 0; j < nb_neighbors; j++) {
				fp << neighbors[j] + 1;
				if (j < nb_neighbors - 1) {
					fp << ",";
					}
				}
			fp << "}";
			if (i < nb_points - 1) {
				fp << ", " << endl;
				}
			}

		FREE_INT(neighbors);
		
		fp << "]>;" << endl;

//> G := Graph< 9 | [ {4,5,6,7,8,9}, {4,5,6,7,8,9}, {4,5,6,7,8,9},
//>                   {1,2,3,7,8,9}, {1,2,3,7,8,9}, {1,2,3,7,8,9},
//>                   {1,2,3,4,5,6}, {1,2,3,4,5,6}, {1,2,3,4,5,6} ]>;


	}
	cout << "Written file " << fname << " of size " << file_size(fname) << endl;

	if (f_v) {
		cout << "colored_graph::export_to_magma" << endl;
		}
}

void colored_graph::export_to_file(const BYTE *fname, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, j;

	if (f_v) {
		cout << "colored_graph::export_to_file" << endl;
		}
	{
		ofstream fp(fname);

		fp << "[" << endl;
		for (i = 0; i < nb_points; i++) {



			fp << "[";
			for (j = 0; j < nb_points; j++) {
				if (is_adjacent(i, j)) {
					fp << "1";
					}
				else {
					fp << "0";
					}
				if (j < nb_points - 1) {
					fp << ",";
					}
				}
			fp << "]";
			if (i < nb_points - 1) {
				fp << ", " << endl;
				}
			}
		fp << "];" << endl;

		


	}
	cout << "Written file " << fname << " of size " << file_size(fname) << endl;

	if (f_v) {
		cout << "colored_graph::export_to_file" << endl;
		}
}

void colored_graph::export_to_file_matlab(const BYTE *fname, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, j;

	if (f_v) {
		cout << "colored_graph::export_to_file_matlab" << endl;
		}
	{
		ofstream fp(fname);

		fp << "A = [" << endl;
		for (i = 0; i < nb_points; i++) {



			//fp << "[";
			for (j = 0; j < nb_points; j++) {
				if (is_adjacent(i, j)) {
					fp << "1";
					}
				else {
					fp << "0";
					}
				if (j < nb_points - 1) {
					fp << ",";
					}
				}
			//fp << "]";
			if (i < nb_points - 1) {
				fp << "; " << endl;
				}
			}
		fp << "]" << endl;

		


	}
	cout << "Written file " << fname << " of size " << file_size(fname) << endl;

	if (f_v) {
		cout << "colored_graph::export_to_file" << endl;
		}
}


void colored_graph::early_test_func_for_clique_search(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	INT j, a, pt;

	if (f_v) {
		cout << "colored_graph::early_test_func_for_clique_search checking set ";
		print_set(cout, len, S);
		cout << endl;
		cout << "candidate set of size " << nb_candidates << ":" << endl;
		INT_vec_print(cout, candidates, nb_candidates);
		cout << endl;
		}
	if (len == 0) {
		nb_good_candidates = nb_candidates;
		INT_vec_copy(candidates, good_candidates, nb_candidates);
		return;
		}

	pt = S[len - 1];

	nb_good_candidates = 0;
	for (j = 0; j < nb_candidates; j++) {
		a = candidates[j];
		
		if (is_adjacent(pt, a)) {
			good_candidates[nb_good_candidates++] = a;
			}
		} // next j
	
}

void colored_graph::early_test_func_for_coclique_search(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	INT j, a, pt;

	if (f_v) {
		cout << "colored_graph::early_test_func_for_coclique_search checking set ";
		print_set(cout, len, S);
		cout << endl;
		cout << "candidate set of size " << nb_candidates << ":" << endl;
		INT_vec_print(cout, candidates, nb_candidates);
		cout << endl;
		}
	if (len == 0) {
		nb_good_candidates = nb_candidates;
		INT_vec_copy(candidates, good_candidates, nb_candidates);
		return;
		}

	pt = S[len - 1];

	nb_good_candidates = 0;
	for (j = 0; j < nb_candidates; j++) {
		a = candidates[j];
		
		if (!is_adjacent(pt, a)) {
			good_candidates[nb_good_candidates++] = a;
			}
		} // next j
	
}

void colored_graph::early_test_func_for_path_and_cycle_search(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	INT i, j, a, b, pt, x, y;
	INT *v;

	if (f_v) {
		cout << "colored_graph::early_test_func_for_path_and_cycle_search checking set ";
		print_set(cout, len, S);
		cout << endl;
		cout << "candidate set of size " << nb_candidates << ":" << endl;
		INT_vec_print(cout, candidates, nb_candidates);
		cout << endl;
		}
	if (len == 0) {
		nb_good_candidates = nb_candidates;
		INT_vec_copy(candidates, good_candidates, nb_candidates);
		return;
		}

	v = NEW_INT(nb_points);
	INT_vec_zero(v, nb_points);

	pt = S[len - 1];

	for (i = 0; i < len; i++) {
		a = S[i];
		b = list_of_edges[a];
		k2ij(b, x, y, nb_points);
		v[x]++;
		v[y]++;
		}

	nb_good_candidates = 0;
	for (j = 0; j < nb_candidates; j++) {
		a = candidates[j];
		b = list_of_edges[a];
		k2ij(b, x, y, nb_points);
		
		if (v[x] < 2 && v[y] < 2) {
			good_candidates[nb_good_candidates++] = a;
			}
		} // next j
	
	FREE_INT(v);
}

INT colored_graph::is_cycle(INT nb_e, INT *edges, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, a, b, x, y;
	INT *v;
	INT ret = TRUE;

	if (f_v) {
		cout << "colored_graph::is_cycle" << endl;
		}
	v = NEW_INT(nb_points);
	INT_vec_zero(v, nb_points);
	
	for (i = 0; i < nb_e; i++) {
		a = edges[i];
		b = list_of_edges[a];
		k2ij(b, x, y, nb_points);
		v[x]++;
		v[y]++;
		}

	ret = TRUE;
	for (i = 0; i < nb_points; i++) {
		if (v[i] != 0 && v[i] != 2) {
			ret = FALSE;
			break;
			}
		}
	

	FREE_INT(v);	
	if (f_v) {
		cout << "colored_graph::is_cycle done" << endl;
		}
	return TRUE;
}




void colored_graph_all_cliques(const BYTE *fname, INT f_output_solution_raw, 
	INT f_draw, INT xmax_in, INT ymax_in, INT xmax_out, INT ymax_out, 
	INT f_output_fname, const BYTE *output_fname, 
	INT f_maxdepth, INT maxdepth, 
	INT f_tree, INT f_decision_nodes_only, const BYTE *fname_tree,  
	INT print_interval, 
	INT &search_steps, INT &decision_steps, INT &nb_sol, INT &dt, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	colored_graph CG;
	BYTE fname_sol[1000];
	BYTE fname_draw[1000];

	if (f_v) {
		cout << "colored_graph_all_cliques" << endl;
		}
	CG.load(fname, verbose_level - 1);
	if (f_output_fname) {
		sprintf(fname_sol, "%s", output_fname);
		sprintf(fname_draw, "%s_graph", output_fname);
		}
	else {
		sprintf(fname_sol, "%s_sol.txt", CG.fname_base);
		sprintf(fname_draw, "%s_graph", CG.fname_base);
		}

	//CG.print();

	{
	ofstream fp(fname_sol);

	if (f_draw) {
		CG.draw_partitioned(fname_draw, xmax_in, ymax_in, xmax_out, ymax_out, verbose_level);
		}
	CG.all_rainbow_cliques(&fp, f_output_solution_raw, 
		f_maxdepth, maxdepth, 
		f_tree, f_decision_nodes_only, fname_tree,  
		print_interval, 
		search_steps, decision_steps, nb_sol, dt, 
		verbose_level - 1);
		fp << -1 << " " << nb_sol << " " << search_steps 
			<< " " << decision_steps << " " << dt << endl;
	}
	if (f_v) {
		cout << "colored_graph_all_cliques done" << endl;
		}
}

void colored_graph_all_cliques_list_of_cases(INT *list_of_cases, INT nb_cases, INT f_output_solution_raw, 
	INT f_draw, INT xmax_in, INT ymax_in, INT xmax_out, INT ymax_out, 
	const BYTE *fname_template, 
	const BYTE *fname_sol, const BYTE *fname_stats, 
	INT f_maxdepth, INT maxdepth, 
	INT f_prefix, const BYTE *prefix, 
	INT print_interval, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, c;
	INT Search_steps = 0, Decision_steps = 0, Nb_sol = 0, Dt = 0;
	INT search_steps, decision_steps, nb_sol, dt;
	BYTE fname[1000];
	BYTE fname_tmp[1000];

	{
	ofstream fp(fname_sol);
	ofstream fp_stats(fname_stats);
	
	fp_stats << "i,Case,Nb_sol,Nb_vertices,search_steps,decision_steps,dt" << endl;
	for (i = 0; i < nb_cases; i++) {
		colored_graph CG;

		c = list_of_cases[i];
		if (f_v) {
			cout << "colored_graph_all_cliques_list_of_cases case " << i << " / " << nb_cases << " which is " << c << endl;
			}
		sprintf(fname_tmp, fname_template, c);
		if (f_prefix) {
			sprintf(fname, "%s%s", prefix, fname_tmp);
			}
		else {
			strcpy(fname, fname_tmp);
			}
		CG.load(fname, verbose_level - 2);

		//CG.print();

		fp << "# start case " << c << endl;
		if (f_draw) {
			BYTE fname_draw[1000];
			
			if (f_prefix) {
				sprintf(fname_draw, "%s%s_graph", prefix, fname_tmp);
				}
			else {
				sprintf(fname_draw, "%s_graph", fname_tmp);
				}
			cout << "colored_graph_all_cliques_list_of_cases fname_draw=" << fname_draw << endl;
			CG.draw_partitioned(fname_draw, xmax_in, ymax_in, xmax_out, ymax_out, verbose_level);
			
			}

		CG.all_rainbow_cliques(&fp, f_output_solution_raw, 
			f_maxdepth, maxdepth, 
			FALSE, FALSE, NULL, 
			print_interval, 
			search_steps, decision_steps, nb_sol, dt, 
			verbose_level - 1);
		fp << "# end case " << c << " " << nb_sol << " " << search_steps 
				<< " " << decision_steps << " " << dt << endl;
		fp_stats << i << "," << c << "," << nb_sol << "," << CG.nb_points << "," << search_steps << "," << decision_steps << "," << dt << endl;
		Search_steps += search_steps;
		Decision_steps += decision_steps;
		Nb_sol += nb_sol;
		Dt += dt;
		}
	fp << -1 << " " << Nb_sol << " " << Search_steps 
				<< " " << Decision_steps << " " << Dt << endl;
	fp_stats << "END" << endl;
	}
	if (f_v) {
		cout << "colored_graph_all_cliques_list_of_cases done Nb_sol=" << Nb_sol << endl;
		}
}


