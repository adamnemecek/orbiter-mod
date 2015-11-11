// layered_graph.C
// 
// Anton Betten
// December 30, 2013
//
//
// 
//
//

#include "galois.h"

double norm_of_vector(INT x1, INT x2, INT y1, INT y2);

layered_graph::layered_graph()
{
	null();
}

layered_graph::~layered_graph()
{
	freeself();
}

void layered_graph::null()
{
	nb_nodes_total = 0;
	L = NULL;
	data1 = -1;
}

void layered_graph::freeself()
{
	if (L) {
		delete [] L;
		}
	null();
}

void layered_graph::init(INT nb_layers, INT *Nb_nodes_layer, 
	const BYTE *fname_base, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;
	
	if (f_v) {
		cout << "layered_graph::init" << endl;
		}
	layered_graph::nb_layers = nb_layers;
	strcpy(layered_graph::fname_base, fname_base);
	L = new graph_layer[nb_layers];
	id_of_first_node = 0;
	for (i = 0; i < nb_layers; i++) {
		if (f_v) {
			cout << "layered_graph::init before L[i].init, i=" << i << endl;
			}
		L[i].init(Nb_nodes_layer[i], id_of_first_node, verbose_level);
		id_of_first_node += Nb_nodes_layer[i];
		}
	nb_nodes_total = id_of_first_node;
	if (f_v) {
		cout << "layered_graph::init done" << endl;
		}
}

void layered_graph::place(INT verbose_level)
{
	double dy, dy2;
	INT i;

	dy = 1. / (double) nb_layers;
	dy2 = dy * .5;
	for (i = 0; i < nb_layers; i++) {
		L[i].y_coordinate = 1. - i * dy - dy2;
		L[i].place(verbose_level);
		}
}

void layered_graph::place_with_y_stretch(double y_stretch, INT verbose_level)
{
	double dy, dy2;
	INT i;

	dy = y_stretch / (double) nb_layers;
	dy2 = dy * .5;
	for (i = 0; i < nb_layers; i++) {
		L[i].y_coordinate = 1. - i * dy - dy2;
		L[i].place(verbose_level);
		}
}

void layered_graph::place_with_grouping(INT **Group_sizes, INT *Nb_groups, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	double dy, dy2;
	INT i;

	if (f_v) {
		cout << "layered_graph::place_with_grouping" << endl;
		}
	dy = 1. / (double) nb_layers;
	dy2 = dy * .5;
	for (i = 0; i < nb_layers; i++) {
		if (f_v) {
			cout << "layered_graph::place_with_grouping layer " << i << endl;
			}
		L[i].y_coordinate = 1. - i * dy - dy2;
		L[i].place_with_grouping(Group_sizes[i], Nb_groups[i], verbose_level);
		if (f_v) {
			cout << "layered_graph::place_with_grouping layer " << i << " done" << endl;
			}
		}
	if (f_v) {
		cout << "layered_graph::place_with_grouping done" << endl;
		}
}

void layered_graph::add_edge(INT l1, INT n1, INT l2, INT n2, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT id1, id2;

	if (f_v) {
		cout << "layered_graph::add_edge l1=" << l1 << " n1=" << n1 << " l2=" << l2 << " n2=" << n2 << endl;
		}
	if (n1 < 0) {
		cout << "layered_graph::add_edge n1 is negative, n1=" << n1 << endl;
		}
	if (n2 < 0) {
		cout << "layered_graph::add_edge n2 is negative, n2=" << n2 << endl;
		}
	if (n1 >= L[l1].nb_nodes) {
		cout << "layered_graph::add_edge n1 >= L[l1].nb_nodes" << endl;
		cout << "l1 = " << l1 << endl;
		cout << "n1 = " << n1 << endl;
		cout << "L[l1].nb_nodes = " << L[l1].nb_nodes << endl;
		exit(1);
		}
	id1 = L[l1].Nodes[n1].id;
	if (n2 >= L[l2].nb_nodes) {
		cout << "layered_graph::add_edge n2 >= L[l2].nb_nodes" << endl;
		cout << "l2 = " << l2 << endl;
		cout << "n2 = " << n2 << endl;
		cout << "L[l2].nb_nodes = " << L[l2].nb_nodes << endl;
		exit(1);
		}
	id2 = L[l2].Nodes[n2].id;
	L[l1].Nodes[n1].add_neighbor(l2, n2, id2);
	L[l2].Nodes[n2].add_neighbor(l1, n1, id1);
	if (f_v) {
		cout << "layered_graph::add_edge done" << endl;
		}
}

void layered_graph::add_text(INT l, INT n, const BYTE *text, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "layered_graph::add_text l=" << l << " n=" << n << endl;
		}
	L[l].Nodes[n].add_text(text);
}

void layered_graph::add_data1(INT data, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "layered_graph::add_data1" << endl;
		}
	data1 = data;
}

void layered_graph::add_node_vec_data(INT l, INT n, INT *v, INT len, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "layered_graph::add_node_vec_data l=" << l << " n=" << n << endl;
		}
	L[l].Nodes[n].add_vec_data(v, len);
}

void layered_graph::set_distinguished_element_index(INT l, INT n, INT index, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "layered_graph::set_distinguished_element_index l=" << l << " n=" << n << endl;
		}
	L[l].Nodes[n].set_distinguished_element(index);
}


void layered_graph::add_node_data1(INT l, INT n, INT data, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "layered_graph::add_node_data1 l=" << l << " n=" << n << endl;
		}
	L[l].Nodes[n].add_data1(data);
}

void layered_graph::add_node_data2(INT l, INT n, INT data, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "layered_graph::add_node_data2 l=" << l << " n=" << n << endl;
		}
	L[l].Nodes[n].add_data2(data);
}

void layered_graph::add_node_data3(INT l, INT n, INT data, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "layered_graph::add_node_data3 l=" << l << " n=" << n << endl;
		}
	L[l].Nodes[n].add_data3(data);
}

void layered_graph::draw(const char *fname, INT xmax, INT ymax, INT x_max, INT y_max, INT rad, 
	INT f_circle, INT f_corners, INT f_nodes_empty, 
	INT f_select_layers, INT nb_layer_select, INT *layer_select, 
	INT f_has_draw_begining_callback, 
	void (*draw_begining_callback)(layered_graph *LG, mp_graphics *G, INT x_max, INT y_max, INT dx, INT dy), 
	INT f_has_draw_ending_callback, 
	void (*draw_ending_callback)(layered_graph *LG, mp_graphics *G, INT x_max, INT y_max, INT dx, INT dy), 
	INT f_has_draw_vertex_callback, 
	void (*draw_vertex_callback)(layered_graph *LG, mp_graphics *G, INT layer, INT node, INT x, INT y, INT dx, INT dy), 
	INT f_show_level_info, 
	INT f_embedded,
	INT f_label_edges, 
	double tikz_global_scale, double tikz_global_line_width
	)
{
	INT f_v = FALSE;
	INT x_min = 0; //, x_max = 10000;
	INT y_min = 0; //, y_max = 10000;
	INT factor_1000 = 1000;
	BYTE fname_full[1000];
	double move_out = 0.005;
	INT edge_label = 1;
	
	strcpy(fname_full, fname);
	strcat(fname_full, ".mp");
	{
	mp_graphics G(fname_full, x_min, y_min, x_max, y_max, f_embedded);
	G.out_xmin() = 0;
	G.out_ymin() = 0;
	G.out_xmax() = xmax;
	G.out_ymax() = ymax;
	//cout << "xmax/ymax = " << xmax << " / " << ymax << endl;
	
	G.tikz_global_scale = tikz_global_scale;
	G.tikz_global_line_width = tikz_global_line_width;

	G.header();
	G.begin_figure(factor_1000);

	if (f_corners) {
		G.frame(move_out);
		}
	
	INT i, j, h, id, n, l;
	//INT rad = 50;
	INT rad_x_twice, rad_y_twice;
	INT x, y, x2, y2;
	INT Px[10], Py[10];
	INT threshold = 50000;
	BYTE text[1000];	
	INT xoffset = 3 * rad / 2;
	INT yoffset = 0;
	INT own_id;
	
	rad_x_twice = rad >> 3;
	rad_y_twice = rad >> 3;
	G.user2dev_dist_x(rad_x_twice);
	G.user2dev_dist_y(rad_y_twice);

	G.sl_thickness(30); // 100 is normal
	

	if (f_has_draw_begining_callback) {
		(*draw_begining_callback)(this, &G, x_max, y_max, rad * 4, rad * 4);
		}


	
	// draw the edges first:
	for (i = 0; i < nb_layers; i++) {

		if (f_select_layers) {
			INT idx;
			
			if (!INT_vec_search_linear(layer_select, nb_layer_select, i, idx)) {
				continue;
				}
			
			}

		if (f_v) {
			cout << "layered_graph::draw drawing edges in layer " << i << "  with " << L[i].nb_nodes << " nodes:" << endl;
			}
		
		for (j = 0; j < L[i].nb_nodes; j++) {
			if (f_v) {
				cout << "layered_graph::draw drawing edges in layer " << i << " node " << j << " neighbors = " << L[i].Nodes[j].nb_neighbors << endl;
				}
			if (f_v) {
				cout << "Vertex " << i << " " << j << " at (" << L[i].Nodes[j].x_coordinate << "," << L[i].y_coordinate << ")" << endl;
				}

			if (L[i].nb_nodes > threshold) {
				if (j > 0 && j < L[i].nb_nodes - 1) {
					if (f_v) {
						cout << "skipping node " << j << " in layer " << i << endl;
						}
					continue;
					}
				}
			coordinates(L[i].Nodes[j].id, x_max, y_max, x, y);
			//G.circle(x, y, rad);


			own_id = L[i].Nodes[j].id;

			INT *up;
			INT *down;
			INT nb_up, nb_down;

			up = NEW_INT(L[i].Nodes[j].nb_neighbors);
			down = NEW_INT(L[i].Nodes[j].nb_neighbors);
			nb_up = 0;
			nb_down = 0;

			for (h = 0; h < L[i].Nodes[j].nb_neighbors; h++) {
				id = L[i].Nodes[j].neighbor_list[h];
				if (f_v) {
					cout << "layered_graph::draw drawing edges in layer " << i << " node " << j << " neighbor = " << h << " / " << L[i].Nodes[j].nb_neighbors << " own_id=" << own_id << " id=" << id << endl;
					}
				if (id < own_id) {
					continue;
					}
				find_node_by_id(id, l, n);
				if (f_v) {
					cout << "is in layer " << l << " mode " << n << endl;
					}
				if (f_select_layers) {
					INT idx;
			
					if (!INT_vec_search_linear(layer_select, nb_layer_select, l, idx)) {
						continue;
						}			
					}
				if (l < i) {
					up[nb_up++] = id;
					if (f_v) {
						cout << "added an up link" << endl;
						}
					}
				else {
					down[nb_down++] = id;
					if (f_v) {
						cout << "added a down link" << endl;
						}
					}
				}


			if (f_v) {
				cout << "layered_graph::draw drawing edges, node " << j << ", nb_up = " << nb_up << endl;
				}
			if (nb_up > threshold) {
				if (f_v) {
					cout << "layered_graph::draw drawing edges nb_up > threshold" << endl;
					}
				for (h = 0; h < nb_up; h++) {
					id = up[h];
					find_node_by_id(id, l, n);
					coordinates(id, x_max, y_max, x2, y2);
					if (h > 0 && h < nb_up - 1) {
#if 1
						Px[0] = x;
						Px[1] = (INT)(x + ((double)(x2 - x)) / norm_of_vector(x, x2, y, y2) * rad_x_twice);
						Py[0] = y;
						Py[1] = (INT)(y + ((double)(y2 - y)) / norm_of_vector(x, x2, y, y2) * rad_y_twice);
#endif
						}
					else {
						Px[0] = x;
						Px[1] = x2;
						Py[0] = y;
						Py[1] = y2;
						}
					G.polygon2(Px, Py, 0, 1);

					if (f_label_edges) {
						Px[2] = (Px[0] + Px[1]) >> 1;
						Py[2] = (Py[0] + Py[1]) >> 1;
						sprintf(text, "%ld", edge_label);
						G.aligned_text_with_offset(Px[2], Py[2], xoffset, yoffset, "", text);
						edge_label++;
						}
					}
				}
			else {
				for (h = 0; h < nb_up; h++) {
					id = up[h];
					find_node_by_id(id, l, n);
					coordinates(id, x_max, y_max, x2, y2);
					Px[0] = x;
					Px[1] = x2;
					Py[0] = y;
					Py[1] = y2;
					G.polygon2(Px, Py, 0, 1);
					if (f_label_edges) {
						Px[2] = (Px[0] + Px[1]) >> 1;
						Py[2] = (Py[0] + Py[1]) >> 1;
						sprintf(text, "%ld", edge_label);
						G.aligned_text_with_offset(Px[2], Py[2], xoffset, yoffset, "", text);
						edge_label++;
						}
					if (l > i) {
						if (f_v) {
							cout << "edge " << i << " " << j << " to " << l << " " << n << endl;
							}
						}
					}
				}

			if (f_v) {
				cout << "layered_graph::draw drawing edges, node " << j << ", nb_down = " << nb_down << endl;
				}
			if (nb_down > threshold) {
				if (f_v) {
					cout << "layered_graph::draw drawing edges nb_down > threshold" << endl;
					}
				for (h = 0; h < nb_down; h++) {
					id = down[h];
					find_node_by_id(id, l, n);
					coordinates(id, x_max, y_max, x2, y2);
					if (h > 0 && h < nb_down - 1) {
#if 1
						Px[0] = x;
						Px[1] = x + ((double)(x2 - x)) / norm_of_vector(x, x2, y, y2) * rad_x_twice;
						Py[0] = y;
						Py[1] = y + ((double)(y2 - y)) / norm_of_vector(x, x2, y, y2) * rad_y_twice;
#endif
						}
					else {
						Px[0] = x;
						Px[1] = x2;
						Py[0] = y;
						Py[1] = y2;
						}
					G.polygon2(Px, Py, 0, 1);
					if (f_label_edges) {
						Px[2] = (Px[0] + Px[1]) >> 1;
						Py[2] = (Py[0] + Py[1]) >> 1;
						sprintf(text, "%ld", edge_label);
						G.aligned_text_with_offset(Px[2], Py[2], xoffset, yoffset, "", text);
						edge_label++;
						}
					if (l > i) {
						if (f_v) {
							cout << "edge " << i << " " << j << " to " << l << " " << n << endl;
							}
						}
					}
				}
			else {
				for (h = 0; h < nb_down; h++) {
					id = down[h];
					find_node_by_id(id, l, n);
					coordinates(id, x_max, y_max, x2, y2);
					Px[0] = x;
					Px[1] = x2;
					Py[0] = y;
					Py[1] = y2;
					G.polygon2(Px, Py, 0, 1);
					if (f_label_edges) {
						Px[2] = (Px[0] + Px[1]) >> 1;
						Py[2] = (Py[0] + Py[1]) >> 1;
						sprintf(text, "%ld", edge_label);
						G.aligned_text_with_offset(Px[2], Py[2], xoffset, yoffset, "", text);
						edge_label++;
						}
					if (l > i) {
						if (f_v) {
							cout << "edge " << i << " " << j << " to " << l << " " << n << endl;
							}
						}
					}
				}

			FREE_INT(up);
			FREE_INT(down);


#if 0
			for (h = 0; h < L[i].Nodes[j].nb_neighbors; h++) {
				id = L[i].Nodes[j].neighbor_list[h];
				find_node_by_id(id, l, n);
				coordinates(id, x_max, y_max, x2, y2);
				Px[0] = x;
				Px[1] = x2;
				Py[0] = y;
				Py[1] = y2;
				G.polygon2(Px, Py, 0, 1);
				if (l > i) {
					if (f_v) {
						cout << "edge " << i << " " << j << " to " << l << " " << n << endl;
						}
					}
				}
#endif

			}
		}

	// now draw the vertices:
	for (i = 0; i < nb_layers; i++) {

		if (f_select_layers) {
			INT idx;
			
			if (!INT_vec_search_linear(layer_select, nb_layer_select, i, idx)) {
				continue;
				}
			
			}

		if (f_v) {
			cout << "layered_graph::draw drawing nodes in layer " << i << "  with " << L[i].nb_nodes << " nodes:" << endl;
			}

		if (L[i].nb_nodes > threshold) {
			coordinates(L[i].Nodes[0].id, x_max, y_max, x, y);
			Px[0] = x;
			Py[0] = y;
			coordinates(L[i].Nodes[L[i].nb_nodes - 1].id, x_max, y_max, x, y);
			Px[1] = x;
			Py[1] = y;
			G.polygon2(Px, Py, 0, 1);
			}
		for (j = 0; j < L[i].nb_nodes; j++) {
			if (f_v) {
				cout << "Vertex " << i << " " << j << " at (" << L[i].Nodes[j].x_coordinate << "," << L[i].y_coordinate << ")" << endl;
				}
			if (L[i].nb_nodes > threshold) {
				if (j > 0 && j < L[i].nb_nodes - 1) {
					continue;
					}
				}
			coordinates(L[i].Nodes[j].id, x_max, y_max, x, y);
			if (L[i].Nodes[j].label) {
				G.nice_circle(x, y, rad * 4);

				if (f_nodes_empty) {
					}
				else {
					if (f_has_draw_vertex_callback) {
						(*draw_vertex_callback)(this, &G, i, j, x, y, rad * 4, rad * 4);
						}
					else {
						//G.circle_text(x, y, L[i].Nodes[j].label);
						G.aligned_text(x, y, "", L[i].Nodes[j].label);
						}
					}
				}
			else {
				G.circle(x, y, rad);
				}
			}
		}


	if (f_has_draw_ending_callback) {
		(*draw_ending_callback)(this, &G, x_max, y_max, rad * 4, rad * 4);
		}


	if (f_show_level_info) {
		// draw depth labels at the side:
		coordinates(L[0].Nodes[0].id, x_max, y_max, x, y);
		Px[0] = 1 * rad;
		Py[0] = y + 4 * rad;
		G.aligned_text(Px[0], Py[0], "", "Level");
		for (i = 0; i < nb_layers - 1; i++) {
			coordinates(L[i].Nodes[0].id, x_max, y_max, x, y);
			Px[0] = 2 * rad;
			Py[0] = y;
			coordinates(L[i + 1].Nodes[0].id, x_max, y_max, x, y);
			Px[1] = 2 * rad;
			Py[1] = y;
			G.polygon2(Px, Py, 0, 1);
			}
		for (i = 0; i < nb_layers; i++) {
			coordinates(L[i].Nodes[0].id, x_max, y_max, x, y);
			Px[0] = 1 * rad;
			Py[0] = y;
			Px[1] = 3 * rad;
			Py[1] = y;
			G.polygon2(Px, Py, 0, 1);
			}
		for (i = 0; i < nb_layers; i++) {
			BYTE str[1000];
			
			coordinates(L[i].Nodes[0].id, x_max, y_max, x, y);
			Px[0] = 0;
			Py[0] = y;
			//G.nice_circle(Px[0], Py[0], rad * 4);
			sprintf(str, "%ld", i);
			G.aligned_text(Px[0], Py[0], "", str);
			}
		}




	G.draw_boxes_final();
	G.end_figure();
	G.footer();
	}
	if (f_v) {
		cout << "layered_graph::draw written file " << fname_full << " of size " << file_size(fname_full) << endl;
		}
	
}

void layered_graph::coordinates_direct(double x_in, double y_in, INT x_max, INT y_max, INT &x, INT &y)
{
	x =  (INT)(x_in * x_max);
	y =  (INT)(y_in * y_max);
}

void layered_graph::coordinates(INT id, INT x_max, INT y_max, INT &x, INT &y)
{
	INT l, n;

	find_node_by_id(id, l, n);
	x = (INT)(L[l].Nodes[n].x_coordinate * x_max);
	y = (INT)(L[l].y_coordinate * y_max);
}

void layered_graph::find_node_by_id(INT id, INT &l, INT &n)
{
	INT i, id0;
	
	id0 = 0;
	for (i = 0; i < nb_layers; i++) {
		if (id >= id0 && id < id0 + L[i].nb_nodes) {
			l = i;
			n = id - id0;
			return;
			}
		id0 += L[i].nb_nodes;
		}
	cout << "layered_graph::find_node_by_id did not find node with id " << id << endl;
	exit(1);
}


void layered_graph::write_file(BYTE *fname, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	memory_object M;
	
	if (f_v) {
		cout << "layered_graph::write_file" << endl;
		}
	M.alloc(1024 /* length */, verbose_level - 1);
	M.used_length = 0;
	M.cur_pointer = 0;
	write_memory_object(&M, verbose_level - 1);
	M.write_file(fname, verbose_level - 1);
	if (f_v) {
		cout << "layered_graph::write_file done" << endl;
		}
}

void layered_graph::read_file(const BYTE *fname, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	memory_object M;
	
	if (f_v) {
		cout << "layered_graph::read_file reading file " << fname << " of size " << file_size(fname) << endl;
		}
	M.read_file(fname, verbose_level - 1);
	if (f_v) {
		cout << "layered_graph::read_file read file " << fname << endl;
		}
	M.cur_pointer = 0;
	read_memory_object(&M, verbose_level - 1);
	if (f_v) {
		cout << "layered_graph::read_file done" << endl;
		}
}

void layered_graph::write_memory_object(memory_object *m, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT i;
	
	if (f_v) {
		cout << "layered_graph::write_memory_object" << endl;
		}
	m->write_int(1); // version number of this file format
	if (f_vv) {
		cout << "after m->write_int(1), m->used_length = " << m->used_length << endl;
		}
	m->write_int(nb_layers);
	if (f_vv) {
		cout << "after m->write_int(nb_layers), nb_layers=" << nb_layers << " m->used_length = " << m->used_length << endl;
		}
	m->write_int(nb_nodes_total);
	m->write_int(id_of_first_node);

	//cout << "layered_graph::write_memory_object data1=" << data1 << endl;
	m->write_int(data1);
	for (i = 0; i < nb_layers; i++) {
		L[i].write_memory_object(m, verbose_level - 1);
		}
	m->write_string(fname_base);
	m->write_int(MAGIC_SYNC); // a check to see if the file is not corrupt
	if (f_v) {
		cout << "layered_graph::write_memory_object finished, data size (in bytes) = " << m->used_length << endl;
		}
}

void layered_graph::read_memory_object(memory_object *m, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;
	INT version, a;
	BYTE *p;
	
	if (f_v) {
		cout << "layered_graph::read_memory_object" << endl;
		}

	freeself();
	
	m->read_int(&version); // version number of this file format
	if (version != 1) {
		cout << "layered_graph::read_memory_object unknown version: version = " << version << endl;
		exit(1);
		}
	m->read_int(&nb_layers);
	m->read_int(&nb_nodes_total);
	m->read_int(&id_of_first_node);
	m->read_int(&data1);

	//cout << "layered_graph::read_memory_object data1=" << data1 << endl;
	
	L = new graph_layer[nb_layers];

	for (i = 0; i < nb_layers; i++) {
		L[i].read_memory_object(m, verbose_level - 1);
		}
	
	m->read_string(p);
	strcpy(fname_base, p);
	FREE_BYTE(p);

	m->read_int(&a);
	if (a != MAGIC_SYNC) {
		cout << "layered_graph::read_memory_object unknown the file seems to be corrupt" << endl;
		exit(1);
		}
	if (f_v) {
		cout << "layered_graph::read_memory_object finished" << endl;
		}
}


double norm_of_vector(INT x1, INT x2, INT y1, INT y2)
{
	return sqrt((double)(x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}


