// tree.C
//
// Anton Betten
// February 7, 2003

#include "galois.h"

tree::tree()
{
	root = NULL;
	nb_nodes = 0;
	max_depth = 0;
	path = NULL;
	f_count_leaves = FALSE;
}

tree::~tree()
{
}

#define TREEPATHLEN 10000
#define BUFSIZE_TREE 100000

void tree::init(const BYTE *fname, INT xmax, INT ymax, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 1);
	BYTE *buf;
	BYTE *p_buf;
	INT l, a, i, i_data, nb_nodes;
	BYTE *c_data;
	INT path[TREEPATHLEN];
	
	if (f_v) {
		cout << "reading tree from file " << fname << endl;
		}
	nb_nodes = 0;
	buf = NEW_BYTE(BUFSIZE_TREE);
	{
	ifstream f(fname);
	//f.getline(buf, BUFSIZE_TREE);
	while (TRUE) {
		if (f.eof()) {
			cout << "premature end of file" << endl;
			exit(1);
			}
		f.getline(buf, BUFSIZE_TREE);

		if (f_vv) {
			cout << "read line '" << buf << "'" << endl;
			}

		p_buf = buf;
		if (buf[0] == '#')
			continue;
		s_scan_int(&p_buf, &a);
		if (a == -1)
			break;
		nb_nodes++;
		}
	//s_scan_int(&p_buf, &nb_nodes);
	}
	if (f_v) {
		cout << "found " << nb_nodes << " nodes in file " << fname << endl;
		}
	
	if (f_v) {
		cout << "calling root->init" << endl;
		}
	root = new tree_node;
	root->init(0 /* depth */, NULL, FALSE, 0, FALSE, 0, NULL, verbose_level - 1);
	
	if (f_v) {
		cout << "reading the file again" << endl;
		}
	{
	ifstream f(fname);
	//f.getline(buf, BUFSIZE_TREE);
	while (TRUE) {
		if (f.eof()) {
			cout << "premature end of file" << endl;
			exit(1);
			}
		f.getline(buf, BUFSIZE_TREE);
		p_buf = buf;
		if (buf[0] == '#')
			continue;
		s_scan_int(&p_buf, &l);
		if (l == -1)
			break;
		if (l >= TREEPATHLEN) {
			cout << "tree::init overflow, please increase the value of TREEPATHLEN" << endl;
			cout << "l=" << l << endl;
			exit(1);
			}
		if (f_vv) {
			cout << "reading entry at depth " << l << endl;
			}
		for (i = 0; i < l; i++) {
			s_scan_int(&p_buf, &path[i]);
			}
		s_scan_int(&p_buf, &i_data);
		while (*p_buf == ' ')
			p_buf++;
		c_data = p_buf;
		for (i = 0; c_data[i]; i++) {
			if (c_data[i] == '#') {
				c_data[i] = 0;
				break;
				}
			}
		if (f_vv) {
			cout << "trying to add node: " << buf << endl;
			}
		root->add_node(l, 0, path, i_data, c_data, 0/*verbose_level - 1*/);
		if (f_vv) {
			cout << "node added: " << buf << endl;
			}
		tree::nb_nodes++;
		max_depth = MAXIMUM(max_depth, l);
		}
	}
	if (f_v) {
		cout << "finished adding nodes, max_depth = " << max_depth << endl;
		cout << "tree::nb_nodes=" << tree::nb_nodes << endl;
		}
	
	if (f_vv) {
		root->print_depth_first();
		}
	tree::path = NEW_INT(max_depth + 1);

	INT my_nb_nodes;
	
	compute_DFS_ranks(my_nb_nodes, verbose_level);
	
	root->calc_weight();
	root->place_xy(0, xmax, ymax, max_depth);
	if (f_v) {
		root->print_depth_first();
		}
	FREE_BYTE(buf);

}

void tree::draw(char *fname, INT xmax_in, INT ymax_in, INT xmax, INT ymax, INT rad, 
	INT f_circle, INT f_circletext, INT f_i, INT f_edge_labels, 
	INT f_has_draw_vertex_callback, 
	void (*draw_vertex_callback)(tree *T, mp_graphics *G, INT *v, INT layer, tree_node *N, 
		INT x, INT y, INT dx, INT dy), 
	INT f_embedded, INT f_on_circle, 
	double tikz_global_scale, double tikz_global_line_width
	)
{
	INT x_min = 0, x_max = xmax_in;
	INT y_min = 0, y_max = ymax_in;
	INT factor_1000 = 1000;
	BYTE fname_full[1000];
	
	strcpy(fname_full, fname);

#if 0
	if (f_edge_labels) {
		strcat(fname_full, "e");
		}
#endif


	strcat(fname_full, ".mp");
	{
	mp_graphics G(fname_full, x_min, y_min, x_max, y_max, f_embedded);
	G.out_xmin() = 0;
	G.out_ymin() = 0;
	G.out_xmax() = xmax;
	G.out_ymax() = ymax;
	cout << "xmax/ymax = " << xmax << " / " << ymax << endl;
	
	G.tikz_global_scale = tikz_global_scale;
	G.tikz_global_line_width = tikz_global_line_width;
	G.header();
	G.begin_figure(factor_1000);
	
	
	//G.frame(0.05);
	
#if 0
	INT x = 500000, y;
	calc_y_coordinate(y, 0, max_depth);
	
	if (f_circletext) {
		G.circle_text(x, y, "$\\emptyset$");
		}
	else {
		G.circle(x, y, 5);
		}
#endif

	//root->draw_sideways(G, f_circletext, f_i, FALSE, 10000 - 0, 10000 - 0, max_depth, f_edge_labels);

	INT *radii;
	INT x0, y0;
	
	if (f_on_circle) {
		INT l;

#if 1
		G.sl_thickness(200); // 100 is normal
		circle_center_and_radii(xmax_in, ymax_in, max_depth, x0, y0, radii);
		for (l = 1; l <= max_depth; l++) {
			G.circle(x0, y0, radii[l]);
			}
#endif
		}


	G.sl_thickness(100); // 100 is normal



	leaf_count = 0;


	root->draw_edges(G, rad, f_circle, f_circletext, f_i, FALSE, 0, 0, max_depth, f_edge_labels, f_has_draw_vertex_callback, draw_vertex_callback, this);


	G.sl_thickness(10); // 100 is normal


	root->draw_vertices(G, rad, f_circle, f_circletext, f_i, FALSE, 0, 0, max_depth, f_edge_labels, f_has_draw_vertex_callback, draw_vertex_callback, this);

	if (f_on_circle) {
		FREE_INT(radii);
		}


	G.draw_boxes_final();
	G.end_figure();
	G.footer();
	}
	cout << "written file " << fname_full << " of size " << file_size(fname_full) << endl;
	
}

void tree::circle_center_and_radii(INT xmax, INT ymax, INT max_depth, INT &x0, INT &y0, INT *&rad)
{
	INT l, dy;
	double y;

	x0 = xmax * 0.5;
	y0 = ymax * 0.5;
	rad = NEW_INT(max_depth + 1);
	for (l = 0; l <= max_depth; l++) {
		dy = (INT)((double)ymax / (double)(max_depth + 1));
		y = tree_node_calc_y_coordinate(ymax, l, max_depth);
		y = ymax - y;
		y -= dy * 0.5;
		y /= ymax;
		rad[l] = y * xmax * 0.5;
		}
}

void tree::compute_DFS_ranks(INT &nb_nodes, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT rk;

	if (f_v) {
		cout << "tree::compute_DFS_ranks" << endl;
		}
	rk = 0;
	root->compute_DFS_rank(rk);
	nb_nodes = rk;
	if (f_v) {
		cout << "tree::compute_DFS_ranks done" << endl;
		}
}


