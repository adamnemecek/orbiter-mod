// sub_generator.C
//
// Jayant Apte
// November 16, 2015

#include "orbiter.h"
#include "discreta.h"
#include "sub.h"

// ##################################################################################################
// start of class sub_generator
// ##################################################################################################


void sub_generator::read_arguments(int argc, const char **argv)
{
	INT i,j;
	INT f_n = FALSE;
	INT f_klist = FALSE;
	INT f_q = FALSE;
	INT f_d = FALSE;


	if (argc <= 4) {
		print_usage();
		exit(1);
		}

	gen->read_arguments(argc, argv, 0);

	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-lex") == 0) {
			f_lex = TRUE;
			cout << "-lex " << endl;
			}
		else if (strcmp(argv[i], "-debug") == 0) {
			f_debug = TRUE;
			cout << "-debug " << endl;
			}
		else if (strcmp(argv[i], "-schreier_depth") == 0) {
			schreier_depth = atoi(argv[++i]);
			cout << "-schreier_depth " << schreier_depth << endl;
			}
		else if (strcmp(argv[i], "-n") == 0) {
			f_n = TRUE;
			n = atoi(argv[++i]);
			cout << "-n " << n << endl;
			}
		else if (strcmp(argv[i], "-klist") == 0) {
			f_klist = TRUE;
			nk = atoi(argv[++i]);
            klist=NEW_INT(nk);
            for (j=0; j<nk; j++)
            {
                klist[j]=atoi(argv[++i]);
                cout << "-k"<< j+1<<"\t" << klist[j] << endl;
            }
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
		else if (strcmp(argv[i], "-draw_poset") == 0) {
			f_draw_poset = TRUE;
			cout << "-draw_poset " << endl;
			}
		else if (strcmp(argv[i], "-print_data_structure") == 0) {
			f_print_data_structure = TRUE;
			cout << "-print_data_structure " << endl;
			}
		}

	if (!f_n) {
		cout << "Please use option -n <n> to specify n" << endl;
		exit(1);
		}
	if (!f_klist) {
		cout << "Please use option -klist <k_1 k_2...k_i> to specify k values" << endl;
		exit(1);
		}
	if (!f_q) {
		cout << "Please use option -q <q> to specify q" << endl;
		exit(1);
		}
	if (!f_d) {
		cout << "Please use option -d <d> to specify d" << endl;
		exit(1);
		}
	cout << "n=" << n << endl;
	cout << "klist=";
    INT_vec_print(cout,klist,nk);
	cout << "q=" << q << endl;
	cout << "d=" << d << endl;

	f_irreducibility_test = TRUE;

	INT p, h;

	is_prime_power(q, p, h);
	if (h > 1) {
		f_semilinear = TRUE;
		}
	else {
		f_semilinear = FALSE;
		}

}

sub_generator::sub_generator()
{
	null();
}

sub_generator::~sub_generator()
{
	freeself();
}

void sub_generator::null()
{
	gen = NULL;
	F = NULL;
	A = NULL;
	schreier_depth = 1000;
	f_use_invariant_subset_if_available = TRUE;
	f_debug = FALSE;
	f_lex = FALSE;
	f_draw_poset = FALSE;
	f_print_data_structure = FALSE;
}

void sub_generator::freeself()
{
	if (A) {
		delete A;
		}
	if (F) {
		delete F;
		}
	if (gen) {
		delete gen;
		}
    if (klist) {
        FREE_INT(klist);
        }
	null();
}

void sub_generator::init(int argc, const char **argv)
{
	F = new finite_field;
  Asub = new action_on_subspaces;
	A_ind=new action;
	A = new action;
	gen = new generator;
	INT f_basis = TRUE;
    INT i;

	read_arguments(argc, argv);

	INT verbose_level = gen->verbose_level;

	INT f_v = (verbose_level >= 1);
	if (f_v) {
		cout << "sub_generator::init" << endl;
		}
	//nmk = n - k;

    char k_str[1000];
	sprintf(gen->fname_base, "%scodes_%ld_%ld_%ld_%ld ", gen->prefix, n, q, d, nk);
    for (i=0;i<nk;i++)
    {
        sprintf(k_str,"%ld ",klist[i]);
        strcat(gen->fname_base,k_str);
    }
	F->init(q, 0);
	if (f_v) {
		cout << "sub_generator::init calling init_matrix_group" << endl;
		}
	A->init_projective_group(d, F,
		f_semilinear,
		f_basis,
		verbose_level - 2);

	if (f_v) {
		cout << "sub_generator::init finished with init_matrix_group" << endl;
		}
	if (f_v) {
		cout << "sub_generator::init calling init_matrix_group_strong_generators_builtin_projective" << endl;
		}

	A->init_matrix_group_strong_generators_builtin(A->G.matrix_grp, verbose_level);

	if (f_v) {
		cout << "sub_generator::init finished with init_matrix_group_strong_generators_builtin_projective" << endl;
		}

	if (f_v) {
		cout << "sub_generator::init group set up" << endl;
		}

	if (f_v) {
		cout << "arc_generator::init creating action on sets of subspaces" << endl;
		}
    subs = new subspaces;
	subs->init( d/*d*/, F, verbose_level - 2);
    Asub->init(*A, subs, verbose_level - 2);

	A_ind->induced_action_on_subspaces(A, Asub,
    	FALSE /*f_induce_action*/, NULL /*sims *old_G */,
	  MINIMUM(verbose_level - 2, 2));
	if (f_v) {
		cout << "action A_ind created: ";
		A_ind->print_info();
		}


	gen->depth = n;

	if (f_v) {
		cout << "A->f_has_strong_generators=" << A->f_has_strong_generators << endl;
		}

	gen->init(A, A_ind, A->Strong_gens, gen->depth /* sz */, verbose_level);

    if (f_v) {
		cout << "sub_generator::init group set up, calling gen->init_check_func" << endl;
		}

	gen->init_check_func(check_klist, (void*)this /* candidate_check_data */);

	if (f_v) {
		cout << "sub_generator::init group set up, calling gen->init_early_test_func" << endl;
		}
	//gen->init_early_test_func(
		//check_klist_early_test_func,
		//this,
		//verbose_level);
	gen->f_its_OK_to_not_have_an_early_test_func = TRUE;


	//rc.init(F, nmk, n, d);

	gen->f_print_function = TRUE;
	gen->print_function = print_code;
	gen->print_function_data = this;


	INT nb_oracle_nodes = ONE_MILLION;

	if (f_v) {
		cout << "sub_generator::init group set up, calling gen->init_oracle" << endl;
		}

	gen->init_oracle(nb_oracle_nodes, verbose_level - 1);

	if (f_v) {
		cout << "sub_generator::init group set up, calling gen->root[0].init_root_node" << endl;
		}

	gen->root[0].init_root_node(gen, gen->verbose_level - 2);
}


void sub_generator::print(ostream &ost, INT len, INT *S)
{
	INT j,rk_k;

	if (len == 0) {
		return;
		}
	cout << "subspaces:" << endl;
	for (j = 0; j < len; j++) {
		//PG_element_unrank_modified(*F, rc.M1 + j, len /* stride */, nmk /* len */, S[j]);
        rk_k=subs->unrank_k(S[j],-2);
        subs->unrank_INT(S[j],-2);
        print_integer_matrix(ost,subs->G[rk_k]->M , rk_k, n);
        ost << "------------------------------"<<endl;
		}
}


void sub_generator::main()
{
	INT depth;
	INT f_embedded = TRUE;
	INT verbose_level = 0;

	depth = gen->main(t0,
		schreier_depth,
		f_use_invariant_subset_if_available,
		f_lex,
		f_debug,
		gen->verbose_level);


	INT *Table;
	INT nb_rows, nb_cols;

	gen->get_table_of_nodes(Table, nb_rows, nb_cols, verbose_level);

	INT_matrix_write_csv("data.csv", Table, nb_rows, nb_cols);


	FREE_INT(Table);

	if (f_draw_poset) {
		gen->draw_poset(gen->fname_base, depth, 0 /* data1 */, f_embedded, gen->verbose_level);
		}
	if (f_print_data_structure) {
		gen->print_data_structure_tex(depth, gen->verbose_level);
		}
}




// ##################################################################################################
// callback functions
// ##################################################################################################


INT check_klist(INT len, INT *S, void *data, INT verbose_level)
{
	sub_generator *sg = (sub_generator *) data;
    INT i,j;
	INT f_OK = TRUE;
	INT f_v = (verbose_level >= 1);
    INT k_i,f_goodk;
	if (f_v) {
		cout << "checking set ";
		print_set(cout, len, S);
		}
    for (i=0;i<len;i++)
    {
        k_i=sg->subs->unrank_k(S[i], verbose_level - 1);
        f_goodk=FALSE;
        for (j=0;j<sg->nk;j++)
        {
            if(k_i==sg->klist[j]){
            f_goodk=TRUE;
            break;
            }
        }
        if (f_goodk==FALSE)
        {
            f_OK=FALSE;
            break;
        }
    }

	if (f_OK) {
		if (f_v) {
			cout << "OK" << endl;
			}
		return TRUE;
		}
	else {
		return FALSE;
		}
}



void print_code(ostream &ost, INT len, INT *S, void *data)
{
	sub_generator *sg = (sub_generator *) data;

	sg->print(ost, len, S);
}
