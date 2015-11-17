// subspace_generator.C
//
// Anton Betten
//
// moved here from subspace.C: May 18, 2009
//
// December 30, 2003

#include "orbiter.h"
#include "discreta.h"
 
int main(int argc, const char **argv)
{
	INT i,vl,rk,nbx,nb,k;
	vl=1;
	INT* arr;
	subspaces* S;
	grassmann* G;
	G=new grassmann;
	S=new subspaces;
	finite_field* F;
	F=new finite_field;
	F->init(2,0);
	G->init(5,1,F,3);
	S->init(5,F,3);
	cout<< "total subspaces= "<<S->nb_of_subspaces(3) <<endl;
    action_on_grassmannian* agrass;
    agrass=new action_on_grassmannian;
    action* A;
    A=new action;
    A->init_projective_group(5,F,FALSE,TRUE,3);
    A->init_matrix_group_strong_generators_builtin(A->G.matrix_grp, 3);
    agrass->init(*A,S->G[2],3);
    //for (i=0;i<S->nb_of_subspaces(3);i++)
	//{
		//cout<< ">>bound "<<S->nb_of_subspaces(3)<<endl;
		//cout<< ">>unranking "<<i<<endl;
		//S->unrank_INT(i,3);
        //S->print(S->unrank_k(i,3));
		//cout<< ">>ranking back"<<endl;
		//cout<<">>i = "<< i << " rank ="<< S->rank_INT(S->unrank_k(i,3),3)<< endl;

	//}
    delete A;
    delete G;
    delete F;
	delete S;
}

/* test grassmann.C

INT i,vl,rk,nbx,nb;
vl=1;
INT* arr;
grassmann* G;
finite_field* F;
F=new finite_field;
F->init(2,0);
G = new grassmann;
G->init(5,2,F,vl);//n,k,F,verbose_level
//for(i = 0; i <= G->nb_of_subspaces(vl); i++)
// G->unrank_INT(i,vl);
nb=0;
for(i=0;i<=5;i++)
{
	cout << i<<"  " << nb <<endl;
	nb=nb+generalized_binomial(5,i,2);
}
nbx=0;
rk=373;
for (i=0;i<=5;i++)
{
	nb =  generalized_binomial(5,i,2);
	nbx = nbx + nb;
	if (rk < nbx)
	break;
}
cout << "i" << i <<endl;
delete G;
arr=NEW_INT(2);
G->points_covered(arr,2);
for(i=0;i<1;i++)
cout<< "i= "<< i <<" p=" << arr[i] <<endl;
*/
