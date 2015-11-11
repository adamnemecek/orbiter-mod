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
	int i,vl,rk,nbx,nb;
	vl=1;
	grassmann* G;
	finite_field* F;
	F=new finite_field;
	F->init(2,0);
	G = new grassmann;
	G->init(5,2,F,vl);//n,k,F,verbose_level
	//for(i = 0; i <= G->nb_of_subspaces(vl); i++)
	// G->unrank_INT(i,vl);
	/*nb=0;
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
	*/
	INT* arr;
	arr=NEW_INT(2);
	G->points_covered(arr,2);
	for(i=0;i<1;i++)
	cout<< "i= "<< i <<" p=" << arr[i] <<endl;
	delete G;
}
