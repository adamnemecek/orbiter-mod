// vector_ge.C
//
// Anton Betten
// December 9, 2003

#include "galois.h"
#include "action.h"

#undef PRINT_WITH_TYPE
#define RANGE_CHECKING

INT vector_ge::cntr_new = 0;
INT vector_ge::cntr_objects = 0;
INT vector_ge::f_debug_memory = FALSE;
INT vector_ge::allocation_id = 0;
void *vector_ge::allocated_objects = NULL;

void *vector_ge::operator new(size_t bytes)
{
	void *p;
	void **pp;
	INT *pi;
	
	allocation_id++;
	cntr_new++;
	cntr_objects++;
	if (f_debug_memory) {
		cout << "vector_ge::operator new allocation_id=" << allocation_id 
			<< " bytes=" << bytes 
			<< " cntr_new=" << cntr_new 
			<< " cntr_objects=" << cntr_objects 
			<< endl;
		}
	p = malloc(bytes + 2 * sizeof(void *) + sizeof(INT));
	pp = (void **) p;
	pp[0] = NULL;
	pp[1] = allocated_objects;
	if (allocated_objects) {
		void **next;
		next = (void **)allocated_objects;
		next[0] = p;
		}
	allocated_objects = p;
	pi = (INT *)&pp[2];
	*pi =  allocation_id;
	return pi + 1;
}

void *vector_ge::operator new[](size_t bytes)
{
	void *p;
	void **pp;
	INT *pi;
	INT n;
	
	n = bytes / sizeof(vector_ge);
	allocation_id++;
	cntr_new++;
	cntr_objects += n;
	if (f_debug_memory) {
		cout << "vector_ge::operator new[] n=" << n 
			<< " allocation_id=" << allocation_id 
			<< " bytes=" << bytes 
			<< " cntr_new=" << cntr_new 
			<< " cntr_objects=" << cntr_objects 
			<< endl;
		}
	p = malloc(bytes + 2 * sizeof(void *) + sizeof(INT));
	pp = (void **) p;
	pp[0] = NULL;
	pp[1] = allocated_objects;
	if (allocated_objects) {
		void **next;
		next = (void **)allocated_objects;
		next[0] = p;
		}
	allocated_objects = p;
	pi = (INT *)&pp[2];
	*pi =  allocation_id;
	return pi + 1;
}

void vector_ge::operator delete(void *ptr, size_t bytes)
{
	void **pp;
	void **prev;
	void **next;
	INT *pi;

	pi = (INT *)ptr;
	pi--;
	pp = (void **)pi;
	pp--;
	pp--;
	prev = (void **)pp[0];
	next = (void **)pp[1];
	if (prev == NULL) {
		allocated_objects = next;
		}
	else {
		prev[1] = next;
		}
	if (next == NULL) {
		}
	else {
		next[0] = prev;
		}
	if (f_debug_memory) {
		cout << "vector_ge::operator delete allocation_id=" << *pi 
			<< " bytes=" << bytes 
			<< " cntr_new=" << cntr_new 
			<< " cntr_objects=" << cntr_objects 
			<< endl;
		}
	cntr_new--;
	cntr_objects--;
	return free(pp);
}

void vector_ge::operator delete[](void *ptr, size_t bytes)
{
	INT n;
	n = bytes / sizeof(vector_ge);
	void **pp;
	void **prev;
	void **next;
	INT *pi;

	pi = (INT *)ptr;
	pi--;
	pp = (void **)pi;
	pp--;
	pp--;
	prev = (void **)pp[0];
	next = (void **)pp[1];
	if (prev == NULL) {
		allocated_objects = next;
		}
	else {
		prev[1] = next;
		}
	if (next == NULL) {
		}
	else {
		next[0] = prev;
		}
	if (f_debug_memory) {
		cout << "vector_ge::operator delete allocation_id=" << *pi 
			<< " n=" << n 
			<< " bytes=" << bytes 
			<< " cntr_new=" << cntr_new 
			<< " cntr_objects=" << cntr_objects 
			<< endl;
		}
	cntr_new--;
	cntr_objects -= n;
	return free(pp);
}


vector_ge::vector_ge()
{
	null();
}

vector_ge::vector_ge(action *A)
{
	null();
	vector_ge::A = A;
}

vector_ge::~vector_ge()
{
	//cout << "in ~vector_ge data = " << data << endl;
	freeself();
}

void vector_ge::null()
{
	vector_ge::A = NULL;
	data = NULL;
	len = 0;
}

void vector_ge::freeself()
{
	if (data) {
		FREE_INT(data);
		data = NULL;
		}
}

void vector_ge::init(action *A)
{
	//cout << "vector_ge::init()" << endl;
	freeself();
	vector_ge::A = A;
	data = NULL;
	len = 0;
}

void vector_ge::init_by_hdl(action *A, INT *gen_hdl, INT nb_gen)
{
	INT i;
	
	init(A);
	allocate(nb_gen);
	for (i = 0; i < nb_gen; i++) {
		A->element_retrieve(gen_hdl[i], ith(i), 0);
		}
}

void vector_ge::init_single(action *A, INT *Elt)
{
	init(A);
	allocate(1);
	A->element_move(Elt, ith(0), 0);
}

INT *vector_ge::ith(INT i)
{
#ifdef RANGE_CHECKING
	if (i < 0 || i >= len) {
		cout << "vector_ge::ith() access error i = " << i << " len = " << len << endl;
		exit(1);
		}
#endif
	return data + i * A->elt_size_in_INT; 
}

ostream& vector_ge::print(ostream& ost)
{
	INT i;
	
	ost << "(" << endl;
	//ost << "len=" << len << " A->elt_size_in_INT=" << A->elt_size_in_INT << " data=" << data << endl;
	for (i = 0; i < len; i++) {
		if (data == NULL) {
			cout << "vector_ge::print fatal: data == NULL" << endl;
			exit(1);
			}
		A->element_print(ith(i), ost);
		if (i < len - 1) {
			ost << ", " << endl;
			}
		}
	ost << ")" << endl;
	return ost;
}

ostream& vector_ge::print_quick(ostream& ost)
{
	INT i;
	
	ost << "(" << endl;
	//ost << "len=" << len << " A->elt_size_in_INT=" << A->elt_size_in_INT << " data=" << data << endl;
	for (i = 0; i < len; i++) {
		if (data == NULL) {
			cout << "vector_ge::print fatal: data == NULL" << endl;
			exit(1);
			}
		A->element_print_quick(ith(i), ost);
		if (i < len - 1) {
			ost << ", " << endl;
			}
		}
	ost << ")" << endl;
	return ost;
}

ostream& vector_ge::print_tex(ostream& ost)
{
	INT i;
	
	//ost << "(" << endl;
	//ost << "len=" << len << " A->elt_size_in_INT=" << A->elt_size_in_INT << " data=" << data << endl;
	for (i = 0; i < len; i++) {
		if (data == NULL) {
			cout << "vector_ge::print fatal: data == NULL" << endl;
			exit(1);
			}
		A->element_print_latex(ith(i), ost);
		if (i < len - 1) {
			ost << ", " << endl;
			}
		}
	//ost << ")" << endl;
	return ost;
}

ostream& vector_ge::print_as_permutation(ostream& ost)
{
	INT i;
	
	ost << "(" << endl;
	for (i = 0; i < len; i++) {
		A->element_print(ith(i), ost);
		cout << "is the permutation" << endl;
		A->element_print_as_permutation(ith(i), ost);
		if (i < len - 1) {
			ost << ", " << endl;
			}
		}
	ost << ")" << endl;
	return ost;
}

void vector_ge::allocate(INT length)
{
	if (data) {
		FREE_INT(data);
		//cout << "vector_ge::allocate warning, data != NULL, we seem to be having a memory leak here" << endl;
		}
	len = length;
	data = NEW_INT(length * A->elt_size_in_INT);
}

void vector_ge::reallocate(INT new_length)
{
	INT *data2 = NEW_INT(new_length * A->elt_size_in_INT);
	INT *elt, *elt2, i;
	
	for (i = 0; i < len; i++) {
		elt = ith(i);
		elt2 = data2 + i * A->elt_size_in_INT;
		A->element_move(elt, elt2, FALSE);
		}
	if (data) {
		FREE_INT(data);
		data = NULL;
		}
	data = data2;
	len = new_length;
};

void vector_ge::reallocate_and_insert_at(INT position, INT *elt)
{
	INT *data2 = NEW_INT((len + 1) * A->elt_size_in_INT);
	INT *elt1, *elt2, i;
	
	for (i = 0; i < len; i++) {
		elt1 = ith(i);
		if (i >= position) {
			elt2 = data2 + (i + 1) * A->elt_size_in_INT;
			}
		else {
			elt2 = data2 + i * A->elt_size_in_INT;
			}
		A->element_move(elt1, elt2, FALSE);
		}
	if (data) {
		FREE_INT(data);
		data = NULL;
		}
	data = data2;
	len = len + 1;

	copy_in(position, elt);
}

void vector_ge::insert_at(INT length_before, INT position, INT *elt)
// does not reallocate, but shifts elements up to make space.
// the last element might be lost if there is no space.
{
	INT *elt1, *elt2, i;
	
	for (i = length_before; i >= position; i--) {
		if (i + 1 >= len)
			continue;
		
		elt1 = ith(i);
		elt2 = ith(i + 1);
		A->element_move(elt1, elt2, FALSE);
		}
	copy_in(position, elt);
}

void vector_ge::append(INT *elt)
{
	reallocate_and_insert_at(len, elt);
}

void vector_ge::copy_in(INT i, INT *elt)
{
	INT *elt2 = ith(i);
	A->element_move(elt, elt2, FALSE);
};

void vector_ge::copy_out(INT i, INT *elt)
{
	INT *elt2 = ith(i);
	A->element_move(elt2, elt, FALSE);
};

void vector_ge::conjugate_svas(INT *Elt)
{
	INT i;
	INT *Elt1, *Elt2, *Elt3;

	Elt1 = NEW_INT(A->elt_size_in_INT);
	Elt2 = NEW_INT(A->elt_size_in_INT);
	Elt3 = NEW_INT(A->elt_size_in_INT);
	A->invert(Elt, Elt1);
	for (i = 0; i < len; i++) {
		A->element_mult(Elt1, ith(i), Elt2, FALSE);
		A->element_mult(Elt2, Elt, Elt3, FALSE);
		A->element_move(Elt3, ith(i), FALSE);
		}
	FREE_INT(Elt1);
	FREE_INT(Elt2);
	FREE_INT(Elt3);
};

void vector_ge::conjugate_sasv(INT *Elt)
{
	INT i;
	INT *Elt1, *Elt2, *Elt3;

	Elt1 = NEW_INT(A->elt_size_in_INT);
	Elt2 = NEW_INT(A->elt_size_in_INT);
	Elt3 = NEW_INT(A->elt_size_in_INT);
	A->invert(Elt, Elt1);
	for (i = 0; i < len; i++) {
		A->element_mult(Elt, ith(i), Elt2, FALSE);
		A->element_mult(Elt2, Elt1, Elt3, FALSE);
		A->element_move(Elt3, ith(i), FALSE);
		}
	FREE_INT(Elt1);
	FREE_INT(Elt2);
	FREE_INT(Elt3);
};

void vector_ge::print_with_given_action(ostream &ost, action *A2)
{
	INT i, l;

	l = len;
	for (i = 0; i < l; i++) {
		ost << "generator " << i << ":" << endl;
		A->element_print_quick(ith(i), ost);
		ost << endl;
		A2->element_print_as_permutation(ith(i), ost);
		}
}

void vector_ge::print(ostream &ost, INT f_print_as_permutation, 
	INT f_offset, INT offset, INT f_do_it_anyway_even_for_big_degree, 
	INT f_print_cycles_of_length_one)
{
	INT i, l;
	
	l = len;
	if (!f_offset)
		offset = 0;
	ost << "Strong generators: (" << l << " of them)" << endl;
	ost << "f_print_as_permutation=" << f_print_as_permutation << endl;
	for (i = 0; i < l; i++) {
		ost << "generator " << i << ":" << endl;
		A->element_print_quick(ith(i), ost);
		ost << endl;
		if (f_print_as_permutation) {
			//A->element_print_as_permutation(ith(i), ost);
			A->element_print_as_permutation_with_offset(ith(i), ost, 
				offset, f_do_it_anyway_even_for_big_degree, 
				f_print_cycles_of_length_one, 0/*verbose_level*/);
			}
		}
}

void vector_ge::write_to_memory_object(memory_object *m, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;

	if (f_v) {
		cout << "vector_ge::write_to_memory_object" << endl;
		}
	m->write_int(len);
	for (i = 0; i < len; i++) {
		A->element_write_to_memory_object(ith(i), m, 0);
		}
}

void vector_ge::read_from_memory_object(memory_object *m, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, l;

	if (f_v) {
		cout << "vector_ge::read_from_memory_object" << endl;
		}
	m->read_int(&l);
	allocate(l);
	for (i = 0; i < len; i++) {
		A->element_read_from_memory_object(ith(i), m, 0);
		}
}

void vector_ge::write_to_file_binary(ofstream &fp, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;

	if (f_v) {
		cout << "vector_ge::write_to_file_binary" << endl;
		}
	fp.write((char *) &len, sizeof(INT));
	for (i = 0; i < len; i++) {
		A->element_write_to_file_binary(ith(i), fp, 0);
		}
}

void vector_ge::read_from_file_binary(ifstream &fp, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, l;

	if (f_v) {
		cout << "vector_ge::read_from_file_binary" << endl;
		}
	fp.read((char *) &l, sizeof(INT));
	allocate(l);
	for (i = 0; i < len; i++) {
		A->element_read_from_file_binary(ith(i), fp, 0);
		}
}




