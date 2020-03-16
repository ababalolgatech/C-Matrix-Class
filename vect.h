#ifndef VECT_H
#define VECT_H


template <class T>
class vect {
private:
	int start_row;       // start index
	int end_row;       // end index
	int nr;       // number of rows
	mutable int ind;

	

	void init()
	{
		nr = end_row - start_row + 1;
		ets = new T[nr];
	};
public:
	T* ets;     // entries of the matrix
	//vect(int, double = 0.0);     // construct a matrix
	vect();
	vect(int srt_row, int en_row);
	vect(int inr);
	vect(const vect&rhs); // copy constructor
	~vect();

	void zeros(int srt_row, int en_row);
	void zeros(int inr);
	void rand(int inr);
	void rand(int inr, T seed);
	void ones(int inr);
	void print();
	T& operator() (int irow); // indexing is the input here
	T& operator[] (int irow) const;
	vect<T>& operator= (const vect<T>&); // copy assignment 
	vect<T>& operator+= (const vect<T>&); // v +=v2
	vect<T>& operator-= (const vect<T>&); // v +=v2

	inline int size() const;
	inline int begin() const;
	inline int end() const;

};



template<class T>
vect<T>& vect<T> ::operator= (const vect<T>& rhs) {
		/*		
		if vector and rhs were different sizes, vector
		has been resized to match the size of rhs
		*/

	if (this != &rhs)
	{
		if (nr != rhs.nr)
		{
			if (ets != NULL)
			{
				delete ets;
			}
			start_row = rhs.start_row;
			end_row = rhs.end_row;
			init(); // nr and ets is set here			
		}
	}


	for (int i = 0; i < nr; i++)
	{
		    // ind = i - start_row;
		ets[i] = rhs.ets[i]; // [] is overloaded for rhs.. starts at start_row
	}

	return *this;
};

template<class T>
vect<T>& vect<T>:: operator+= (const vect<T>& rhs) {
	if (nr != rhs.nr)
	{
		merror("bad vector sizes");
	}
	for (int i = start_row; i <= nr; i++)
	{
		      ind = i - start_row;
		ets[ind] += rhs[i]       ;  // [] is overloaded for rhs
	}
	return *this;
}

template<class T>
vect<T>& vect<T>:: operator-= (const vect<T>& rhs) {
	if (nr != rhs.nr)
	{
		merror("bad vector sizes");
	}
	for (int i = start_row; i <= nr; i++)
	{
		      ind = i - start_row;
		ets[ind] -= rhs[i];  //  [] is overloaded for rhs
	}
	return *this;
}

template<class T>
vect<T> operator+(const vect<T>& rhs) {// usage: v1 = + v2;
	return rhs;
}

template<class T>
vect<T> operator-(const vect<T>& rhs) { // usage: v1 = - v2;
	return vect<T>(rhs.nr) - rhs;
}

template<class T>
vect<T> operator+(const vect<T>& lhs, const vect<T>& rhs) {
	if (lhs.nr != rhs.nr)
	{
		merror("bad vector sizes");
	}
	vect<T> sum = lhs; // It would cause problem without copy constructor
	sum += rhs;
	return sum;
}

template<class T>
vect<T> operator-(const vect<T>& lhs, const vect<T>& rhs) {
	if (lhs.nr != rhs.nr)
	{
		merror("bad vtor sizes");
	}
	vect<T> sum = lhs; // It would cause problem without copy constructor
	sum -= rhs;
	return sum;
}

template <class T>
vect<T>operator*(const double scalar, const vect<T>& rhs) {
	
	vect<T> tm(rhs.size());


	for (int i = rhs.begin(); i <= rhs.size(); i++)
	{
		tm[i] = rhs[i] * scalar;  // note, it is overloaded []..default indexing at 1
	}

	// first object (ets) is destroyed when you return
	return tm;
}

template <class T>
vect<T>operator*(const vect<T>& v, const double scalar) {

	return scalar * v ;
}

template<class T>
vect<T> operator*(const vect<T>& lhs, const vect<T>& rhs) {
	int ind;
	int sz = lhs.nr;
	if (sz != rhs.nr) {
		merror("bad vtor sizes");
	}
	vect<T> tm(sz);
	for (int i = rhs.start_row; i <= sz; i++)
	{
		   // ind = i - rhs.start_row  ;  // [] is overloaded for vectors
		tm[i] = lhs[i] * rhs[i];
	}
	return tm;
}


template <class T>
vect<T>operator/(const vect<T>& rhs, const double scalar) {
	if (scalar == 0)
	{
		merror("division by zero in vector-scalar division");
	}
	return (1.0 / scalar)*rhs;
	
}


template <class T>
vect<T>operator/(const double scalar, const vect<T>& rhs) {
	if (scalar == 0)
	{
		merror("division by zero in vector-scalar division");
	}

	vect<T> tm(rhs.size());

	for (int i = rhs.begin(); i <= rhs.size(); i++)
	{
		tm[i] =   scalar/ rhs[i];  // note, it is overloaded []..default indexing at 1
	}

	return tm;
}

template <class T>
vect<T>operator+(const vect<T>& rhs, const T scalar) {
	vect<T> tm(rhs.size());
	for (int i = rhs.begin(); i <= rhs.size(); i++)
	{
		tm[i] = scalar + rhs[i];  // note, it is overloaded []..default indexing at 1
	}

	return tm;
}


template <class T>
vect<T>operator+(const double scalar,const vect<T>& rhs) {
	vect<T> tm(rhs.size());
	for (int i = rhs.begin(); i <= rhs.size(); i++)
	{
		tm[i] = scalar + rhs[i];  // note, it is overloaded []..default indexing at 1
	}

	return tm;
}

template <class T>
vect<T>operator-(const vect<T>& rhs, const double scalar) {
	vect<T> tm(rhs.size());
	for (int i = rhs.begin(); i <= rhs.size(); i++)
	{
		tm[i] = rhs[i] - scalar;  // note, it is overloaded []..default indexing at 1
	}

	return tm;
}


template <class T>
vect<T>operator-(const double scalar, const vect<T>& rhs) {
	vect<T> tm(rhs.size());
	for (int i = rhs.begin(); i <= rhs.size(); i++)
	{
		tm[i] = scalar - rhs[i];  // note, it is overloaded []..default indexing at 1
	}

	return tm;

}


template<class T>
vect<T> ::vect() : nr(0), start_row(0),end_row(0),ind(0), ets(NULL) {};

template<class T> // copy constructor
vect<T>::vect(const vect& rhs) { // you must account for all member variables
	start_row = rhs.start_row;
	  end_row = rhs.end_row;
	      ind = rhs.ind;

//	init(); // nr & ets is set here

	nr = rhs.nr;
	ets =  new T [nr];

	for (int i = 0; i < nr; i++)
	{
		 
		ets[i] = rhs.ets[i];  // note, it is overloaded []..default indexing at 1
	}
}

template<class T>
vect<T>::vect(int inr){

	start_row = 1  ;
	  end_row = inr;
	  init();
};

template<class T>
vect<T>::vect(int srt_row, int en_row)
{

	start_row = strt_row;
	end_row = en_row;
	init();
};

template <class T>
void vect<T>::zeros(int strt_row, int en_row) {
	start_row = strt_row;
	end_row = en_row;

	init();

	for (int i = 0; i < nr; i++)
	{
		ets[i] = 0;
	}
}

template <class T>
void vect<T>::ones(int inr) {
	start_row = 1;
	end_row = inr;

	init();

	for (int i = 0; i < nr; i++)
	{
		ets[i] = 1;
	}
}

template <class T>
void vect<T>::zeros(int inr) {
	start_row = 1;
	end_row = inr;

	init();

	for (int i = 0; i < nr; i++)
	{
		ets[i] = 0;
	}
}

template <class T>
void vect<T>::print() {
	int np;
	if (nr > 100) { np = 100; }
	else
	{ np = nr; }

	for (int i = 0; i < np; i++)
	{
		std::cout << ets[i] << std::endl;
	}

}

template <class T>
T& vect<T>::operator() (int irow) {

	assert(irow >= start_row);
	assert(irow <= end_row);


	ind = irow - start_row;
	return ets[ind];
}

template <class T>
T& vect<T>::operator[] (int irow) const {

	assert(irow >= start_row);
	assert(irow <= end_row);


	ind = irow - start_row;
	return ets[ind];
}


template <class T>
inline int vect<T>::size() const
{
	return nr;
}

template <class T>
inline int vect<T>::begin() const
{
	return start_row;
}

template <class T>
inline int vect<T>::end() const
{
	return end_row;
}

template <class T>
vect<T>::~vect() {  // destructor
	if (ets != NULL) delete[](ets);
	ets = NULL; // why the second Null ?
}



template <class T>
void vect<T> ::rand(int inr, T seed) {
	T p;
	vect<T> ixv;
	start_row = 1     ;
	  end_row = inr   ;	
	   
	
	init();
	ixv.zeros(13);
	ixv(ixv.begin()) = seed;
	for (int i = 0; i < 1000; i++)
	{
		p = acorni(ixv);
	}

	for (int i = 0; i < nr; i++)
	{
		ets[i] = acorni(ixv);
	}
}

template <class T>
void vect<T> ::rand(int inr) {
	T p, seed;
	vect<T> ixv;
	start_row = 1;
	end_row = inr;
	seed = 100000;

	init();
	ixv.zeros(13);
	ixv(ixv.begin()) = seed;
	for (int i = 0; i < 1000; i++)
	{
		p = acorni(ixv);
	}

	for (int i = 0; i < nr; i++)
	{
		ets[i] = acorni(ixv);
	}
}

/*
ISSUES RESOLVED
----------------
1. Link2019 error due to templated friend member function
https://stackoverflow.com/questions/29055696/why-lnk1120-lnk2019-appears-in-case-of-template-and-friend-function

GTSL Bounds check
random number generator
*/


#endif