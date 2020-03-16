#include <iostream> 

void merror(char* v) {
	std::cout << v << "\n";
	exit(1);
}


template <class T>
class vect {
private:
	int nr;      // number of rows in the matrix 
	T* ets;     // entries of the matrix

	void init(int n)
	{
		nr = n;
		ets = new T[nr];
	};
public:
	//vect(int, double = 0.0);     // construct a matrix
	vect();
	vect(int row)
	{
		nr = row;
		init(nr);
	};
	~vect()   // destructor		
	{
		delete[] ets;
	}

	void zeros(int nr);
	T& operator() (int row);
	T& operator[] (int row);

	inline int size() const;

	//	double& operator() (int row, int col) const;
};


template<class T>
vect<T> ::vect() : nr(0), ets(NULL) {};

template <class T>
void vect<T>::zeros(int nr) {
	init(nr);
	if (!ets) perror("no more space");
	for (int i = 0; i < nr; i++) {
		ets[i] = 0;
	}
}

template <class T>
T& vect<T>::operator() (int row) {
	if (row >= nr + 1)
		merror("vect subscript out of bounds");
	else if (row <= 0)
		merror("vect subscript out of bounds");
	return ets[row - 1];
}

template <class T>
T& vect<T>::operator[] (int row) {
	if (row >= nr + 1)
		merror("vect subscript out of bounds");
	else if (row <= 0)
		merror("vect subscript out of bounds");
	return ets[row - 1];
}


template <class T>
inline int vect<T>::size() const
{
	return nr;
}

//-------------------------------------
template <class T>
class mat2D {
private:
	int nr;      // number of rows in the matrix 
	int nc;      // number of columns in the matrix 
	int ind;
	T* ets;    // entries of the matrix

	void init(int row, int col)
	{
		nr = row;
		nc = col;
		ets = new T[nr*nc];
	};

public:
	mat2D();
	mat2D(int nr, int nc);
	~mat2D()
	{                       // destructor		
		delete[] ets;
	}

	T& operator() (int row, int col);
	void zeros(int nr, int nc);
	inline int nrows() const;
	inline int ncols() const;

};



template<class T>
mat2D<T>::mat2D() : nr(0), nc(0), ets(NULL) {};

template <class T>
mat2D<T>::mat2D(int row, int col) // constructor
{
	nr = row;
	nc = col;
	init(nr, nc);

}

template <class T>
void mat2D<T>::zeros(int row, int col) {
	nr = row;
	nc = col;
	//init(int nr, int nc);
	ets = new T[nr*nc];
	for (int i = 0; i < nr; i++) {
		for (int j = 0; j < nc; j++)
		{
			ets[i*nr + j] = 0;
		}
	}
}

template <class T>
T& mat2D<T>::operator() (int row, int col) {
	if (row >= nr + 1 || col >= nc + 1)
		merror("Matrix subscript out of bounds");
	else if (row <= 0 || col <= 0)
		merror("Matrix subscript out of bounds");

	ind = (row - 1) * (nr)+(col - 1);
	return ets[ind];
}



template <class T>
inline int mat2D<T>::nrows() const
{
	return nr;
}

template <class T>
inline int mat2D<T>::ncols() const
{
	return nc;
}


//-----------------------------------------

template <class T>
class mat3D {

private:
	int nr;      // number of rows in the matrix 
	int nc;      // number of columns in the matrix 
	int ndim;
	int ind;
	T* ets;   // entries of the matrix

	void init(int row, int col, int nd)
	{
		nr = row;
		nc = col;
		ndim = nd;
		ets = new T[nr*nc*ndim];
	};

public:
	mat3D();
	mat3D(int nr, int nc, int ndim); // constructor
	~mat3D() {
		delete[] ets;
	};

	T& operator() (int, int, int);
	inline int dim1() const;
	inline int dim2() const;
	inline int dim3() const;
	void zeros(int nr, int nc, int ndim);

};

template<class T>
mat3D<T>::mat3D() : nr(0), nc(0), ndim(0), ets(NULL) {};

template <class T>
mat3D<T>::mat3D(int row, int col, int nd) // constructor
{
	nr = row;
	nc = col;
	ndim = nd;

	init(nr, nc, ndim);

}

template <class T>
void mat3D<T>::zeros(int row, int col, int nd) {
	nr = row;
	nc = col;
	ndim = nd;
	init(nr, nc, ndim);
	if (!ets) perror("no more space");
	for (int iii = 0; iii < nd; iii++) {
		for (int ii = 0; ii < nc; ii++)
		{
			for (int i = 0; i < nr; i++)
			{
				ets[iii*nr*nc + ii*nr + i] = 0;
			}

		}
	}

}



template <class T>
T& mat3D<T>::operator() (int row, int col, int nd)
{
	if (row >= nr + 1 || col >= nc + 1 || nd >= nd + 1)
		merror("Matrix subscript out of bounds");
	else if (row <= 0 || col <= 0 || nd <= 0)
		merror("Matrix subscript out of bounds");

	ind = (nd - 1)* (nr*nc) + (col - 1)*nr + (row - 1);
	return ets[ind];
}



template <class T>
inline int mat3D<T>::dim1() const
{
	return nr;
}

template <class T>
inline int mat3D<T>::dim2() const
{
	return nc;
}

template <class T>
inline int mat3D<T>::dim3() const
{
	return ndim;
}

//-----------------------------------------------------------



template <class T>
class mat4D {

private:
	int nr;      // number of rows in the matrix 
	int nc;      // number of columns in the matrix 
	int ndim1;
	int ndim2;
	int ind;
	T* ets;   // entries of the matrix

	void init(int nr, int nc, int ndim1, int ndim2) {
		ets = new T[nr*nc*ndim1*ndim2];
	};

public:
	mat4D();
	mat4D(int row, int col, int nd1, int nd2);

	~mat4D() {
		delete[] ets;
	};

	T& operator() (int, int, int, int);
	inline int dim1() const;
	inline int dim2() const;
	inline int dim3() const;
	inline int dim4() const;
	void zeros(int nr, int nc, int ndim1, int ndim2);
};


template<class T>
mat4D<T>::mat4D() : nr(0), nc(0), ndim1(0), ndim2(0), ets(NULL) {};

template <class T>
mat4D<T>::mat4D(int row, int col, int nd1, int nd2) // constructor
{
	nr = row;
	nc = col;
	ndim1 = nd1;
	ndim2 = nd2;
	init(nr, nc, ndim1, ndim2);

}

template <class T>
void mat4D<T>::zeros(int row, int col, int nd1, int nd2)
{
	nr = row;
	nc = col;
	ndim1 = nd1;
	ndim2 = nd2;
	init(nr, nc, nd1, nd2);
	if (!ets) perror("no more space");
	for (int iiii = 0; iiii < ndim2; iiii++) {
		for (int iii = 0; iii < ndim1; iii++) {
			for (int ii = 0; ii < nc; ii++) {
				for (int i = 0; i < nr; i++) {
					ets[iiii*nr*nc*ndim1 + iii * nc*nr + ii * nr + i] = 0;
				}
			}

		}
	}

}

template <class T>
T& mat4D<T>::operator() (int row, int col, int nd1, int nd2)
{
	if (row >= nr + 1 || col >= nc + 1 || nd1 >= nd1 + 1 || nd2 >= nd2 + 1)
		merror("Matrix subscript out of bounds");
	else if (row <= 0 || col <= 0 || nd1 <= 0 || nd1 <= 0)
		merror("Matrix subscript out of bounds");

	ind = (nd2 - 1)* (nr*nc*ndim1) + (nd1 - 1)* (nr*nc) + (col - 1)*nr + (row - 1);
	return ets[ind];
}


template <class T>
inline int mat4D<T>::dim1() const
{
	return nr;
}

template <class T>
inline int mat4D<T>::dim2() const
{
	return nc;
}

template <class T>
inline int mat4D<T>::dim3() const
{
	return nd1;
}

template <class T>
inline int mat4D<T>::dim4() const
{
	return nd2;
}


typedef vect<int> ivec;
typedef vect<float> fvec;
typedef mat2D<float> fmat2d;
typedef mat2D<int> imat2d;
typedef mat3D<int> imat3d;
typedef mat3D<short> smat3d;
typedef mat3D<float> fmat3d;
typedef mat4D<float> fmat4d;
typedef mat4D<short> smat4d;


// Todo
// - inint for mat2d
