
#ifndef _MAT3D_H_
#define _MAT3D_H_


template <class T>
class mat3D {

private:
	int start_row;        // start index
	int end_row;          // end index
	int start_col;        // start index
	int end_col;          // end index
	int start_dim;        // start index
	int end_dim;          // end index


	int nr;      // number of rows in the matrix 
	int nc;      // number of columns in the matrix 
	int ndim;
	mutable	int ind;


	void init()	{
		  nr = end_row - start_row + 1;
		  nc = end_col - start_col + 1;
		ndim = end_dim - start_dim + 1;
		 ets = new T[nr*nc*ndim];
	};

public:
	T* ets;   // entries of the matrix
	mat3D();
	mat3D(const mat3D<T>&);  // copy constructor
	mat3D(int srt_row, int en_row, int srt_col, int en_col, int srt_dim, int en_dim); // constructor
	mat3D(int inr, int inc, int idim1); // constructor
	~mat3D() {
		delete[] ets;
	};

	mat3D<T>& operator=(const mat3D<T>&);
	T& operator() (int irow, int icol, int idim);
	void zeros(int srt_row, int en_row, int srt_col, int en_col, int srt_dim, int en_dim);
	void zeros(int inr, int inc, int idim1);
	void print();
	inline int dim1() const;
	inline int dim2() const;
	inline int dim3() const;
	vect<T> col(int dim, int col);
	vect<T> row(int dim, int row);
	mat2D<T> extract2d(int dim);

};

template<class T>
mat3D<T>& mat3D<T>::operator=(const mat3D<T>& mat) {
	if (this != &mat)
	{
		if (nr != mat.nr || nc != mat.nc)
		{
			if (ets != NULL)
			{
				delete ets;
			}
			start_row = mat.start_row;
			start_col = mat.start_col;
			end_row = mat.end_row;
			end_col = mat.end_col;
	  	  start_dim = mat.start_dim;
			end_dim = mat.end_dim;
			init();	 // nr and nc defined here
			for (int iii = 0; iii < ndim; iii++) {
				for (int ii = 0; ii < nc; ii++)	{
					for (int i = 0; i < nr; i++){
		ets[iii*nr*nc + ii * nr + i] = mat.ets[iii*nr*nc + ii * nr + i];
					}

				}
			}

		}
	}
	return *this;
}

template<class T>
mat3D<T>::mat3D(const mat3D<T>& mat) { // copy constructor
	
	start_row = mat.start_row;
	start_col = mat.start_col;
	end_row = mat.end_row;
	end_col = mat.end_col;
	start_dim = mat.start_dim;
	end_dim = mat.end_dim;

	if (this != &mat) {
		init();

		for (int iii = 0; iii < ndim; iii++) {
			for (int ii = 0; ii < nc; ii++)
			{
				for (int i = 0; i < nr; i++)
				{
		ets[iii*nr*nc + ii*nr + i] = mat.ets[iii*nr*nc + ii * nr + i];
				}

			}
		}

	}

	

}

template <class T>
vect<T> mat3D<T>::col(int dim, int col) {

	assert(col >= start_col);
	assert(col <= end_col);
	assert(dim >= start_dim);
	assert(dim <= end_dim);

	dim -= start_dim;
	col -= start_col;
	vect<T> vec ;
	vec.zeros(nc);

		for (int i = 0; i < nc; i++)
		{
			vec.ets[i] = ets[dim*nr*nc + i*nr + col];
		}

	return vec;
};

template <class T>
vect<T> mat3D<T>::row(int dim, int row) {

	assert(row >= start_row);
	assert(row <= end_row);
	assert(dim >= start_dim);
	assert(dim <= end_dim);

	dim -= start_dim;
	row -= start_row;
	vect<T> vec;
	vec.zeros(nr);

	for (int i = 0; i < nr; i++)
	{
	vec.ets[i] = ets[dim*nr*nc + row*nc + i];
	}

	return vec;
};

template <class T>
mat2D<T> mat3D<T>::extract2d(int dim) {
	
	assert(dim >= start_dim);
	assert(dim <= end_dim);

	dim = dim - start_dim;
	mat2D<T> mat;
	mat.zeros(nr,nc);
	
		for (int ii = 0; ii < nc; ii++)
		{
			for (int i = 0; i < nr; i++)
			{
	   mat.ets[ii*nr+i] = ets[dim*nr*nc + ii*nr+i] ;
			}

		}
	
		return mat;
};

template<class T>
mat3D<T>::mat3D() : nr(0), nc(0), ndim(0), ets(NULL) {};

template <class T>
mat3D<T>::mat3D(int srt_row, int en_row, int srt_col, int en_col, int srt_dim, int en_dim) // constructor
{
	start_row = srt_row;
	end_row = en_row;
	start_col = srt_col;
	end_col = en_col;
	start_dim = srt_dim;
	end_dim = en_dim;

	init();

}

template <class T>
mat3D<T>::mat3D(int inr, int inc, int idim) // constructor
{
	start_row = 1;
	end_row = inr;
	start_col = 1;
	end_col = inc;
	start_dim = 1;
	end_dim = idim;

	init();

}

template <class T>
void mat3D<T>::zeros(int srt_row, int en_row, int srt_col, int en_col, int srt_dim, int en_dim) {

	start_row = srt_row;
	end_row = en_row;
	start_col = srt_col;
	end_col = en_col;
	start_dim = srt_dim;
	end_dim = en_dim;

	init();
	if (!ets) perror("no more space");
	for (int iii = 0; iii < ndim; iii++) {
		for (int ii = 0; ii < nc; ii++)
		{
			for (int i = 0; i < nr; i++)
			{
				ets[iii*nr*nc + ii * nr + i] = 0;
			}

		}
	}

}

template <class T>
void mat3D<T>::zeros(int inr, int inc, int idim1) {

	start_row = 1;
	end_row = inr;
	start_col = 1;
	end_col = inc;
	start_dim = 1;
	end_dim = idim1;

	init();
	if (!ets) perror("no more space");
	for (int iii = 0; iii < ndim; iii++) {
		for (int ii = 0; ii < nc; ii++)
		{
			for (int i = 0; i < nr; i++)
			{
				ets[iii*nr*nc + ii * nr + i] = 0;
			}

		}
	}

}

template <class T>
void mat3D<T>::print() {
	int nr_, nc_, ndim_;


	if (nr > 100) { nr_ = 100; }
	else { nr_ = nr; }

	if (nc > 100) { nc_ = 100; }
	else { nc_ = nc; }

	if (ndim > 10) { ndim_ = 10; }
	else { ndim_ = ndim; }

	for (int iii = 0; iii < ndim_; iii++)
	{
		std::cout << " 3rd dimension = " << iii + 1 << std::endl;

		for (int i = 0; i < nr_; i++)
		{
			for (int ii = 0; ii < nc_; ii++)
			{
				ind = (iii*nr*nc) + (i*nc) + ii;
				std::cout << ets[ind] << ' ';
			}
			std::cout << std::endl;
		}


	}

}

template <class T>
T& mat3D<T>::operator() (int irow, int icol, int idim)
{
	int row, col, dim;
	assert(irow >= start_row && irow <= end_row);
	assert(icol >= start_col && icol <= end_col);
	assert(idim >= start_dim && idim <= end_dim);


	row = irow - start_row;
	col = icol - start_col;
	dim = idim - start_dim;

	ind = (nr*nc*dim) + (row *nc) + col;
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


#endif