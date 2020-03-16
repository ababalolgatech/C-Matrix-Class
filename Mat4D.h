
#ifndef _MAT4D_H_
#define _MAT4D_H_

template <class T>
class mat4D {

private:
	int start_row;         // start index
	int end_row;         // end index
	int start_col;         // start index
	int end_col;         // end index
	int start_dim1;        // start index
	int end_dim1;        // end index
	int start_dim2;        // start index
	int end_dim2;        // end index

	int nr;      // number of rows in the matrix 
	int nc;      // number of columns in the matrix 
	int ndim1;
	int ndim2;
	mutable	int ind;
	T* ets;   // entries of the matrix

	void init() {
		nr = end_row - start_row + 1;
		nc = end_col - start_col + 1;
		ndim1 = end_dim1 - start_dim1 + 1;
		ndim2 = end_dim2 - start_dim2 + 1;

		ets = new T[nr*nc*ndim1*ndim2];  // ndim2 is set in the constructor definition
	};

public:
	mat4D();
	mat4D(const mat4D<T>& mat);  // copy constructor
	mat4D<T>& operator=(const mat4D<T>&); // copy assignment
	mat4D(int srt_row, int en_row, int srt_col, int en_col, int srt_dim1, int en_dim1, int srt_dim2, int en_dim2);
	mat4D(int inr, int inc, int idim1, int idim2);
	// note dim3 is specified for the covariance matrix
	~mat4D() { delete[] ets; };

	void zeros(int srt_row, int en_row, int srt_col, int en_col, int srt_dim1, int en_dim1, int srt_dim2, int en_dim2);
	void zeros(int inr, int inc, int idim1, int idim2);
	T& operator() (int irow, int icol, int idim1, int idim2);
	void  print();
	inline int dim1() const;
	inline int dim2() const;
	inline int dim3() const;
	inline int dim4() const;
	vect<T> col(int dim2, int dim1, int row);
	vect<T> row(int dim2, int dim1, int col);
	mat2D<T> extract2d(int dim1, int dim2);
	mat3D<T> extract3d(int dim2);
};

template <class T>
vect<T> mat4D<T>::col(int dim2, int dim1, int col) {

	assert(col >= start_col && col <= end_col);
	assert(dim1 >= start_dim1 && dim1 <= end_dim1);
	assert(dim2 >= start_dim2 && dim2 <= end_dim2);

	// Don't modify when using pre-built accessors
	// col -= start_col;
	//dim1 -=  start_dim1;
	//dim2 -=  start_dim2;

	mat3D<T>tmp3D;
	vect<T>tmp1D;

	tmp3D = this->extract3d(dim2);
	tmp1D = tmp3D.col(dim1, col);

	return tmp1D;
}

template <class T>
vect<T> mat4D<T>::row(int dim2, int dim1, int row) {

	assert(row >= start_row && row <= end_row);
	assert(dim1 >= start_dim1 && dim1 <= end_dim1);
	assert(dim2 >= start_dim2 && dim2 <= end_dim2);

	// Leave cos of pre-built accessors
	// row -= start_row;
	//dim1 -= start_dim1; 
	//dim2 -= start_dim2; 


	mat3D<T>tmp3D;
	vect<T>tmp1D;

	tmp3D = this->extract3d(dim2);
	tmp1D = tmp3D.row(dim1, row);

	return tmp1D;
}

template<class T>
mat2D<T> mat4D<T>::extract2d(int dim2, int dim1) {

	assert(dim1 >= start_dim1 && dim1 <= end_dim1);
	assert(dim2 >= start_dim2 && dim2 <= end_dim2);

	// Don't modify when using pre-built accessors
	//dim1 -=  start_dim1;  
	//dim2 -=  start_dim2;

	mat2D<T> tmp2D;
	mat3D<T>tmp3D;
	tmp3D = this->extract3d(dim2);
	tmp3D.print();
	tmp2D = tmp3D.extract2d(dim1);
	return tmp2D;
}

template<class T>
mat3D<T> mat4D<T>::extract3d(int dim2) {

	assert(dim2 >= start_dim2 && dim2 <= end_dim2);
	mat3D<T> cube;
	cube.zeros(nr, nc, ndim1);

	dim2 -= start_dim2;

	for (int iii = 0; iii < ndim1; iii++) {
		for (int ii = 0; ii < nc; ii++) {
			for (int i = 0; i < nr; i++) {
				cube.ets[iii*nr*nc + ii * nr + i] = ets[dim2*nr*nc*ndim1 + iii * nr*nc + ii * nr + i];
			}
		}
	}

	return cube;
}

template<class T>
mat4D<T>::mat4D(const mat4D<T>& mat) { // copy constructor

	start_row = mat.start_row;
	start_col = mat.start_col;
	end_row = mat.end_row;
	end_col = mat.end_col;
	start_dim1 = mat.start_dim1;
	end_dim1 = mat.end_dim1;
	start_dim2 = mat.start_dim2;
	end_dim2 = mat.end_dim2;

	if (this != &mat)
	{
		init();

		for (int iiii = 0; iiii < ndim2; iiii++) {
			for (int iii = 0; iii < ndim1; iii++) {
				for (int ii = 0; ii < nc; ii++) {
					for (int i = 0; i < nr; i++)
					{
						ets[iiii*nr*nc*ndim1 + iii * nr*nc + ii * nr + i] =
							mat.ets[iiii*nr*nc*ndim1 + iii * nr*nc + ii * nr + i];  // edit
					}
				}

			}
		}

	}



}

template<class T>
mat4D<T>& mat4D<T>::operator=(const mat4D<T>& mat) {
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
			start_dim1 = mat.start_dim1;
			end_dim1 = mat.end_dim1;
			start_dim2 = mat.start_dim2;
			end_dim2 = mat.end_dim2;

			init();	 // nr and nc defined here
			for (int iii = 0; iii < ndim; iii++) {
				for (int ii = 0; ii < nc; ii++) {
					for (int i = 0; i < nr; i++) {
						ets[iiii*nr*nc*ndim1 + iii * nr*nc + ii * nr + i] =
							mat.ets[iiii*nr*nc*ndim1 + iii * nr*nc + ii * nr + i];  // edit
					}

				}
			}

		}
	}
	return *this;
}

template<class T>
mat4D<T>::mat4D() : nr(0), nc(0), ndim1(0), ndim2(0), ets(NULL) {};

template <class T>
mat4D<T>::mat4D(int srt_row, int en_row, int srt_col, int en_col, int srt_dim1, int en_dim1, int srt_dim2, int en_dim2) // constructor
{
	start_row = srt_row;
	end_row = en_row;
	start_col = srt_col;
	end_col = en_col;
	start_dim1 = srt_dim1;
	end_dim1 = en_dim1;
	start_dim2 = srt_dim2;
	end_dim2 = en_dim2;

	init();

}

template <class T>
mat4D<T>::mat4D(int inr, int inc, int ndim1, int ndim2) { // constructor

	start_row = 1;
	end_row = inr;
	start_col = 1;
	end_col = inc;
	start_dim1 = 1;
	end_dim1 = ndim1;
	start_dim2 = 1;
	end_dim2 = ndim2;


	init();
}

template <class T>
void mat4D<T>::zeros(int srt_row, int en_row, int srt_col, int en_col, int srt_dim1, int en_dim1, int srt_dim2, int en_dim2) {
	start_row = srt_row;
	end_row = en_row;
	start_col = srt_col;
	end_col = en_col;
	start_dim1 = srt_dim1;
	end_dim1 = en_dim1;
	start_dim2 = srt_dim2;
	end_dim2 = en_dim2;




	init();
	if (!ets) perror("no more space");
	for (int iiii = 0; iiii < ndim2; iiii++) {
		for (int iii = 0; iii < ndim1; iii++) {
			for (int ii = 0; ii < nc; ii++) {
				for (int i = 0; i < nr; i++) {
					ets[iiii*nr*nc*ndim1 + iii * nr*nc + ii * nr + i] = 0;  // edit
				}
			}

		}
	}

}

template <class T>
void mat4D<T>::zeros(int inr, int inc, int idim1, int idim2)
{
	start_row = 1;
	end_row = inr;
	start_col = 1;
	end_col = inc;
	start_dim1 = 1;
	end_dim1 = idim1;
	start_dim2 = 1;
	end_dim2 = idim2; // note 


//		merror("Matrix subscript of fourth dimension is out of bounds");

	init();
	if (!ets) perror("no more space");
	for (int iiii = 0; iiii < ndim2; iiii++) {
		for (int iii = 0; iii < ndim1; iii++) {
			for (int ii = 0; ii < nc; ii++) {
				for (int i = 0; i < nr; i++) {
					ets[iiii*nr*nc*ndim1 + iii * nr*nc + ii * nr + i] = 0;  // edit
				}
			}

		}
	}

}

template <class T>
T& mat4D<T>::operator() (int irow, int icol, int idim1, int idim2)
{
	int row, col, dim1, dim2;

	assert(irow >= start_row && irow <= end_row);
	assert(icol >= start_col && icol <= end_col);
	assert(idim1 >= start_dim1 && idim1 <= end_dim1);
	assert(idim2 >= start_dim2 && idim2 <= end_dim2);

	row = irow - start_row;
	col = icol - start_col;
	dim1 = idim1 - start_dim1;
	dim2 = idim2 - start_dim2;  // dim2 starts from 0

	ind = (nr*nc*ndim1*dim2) + (nr*nc*dim1) + (row*nc) + col;
	return ets[ind];
}

template <class T>
void mat4D<T>::print() {
	int nr_, nc_, ndim1_, ndim2_;


	if (nr > 100) { nr_ = 100; }
	else { nr_ = nr; }

	if (nc > 100) { nc_ = 100; }
	else { nc_ = nc; }

	if (ndim1 > 5) { ndim1_ = 5; }
	else { ndim1_ = ndim1; }

	if (ndim2 > 5) { ndim2_ = 5; }
	else { ndim2_ = ndim2; }


	for (int iiii = 0; iiii < ndim2_; iiii++)
	{
		std::cout << " 4th dimension = " << iiii + 1 << std::endl;
		std::cout << " -------------- " << std::endl;

		for (int iii = 0; iii < ndim1_; iii++)
		{
			std::cout << " 3rd dimension = " << iii + 1 << std::endl;

			for (int i = 0; i < nr_; i++)
			{
				for (int ii = 0; ii < nc_; ii++)   // print from columns
				{
					ind = (iiii*ndim1_*nr*nc) + (iii*nr*nc) + (i*nc) + ii;
					std::cout << ets[ind] << ' ';
				}
				std::cout << std::endl;
			}


		}
	}

}

template <class T>
inline int mat4D<T>::dim1() const
{
	return nr;
}

template <class T>
inline int mat4D<T>::dim2() const {
	return nc;
}

template <class T>
inline int mat4D<T>::dim3() const {
	return ndim1;
}

template <class T>
inline int mat4D<T>::dim4() const {
	return ndim2;
}



#endif



/*
	To do
	------
	1. copy constructor
	2. copy assigment
*/