
#ifndef _MAT2D_H
#define _MAT2D_H

template <class T>
class mat2D {
private:
	int start_row;         // start index
	int end_row;           // end index
	int start_col;         // start index
	int end_col;           // end index

	int nr;      // number of rows in the matrix 
	int nc;      // number of columns in the matrix 
	mutable	int ind;
	

	void init()	{
		nr = end_row - start_row + 1;
		nc = end_col - start_col + 1;
		ets = new T[nr*nc];
	};

public:
	T* ets;    // entries of the matrix - should be private
	mat2D();
	mat2D(int srt_row, int en_row, int srt_col, int en_col);
	mat2D(int inr, int inc);
	//template<class U>
	mat2D(const mat2D<T>&);  // copy constructor
	~mat2D(){   delete[] ets;	} // destructor

	mat2D<T>& operator=(const mat2D<T>& ); 
	mat2D<T>& operator+=(const mat2D<T>& );
	mat2D<T>& operator-=(const mat2D<T>& at);
	mat2D<T>& operator+();
	mat2D<T>& operator+(const mat2D<T>& mat);
	vect<T> operator*(const vect<T>& v) const; // matrix vector multiply
	template<class U>
	friend mat2D<T> operator-(const mat2D<U>& mat); // unary -, m1 = -m2
	template<class U>
	friend mat2D<T> operator-(const mat2D<U>& mat1, const mat2D<U>& mat2); // binary -
	T& operator() (int irow, int icol) const;
	void zeros(int srt_row, int en_row, int srt_col, int en_col);
	void zeros(int inr, int inc);
	void ones(int srt_row, int en_row, int srt_col, int en_col);
	void ones(int inr, int inc);
	
	void print();
	inline int nrows() const;
	inline int ncols() const;
	inline int row_begin() const;
	inline int row_end() const;
	inline int col_begin() const;
	inline int col_end() const;
	vect<T> col(int);
	vect<T> row(int);
};

template <class T>
vect<T> mat2D<T>::col(int column) {
	assert(column >= start_col);
	assert(column <= end_col);
	
	/*{
		merror(" column is out of bounds");
	}*/

	column -= 1;
	
	vect<T> tm;
	tm.zeros(nc);

	for (int i = 0; i < nc; i++)
	{
		tm.ets[i] += ets[i*nr + column];
	}

	return tm;
};

template <class T>
vect<T> mat2D<T>::row(int row) {
	assert(row >= start_row);
	assert(row <= end_row);

	row = row - start_col;
	vect<T> tm;
	tm.zeros(nc);

	for (int i = 0; i < nr; i++)
	{
		tm.ets[i] = ets[row*nc + i];
	}

	return tm;
};

template <class T>
mat2D<T>::mat2D(const mat2D<T>& mat):start_row(mat.start_row),start_col(mat.start_col),end_row(mat.end_row),end_col(mat.end_col)
{ // copy constructor
	if (this != &mat) 
	{
// if (nr != mat.nr || nc != mat.nc){	merror("bad matrix sizes");		}
		init();
		
		for (int i = 0; i < nr; i++)
		{
			for (int ii = 0; ii < nc; ii++)
			{
				ets[ii*nr + i] = mat.ets[ii*nr + i];
			}
		}
	}
	//return *this;
}
template <class T>
mat2D<T>& mat2D<T>::operator=(const mat2D<T>& mat) {//copy assignment
	if (this != &mat) 
	{
		//int i, j, nel;
		if (nr != mat.nr || nc != mat.nc )
		{
			if (ets != NULL)
			{
				delete ets;
			}

			start_row = mat.start_row;
			start_col = mat.start_col;
			  end_row = mat.end_row;
			  end_col = mat.end_col;

			init();	 // nr and nc defined here

		}
		for (int ii = 0; ii < nc; ii++) 
		{
			for (int i = 0; i < nr; i++)
			{
				ets[ii*nr + i] = mat.ets[ii*nr + i]; 
			}
		}
		
	}
	return *this;
}

template <class T>
mat2D<T>& mat2D<T>::operator+=(const mat2D<T>& mat) {
	if (this != mat)
	{
		if (nrows != mat.nrows || ncols != mat.ncols)
		{
			merror("bad matrix sizes");
		}
		for (int ii = 0; ii < nc; ii++)
		{
			for (int i = 0; i < nr; i++)
			{
				ets[ii*nr + i] += mat.ets[ii*nr + i];
			}
		}

	}
	return *this;
}

template<class T>
mat2D<T>& mat2D<T>::operator-=(const mat2D<T>& mat) {
	if (this != mat)
	{
		if (nrows != mat.nrows || ncols != mat.ncols)
		{
			merror("bad matrix sizes");
		}
		for (int ii = 0; ii < nc; ii++)
		{
			for (int i = 0; i < nr; i++)
			{
				ets[ii*nr + i] -= mat.ets[ii*nr + i];
			}
		}

	}
	return *this;
}

template <class T>
mat2D<T>& mat2D<T>::operator+() {// usage: mat1 = + mat2;
	return *this;
}

template <class T>
mat2D<T>& mat2D<T>::operator+(const mat2D<T>& mat) {
	mat2D<T> sum = *this;
	sum += mat;
	return sum;
}

template <class T>
vect<T> mat2D<T>::operator*(const vect<T>& rhs)  const { // mat-vec mult
	int row, col;
	if (nc != rhs.size())
	{
		merror("matrix and vector sizes do not match");
	}
	vect<T> tm;
	tm.zeros(nr);

	for (int ii = 0; ii < nc; ii++)
	{
		for (int i = 0; i < nr; i++)
		{	
		  	tm.ets[ii] += ets[ii*nr + i] * rhs.ets[i]; 
		}
	}
	return tm;
}

template<class T>
mat2D<T> operator-(const mat2D<T>& mat) {
	// usage: mat1 = - mat2;
	return mat2D<T>(mat.nrows, mat.ncols) - mat;
}

template<class T>
mat2D<T> operator-(const mat2D<T>& mat1, const mat2D<T>& mat2) {
	if (mat1.nrows != mat2.nrows || mat1.ncols != mat2.ncols)
		error("bad matrix sizes");
	mat2D<T> sum = mat1;
	sum -= mat2;
	return sum;
}

template<class T>
mat2D<T>::mat2D() : nr(0), nc(0), ets(NULL) {};

template <class T>
mat2D<T>::mat2D(int srt_row, int en_row, int srt_col, int en_col) // constructor
{
	start_row = srt_row;
	end_row = en_row;
	start_col = srt_col;
	end_col = en_col;
	init();

}

template <class T>
mat2D<T>::mat2D(int inr, int inc) // constructor
{
	start_row = 1;
	end_row = inr;
	start_col = 1;
	end_col = inc;
	init();

}

template <class T>
void mat2D<T>::zeros(int srt_row, int en_row, int srt_col, int en_col) {
	start_row = srt_row;
	end_row = en_row;
	start_col = srt_col;
	end_col = en_col;
	init();
	if (!ets) perror("no more space");
	for (int ii = 0; ii < nc; ii++) {
		for (int i = 0; i < nr; i++)
		{
			ets[ii*nr + i] = 0;
		}
	}
}

template <class T>
void mat2D<T>::ones(int srt_row, int en_row, int srt_col, int en_col) {
	start_row = srt_row;
	end_row = en_row;
	start_col = srt_col;
	end_col = en_col;
	init();
	if (!ets) perror("no more space");
	for (int ii = 0; ii < nc; ii++) {
		for (int i = 0; i < nr; i++)
		{
			ets[ii*nr + i] = 0;
		}
	}
}

template <class T>
void mat2D<T>::zeros(int inr, int inc) {
	start_row = 1;
	end_row = inr;
	start_col = 1;
	end_col = inc;
	init();
	if (!ets) perror("no more space");
	for (int ii = 0; ii < nc; ii++) {
		for (int i = 0; i < nr; i++)
		{
			ets[ii*nr + i] = 0;
		}
	}
}

template <class T>
void mat2D<T>::ones(int inr, int inc) {
	start_row = 1;
	end_row = inr;
	start_col = 1;
	end_col = inc;
	init();
	if (!ets) perror("no more space");
	for (int ii = 0; ii < nc; ii++) {
		for (int i = 0; i < nr; i++)
		{
			ets[ii*nr + i] = 1;
		}
	}
}

template <class T>
void mat2D<T>::print() {
	int nr_, nc_;


	if (nr > 100) { nr_ = 100; }
	else { nr_ = nr; }

	if (nc > 100) { nc_ = 100; }
	else { nc_ = nc; }

	for (int i = 0; i < nr_; i++)
	{
		for (int ii = 0; ii < nc_; ii++)
		{
			ind = (i*nc) + ii;
			std::cout << ets[ind] << ' ';   // print from colums first 
		}
		std::cout << std::endl;
	}

}

template <class T>
T& mat2D<T>::operator() (int irow, int icol)const {
	int row, col;
	assert(irow >= start_row &&  irow <= end_row);
	assert(icol >= start_col && icol <= end_col);

	row = irow - start_row;
	col = icol - start_col;

	ind = (row *nc) + col;
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

template <class T>
inline int mat2D<T>::row_begin() const
{
	return start_row;
}

template <class T>
inline int mat2D<T>::row_end() const
{
	return end_row;
}


template <class T>
inline int mat2D<T>::col_begin() const
{
	return start_col;
}

template <class T>
inline int mat2D<T>::col_end() const
{
	return end_col;
}

#endif