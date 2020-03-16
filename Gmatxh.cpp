#include <iostream>

void error(char* v) {
	std::cout << v << "\n";
	exit(1);
}
class mtx2d {
private:
	int nr;      // number of rows in the matrix 
	int ncs;      // number of columns in the matrix 
	int ind;
	double* ets;   // entries of the matrix
public:
	mtx2d(int, int, double = 0.0);     // construct a matrix
	~mtx2d()
	{                       // destructor		
		delete[] ets;
	}

	double& operator() (int row, int col);
	//	double& operator() (int row, int col) const;
};


mtx2d::mtx2d(int n, int m, double a) {
	nr = n;
	ncs = m;
	ets = new double[nr*ncs];
	if (!ets) perror("no more space");
	for (int i = 0; i < nr; i++) {
		for (int j = 0; j < ncs; j++)
		{
			ets[i*nr + j] = 0;
		}
	}
}

double& mtx2d::operator() (int row, int col) {
	if (row >= nr || col >= ncs)
		error("Matrix subscript out of bounds");
	ind = (row) * (nr)+(col);
	return ets[ind];
}



class mat2d {
private:
	int nr;      // number of rows in the matrix 
	int nc;      // number of columns in the matrix 
	int ind;
	double* ets;   // entries of the matrix
public:
	mat2d(int, int, double = 0.0);     // construct a matrix
	~mat2d()
	{                       // destructor		
		delete[] ets;
	}

	double& operator() (int row, int col);

};


mat2d::mat2d(int n, int m, double a) {
	nr = n;
	nc = m;
	ets = new double[nr*nc];
	if (!ets) perror("no more space");
	for (int i = 0; i < nr; i++) {
		for (int j = 0; j < nc; j++)
		{
			ets[i*nr + j] = ((i* nr) + j);
		}
	}
}

double& mat2d::operator() (int row, int col) {
	if (row >= nr + 1 || col >= nc + 1)
		error("Matrix subscript out of bounds");
	else if (row <= 0 || col <= 0)
		error("Matrix subscript out of bounds");

	ind = (row - 1) * (nr)+(col - 1);
	return ets[ind];
}



class mat3d {

private:
	int nr;      // number of rows in the matrix 
	int nc;      // number of columns in the matrix 
	int ndim;
	int ind;
	double* ets;   // entries of the matrix
public:
	mat3d(int, int, int); // constructor
	~mat3d() {
		delete[] ets;
	};

	double& operator() (int, int, int);

};


mat3d::mat3d(int row, int col, int nd) {
	nr = row;
	nc = col;
	ndim = nd;
	ets = new double[nr* nc * nd];
	if (!ets) perror("no more space");
	for (int iii = 0; iii < nd; iii++) {
		for (int ii = 0; ii < nc; ii++)
		{
			for (int i = 0; i < nr; i++)
			{
				ets[iii*nr*nc + ii * nc + i] = 0;
			}

		}
	}

}

double& mat3d::operator() (int row, int col, int nd)
{
	if (row >= nr + 1 || col >= nc + 1 || nd >= nd + 1)
		error("Matrix subscript out of bounds");
	else if (row <= 0 || col <= 0 || nd <= 0)
		error("Matrix subscript out of bounds");

	ind = (nd - 1)* (nr*nc) + (col - 1)*nc + (row - 1);
	return ets[ind];
}







class mat4d {

private:
	int nr;      // number of rows in the matrix 
	int nc;      // number of columns in the matrix 
	int ndim1;
	int ndim2;
	int ind;
	double* ets;   // entries of the matrix
public:
	mat4d(int, int, int, int); // constructor
	~mat4d() {
		delete[] ets;
	};

	double& operator() (int, int, int, int);

};


mat4d::mat4d(int row, int col, int nd1, int nd2) {
	nr = row;
	nc = col;
	ndim1 = nd1;
	ndim2 = nd2;
	ets = new double[nr* nc * ndim1*ndim2];
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

double& mat4d::operator() (int row, int col, int nd1, int nd2)
{
	if (row >= nr + 1 || col >= nc + 1 || nd1 >= nd1 + 1 || nd2 >= nd2 + 1)
		error("Matrix subscript out of bounds");
	else if (row <= 0 || col <= 0 || nd1 <= 0 || nd1 <= 0)
		error("Matrix subscript out of bounds");

	ind = (nd2 - 1)* (nr*nc*ndim1) + (nd1 - 1)* (nr*nc) + (col - 1)*nc + (row - 1);
	return ets[ind];
}
