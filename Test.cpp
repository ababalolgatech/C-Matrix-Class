
#include"gsim.h"



using namespace std;

int main(){



	std::string fgeom, fvario, fdata, fout, fdebug;
	fgeom  = "geom.par"; 
	fvario = "gvar2.par";
	fdata  = "gcosim3d.dat";
	fout   = "sgs.dat";
	fdebug = "sgs_debug.dat";
	FILE* Debug;


	gsim sgs;
	sgs.simulate(fgeom, fvario, fdata, fout, fdebug);
	
	
	system("pause");
		return 0;
}



/*
1. How do de-allocate memory of the vectors & matrices created
2. Test all the functions 
*/

/*
what is the best way to design a class that has many class members
*/

/*

https://www.geeksforgeeks.org/abstraction-in-c/
https://www.geeksforgeeks.org/encapsulation-in-c/

You only need to pass in variables into member functions when are
using abstractions i.e declaring member variables as private
and modifying them with member functions.
*/

/*
vectors and matrices  - pass in the variables
member-defined data like pt & variogram  *&
*/

/*
To do 
Use enums for the different variogram model
*/