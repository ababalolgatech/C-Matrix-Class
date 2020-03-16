
//#include"DataClass.h"
#include"GstatData.h"  // added in Geostat





#define	E  (2.71828182845904523536)
#define UNKNOWN  (-999999.0)
#define TINY (1e-20)

#define MAX(a,b) ((a>b) ? a : b)
#define MIN(a,b) ((a<b) ? a : b)
#define NINT(a)  ((a>0) ? (int)(a+0.5) :  (int)(a-0.5)) 

#define M 714025
#define IA 1366
#define IC 150889

// point is float & ipt is integers
class pt {
public:
	float x, y, z;
//	~pt() {};
};

class ipt {
public:
	int  x, y, z;

};

class variogram_model {
public:
	int   num_struct;
	float cmax, nugget;
	ivec     type;
	fvec    sill, a;
	fvec    ang1, ang2, ang3;
	
	fmat2d  dir_cosines;
	pt*		anis;

	variogram_model()
	{
		anis = new pt[10];
		dir_cosines.zeros(3,3);
		       type.zeros(10) ;
		       sill.zeros(10) ;
		          a.zeros(10) ;
				  ang1.zeros(10);
				  ang2.zeros(10);
				  ang3.zeros(10);
				
	};

	~variogram_model() {};

	};

class Search {
public:
	Search() {
		cosines.zeros(1,3,1,3);
	}
	fmat2d cosines;
	pt radius;
	int	rotation;
	float sang1, sang2, sang3;
};

class geotemplate {
public:
	ivec	irx, iry, irz, octant, index;
	fvec	distance;

	//~geotemplate() {};
};

class cond_point {
public :
	ivec		ix, iy, iz;
	fvec      	value;
	ivec        index;

//	~cond_point() {};
} ;


