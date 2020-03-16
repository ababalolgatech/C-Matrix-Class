#include"GstatData.h"
#include"GeoMath.h" // you need to add the header before you can inherit






#define	E  (2.71828182845904523536)
#define UNKNOWN  (-999999.0)
#define TINY (1e-20)

#define MAX(a,b) ((a>b) ? a : b)
#define MIN(a,b) ((a<b) ? a : b)
#define NINT(a)  ((a>0) ? (int)(a+0.5) :  (int)(a-0.5)) 

 #define M 714025   // defined in gmain
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
	float    ang1, ang2, ang3;

	fmat2d  dir_cosines;
	pt*		anis; // ranges in the x , y & z direction

	variogram_model()
	{
		anis = new pt[10];
		dir_cosines.zeros(3, 3);
		type.zeros(10);
		sill.zeros(10);
		a.zeros(10);
		//ang1.zeros(10);
		//ang2.zeros(10);
		//ang3.zeros(10);

	};

	~variogram_model() {};

};

class Search {
public:
	Search() {
		cosines.zeros(1, 3, 1, 3);
	}
	fmat2d cosines;
	pt radius; // ranges in the x , y & z direction
	int	rotation;
	float sanis1, sanis2;
	float sang1, sang2, sang3;
};

class geotemplate {
public:
	ivec	irx, iry, irz, octant, index;
	fvec	distance;

};

class cond_point {
public:
	ivec		ix, iy, iz;
	fvec      	value;
	ivec        index;

};

class Geostat :public MathsFun {
public:
	//---------------------------------------------------------------------
	void nscore_transform(fvec& value_org, float tmin, float tmax, int iwt, fvec& wt,fvec& value_trans, int ierror) {
		/*
	-----------------------------------------------------------------------

				  Transform Univariate Data to Normal Scores
				  ******************************************

	 This subroutibe takes "nd" data "value_org(i),i=1,...,nd" possibly weighted
	by "wt(i),i=,...,nd" and returns the normal scores transform N(0,1)
	 as "value_trans(i),i=1,...,nd".  The extra storage array "tmp" is required
	 so that the data can be returned in the same order (just in case
	 there are associated arrays like the coordinate location).



	 INPUT VARIABLES:

	   nd               Number of data (no missing values)
	   value_org(nd)        Data values to be transformed
	   tmin,tmax        data trimming limits
	   iwt              =0, equal weighted; =1, then apply weight
	   wt(nd)           Weight for each data (don't have to sum to 1.0)
	   tmp(nd)          Temporary storage space for sorting
	   lout             if > 0 then transformation table will be written



	 OUTPUT VARIABLES:

	   valueg(nd)          normal scores
	   ierror           error flag (0=error free,1=problem)



	 EXTERNAL REFERENCES:

	   gauinv           Calculates the inverse of a Gaussian cdf
	   sortem           sorts a number of arrays according to a key array



	-----------------------------------------------------------------------
		*/

		float EPSLON = 1.0e-20;
		int nd, ierr;
		fvec tmp;
		float cp, oldcp, twt;

		ierr = 0;

		//  Sort the data in ascending order and calculate total weight :
		nd = value_org.size();
		tmp.zeros(nd);  // allocate
		ierror = 0  ;
		   twt = 0.0;
		for (int i = 1; i < nd; i++){
			tmp(i) = std::abs(i);
			if (value_org(i) > tmin &&  value_org(i) < tmax)
			{
				if (iwt == 0)
				{
					twt = twt + 1;
				}
				else {
					twt = twt + wt(i);
				}
			}
		}
		if (nd < 1 || twt < EPSLON)
		{
			ierror = 1;
			return;
		}

		sortem(1, nd, value_org, 2, tmp,tmp);

		// Compute the cumulative probabilities :

		oldcp = 0.0;
		cp = 0.0;
		for (int i = 1; i < nd; i++)
		{
			   cp = cp + wt(i) / twt  ;
			wt(i) = (cp + oldcp) / 2.0;
			oldcp = cp;

		value_trans[i] = gauinv(double(wt[i]), ierr);

		}
		// Get the arrays back in original order :

		sortem(1, nd, tmp,3, value_org,value_trans);

	}
	//---------------------------------------------------------------------
	float gauinv(double p, int ierr) {
		/*
		----------------------------------------------------------------------
		  Computes the inverse of the standard normal cumulative distribution
		  function with a numerical approximation from : Statistical Computing,
		  by W.J.Kennedy, Jr. and James E.Gentle, 1980, p. 95.

		 INPUT / OUTPUT :

			  p = float precision cumulative probability value : dble(psingle)
			 xp = G^-1 (p)in single precision
		   ierr = 1 - then error situation(p out of range), 0 - OK
		 ---------------------------------------------------------------------
		*/

		//  Coefficients of approximation :  
		float lim, p0, p1, p2, p3, p4, q0, q1, q2, q3, q4, xp, pp, y;
		lim = 1.0E-10;
		p0 = -0.322232431088;
		p1 = -1.0;
		p2 = -0.342242088547;
		p3 = -0.0204231210245;
		p4 = -0.0000453642210148;
		q2 = 0.531103462366;
		q0 = 0.0993484626060;
		q1 = 0.588581570495;
		q3 = 0.103537752850;
		q4 = 0.0038560700634;




		// Check for an error situation:
		ierr = 1;
		if (p < lim)
		{
			xp = -1.0E10;
			return xp;
		}
		if (p > (1.0 - lim))
		{
			xp = 1.0E10;
			return xp;
		}

		ierr = 0;

		// Get k for an error situation:
		pp = p;
		if (p > 0.5) pp = 1 - pp;
		xp = (float)0.0;
		if (p == (double)0.5) return xp;

		// Approximate the function:
		y = sqrt(log(1.0 / (pp*pp)));
		xp = (float)(y + ((((y*p4 + p3)*y + p2)*y + p1)*y + p0) / ((((y*q4 + q3)*y + q2)*y + q1)*y + q0));
		if ((float)p == (float)pp) xp = -xp;

		// Return with G^-1(p):
		return xp;
	}
	//---------------------------------------------------------------------
	float gcum(float x) {
		/*
		Evaluate the standard normal cdf given a normal deviate x.gcum is
		the area under a unit normal curve to the left of x.The results are
		accurate only to about 5 decimal places.
		*/
		float z, t, gcum, e2;
		if (x < 0) { z = -x; }
		else { z = x; }  // float check

		t = 1 / (1 + 0.2316419 * z);
		gcum = t * (0.31938153 + t * (-0.356563782 + t * (1.781477937 +
			t * (-1.821255978 + t * 1.330274429))));


		if (z <= 6) { e2 = exp(-z * z / 2.)*0.3989422803; }
		else { e2 = 0; }

		gcum = 1.0 - e2 * gcum;
		if (x >= 0) { return gcum; }
		else { return 1.0 - gcum; }

	}
	//---------------------------------------------------------------------
	float powint(float xlow, float xhigh, float ylow, float yhigh, float xval, float potencia)
	{
		float EPSLON = 1.0E-20;

		float valor;

		if ((xhigh - xlow) < EPSLON)
			valor = (yhigh + ylow) / float(2.0);
		else
			valor = ylow + (yhigh - ylow)*(float)pow(((xval - xlow) / (xhigh - xlow)), potencia);

		return valor;
	}
	//---------------------------------------------------------------------
	float locate(fvec& xx, int n, int is, int ie, float x) {
		/*--------------------------------------------------------------------- -

			 Given an array "xx" of length "n", and given a value "x", this routine
			 returns a value "j" such that "x" is between xx(j) and xx(j + 1).
			 xx must be monotonic, either increasing or decreasing.
			 j = is - 1 or j = ie is returned to indicate that x is out of range.

			 Bisection Concept From "Numerical Recipes", Press et.al. 1986  pp 90.
			---------------------------------------------------------------------- -
			 Initialize lower and upper methods :
		*/


		int j, jl, ju, jm;

		if (is <= 0) {
			is = 0;
		};
		jl = is - 1;
		ju = ie;

		if (xx(n) <= x) {
			j = ie;
			return j;
		}


		// If we are not done then compute a midpoint :

	L: if (ju - jl > 0)
	{
		jm = int((ju + jl) / 2);
	};

	   //  Replace the lower or upper limit with the midpoint :

	   if ((xx[ie] > xx[is]) == (x > xx[jm]))
	   {
		   jl = jm;
	   }
	   else
	   {
		   ju = jm;
	   }

	   goto L;
	   j = jl; // Return with the array ivargex:
	   return j;
	}
	//---------------------------------------------------------------------
	float locate_c(fvec & xx, int  n, float x) {
		float ju, jm, jl, j;
		int ascnd;

		jl = 0;      // lower
		ju = n + 1;  // upper
		ascnd = std::int16_t(xx[n] >= xx[1]);
		while (ju - jl > 1) {
			jm = int((ju + jl) / 2);

			if (x >= xx[jm] == ascnd)
			{
				jl = jm;
			}
			else
			{
				ju = jm;
			}
		}
		if (x == xx[1]) {
			j = 1;
		}
		else if (x == xx[n])
		{
			j = n - 1;
		}
		else {
			j = jl;
		};


		return j;
	}
	//---------------------------------------------------------------------
	double acorni(fvec& ixv) {
		// I used the idum to know initialize ixv 

		double MAXOP1, MAXINT, KORDEI, temp, Acorni;
		KORDEI = 12;

		MAXOP1 = KORDEI + 1;
		MAXINT = 1073741824;


		for (int i = 0; i < KORDEI; i++) {
			ixv[i + 1] = ixv[i + 1] + ixv[i];
			if (ixv[i + 1] >= MAXINT) {
				ixv[i + 1] = ixv[i + 1] - MAXINT;
			}
		}

		Acorni = double(ixv(KORDEI)) / MAXINT;
		return Acorni;

	}
	//---------------------------------------------------------------------
	ivec genrandpath(int ns) {

		fvec sim;
		ivec order,tmp;

		tmp.zeros(ns);
		  sim.rand(ns);
		order.zeros(ns);
		for (int i = 1; i <= ns; i++)
		{
			order[i] = i;
		}


		//hpsort
		sortem(1, ns, sim,	2, order,tmp);
		return order;

	}

	void cosgs_error(char *err) {
		/* displays an error message and stops program */
		while (*err) printf("%c", *err++);
		std::cout << "" << std::endl;
		system("pause");
		exit(0);
		
	}
};
