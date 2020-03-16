#include"gmain.h"


class gsim :public gmain {

	/* 
	The following set of subroutines correspond to the simulation itself. All
    this subroutines share a number of variables which are explained below.
    This variables are capitalized to distinguish them from local variables 
	*/

   /*  Cond_Pts            Structure that stores the information corresponding
							   to the contioning points retained to krig a
							   particular node
	   Kriged_Mean         Result of OK at each simulated location
	   Kriged_Variance     Kriging variance at each simulated location
	   Kriged_Std          Square root of Kriged_Variance
	   Index_Number        Index number that varies randomly between 1 and
							  the total number of nodes
	   Krig_Indx[]         Intermediate vector used in the LU decomposition of
							   the kriging matrix
	   Krig_Lhs[,]         Matrix that stores the left hand side of the kriging
						   system
	   Krig_Nper           Intermediate value used in the LU decomposition of
							   the kriging matrix
	   Krig_Rhs[]          Vector that stores the right hand side of the kriging
						   system
	   Krig_Vect[]         Intermediate vector used in the LU decomposition of
							   the kriging matrix
	   Mod                 Smallest power of 2 that is larger than Total_Nodes.
							   (used to generate the random path)
	   Sim_Nodes           A 'ipt' that stores the size of the subarea being
							   simulated
	   Total_Nodes         Total number of nodes within simulation sub-area
   */


public:

	 ivec         Krig_Indx;
	 long int     Total_Nodes, Index_Number, Mod;
	 cond_point   Cond_Pts;
	 fvec         Krig_Rhs, Krig_Vect;
	 float        Krig_Nper;
	 ipt          Sim_Nodes;
	 fmat2d		  Krig_Lhs;
	 float        Kriged_Mean, Kriged_Variance, Kriged_Std;
	

	gsim() {};
	gsim(std::string fgeom, std::string fvario, std::string fdata
		, std::string fout, std::string fdebug) {};
	~gsim() {};

	

	void do_simulation(int isim);
	float obtain_SK_estimate(int npoints, int  ivar, fvec&expected_value);
	int find_conditioning_points(cond_point& Cond_Pts, ipt point);
	void build_kriging_system(cond_point& Cond_Pts, ipt point,int npoints, int ivar);
	void solve_kriging_system(int nequations);
	float obtain_SK_kriging_variance(cond_point& Cond_Pts, ipt point,int npoints, int ivar);
	void ludcmp(fmat2d& a, int n, ivec& indx, float d);
	void lubksb(fmat2d& a, int n, ivec& indx, fvec& b);
	float gasdev(long idum);
	long int initialize_random_path();
	void allocate_conditioning_data_matrices(int nvar, int max_data);
	ipt get_next_point(ipt from);
	ipt get_next_point(); // modified from GLSIB
	void simulate(std::string fgeom, std::string fvario, std::string fdata
		,std::string fout, std::string fdebug);
	void write_data(int nvar_sim);
	void display_report(int ireport, int iter, int TotalNodes);
};


void gsim::display_report(int ireport, int iter, int TotalNodes) {

	if (int(iter / ireport)*ireport == iter) 
	{
	std::cout << "Currently on node = " << iter << "  out of  = " << TotalNodes<<  std::endl;
	}
}

void gsim::write_data(int nvar_sim) {
	/* 
   This is an example of a subroutine to output the simulation results. In this
   case each of the output records contains as many columns as variables
   have been simulated  
   */

	std::cout << "WRITING DATA TO FILE = " << fout << std::endl;

	int iv;
	std::ofstream FL ("sgs.dat");

	/* output by layers, then by rows, then by columns */
    /* all variables in the same output record */

		
	for (int iiii = 1; iiii <= nvar_sim; iiii++)
	{
		iv = simulation_var[iiii];
		for (int iz = from.z; iz <= to.z; iz++)
		{
			for (int iy = from.y; iy <= to.y; iy++)
			{
				for (int ix = from.x; ix <= to.x; ix++)
				{
					FL << node(ix, iy, iz, iv) << std::endl;
				}
			}
				
		}
	}

	std::cout << "iv = " << iv << std::endl;

	//FL.close();
}

void gsim::simulate(std::string fgeom, std::string fvario, std::string fdata
	, std::string fout, std::string fdebug) {

	/*
	read_in_geometry(fgeom, seed, nsim, max_data, delta, origin, n_nodes,
	from, to, srch, max_per_octant, dbg);
	*/
	read_in_geometry(fgeom);

	 nvar = read_in_variograms(fvario,variogram);
  	 ncov = nvar * (nvar + 1) / 2;
	ndata = read_in_conditioning_data(fdata, coord);
	// allocating memory

	allocate_nodes(nvar, n_nodes); 
	// returns the allocated nodes

	allocate_covariance_tables(); 
	// returns the allocated cova_table

	ntemp = allocate_template();


	/* initialize variables and arrays */

	initialize_covariance_tables(variogram);
	acov = lookup_table(simulation_var[1],simulation_var[1]);

	initialize_template(acov);

	/* some consistency checks */

	a_few_checks(); 


	/* seed random number generator */
	if (seed > 0) 
	{
		seed = (-seed); 
	}

	rand0(seed);  // this should modify the seed


	FILE* fperr = NULL;
	for (int i = 1; i < nsim; i++)
	{
		initialize_nodes(coord,i);

		/* simulation starts */

		do_simulation(i);

		write_data(nvar_nsim); // note , I did not add headers or No of simulation



	}  

	exit(0); /* normal termination */
}

float gsim::gasdev(long idum) {
	static int iset = 0;
	static float gset;
	float fac, r, v1, v2;

	if (iset == 0)
	{
		do
		{
			v1 = 2.0*rand0(idum) - 1.0;
			v2 = 2.0*rand0(idum) - 1.0;
			r = v1 * v1 + v2 * v2;
		} while (r >= 1.0);
		 fac = sqrt(-2.0*log(r) / r);
		gset = v1 * fac;
		iset = 1;
		return v2 * fac;
	}
	else
	{
		iset = 0;
		return gset;
	}
}

void gsim::lubksb(fmat2d& a, int n, ivec& indx, fvec& b) {
	/* 
	this is the solution of the system LUx=b, in which LU is the LU
	decomposition of matrix 'a' 
	   */

	// don't change the body of this with some fancy curly braces !

	int i, ii = 0, ip, j;
	float sum;

	for (i = 1; i <= n; i++)
	{
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii)
			for (j = ii; j <= i - 1; j++) sum -= a(i,j) * b[j];
		else if (sum) ii = i;
		b[i] = sum;
	}
	for (i = n; i >= 1; i--)
	{
		sum = b[i];
		for (j = i + 1; j <= n; j++) sum -= a(i,j) * b[j];
		b[i] = sum / a(i,i);
	}
}

void gsim::ludcmp(fmat2d& a, int n, ivec& indx, float d) {

	/* this is the LU decomposition itself */

	int i, imax, j, k;
	float big, dum, sum, temp;

	/* 
	since this subroutine is called many times, the allocation of the vector
	'vv' is done outside. This vector was called Krig_Vect and is #define'd here
	as 'vv'.
	*/

#define vv Krig_Vect

	d = 1.0;
	for (i = 1; i <= n; i++)
	{
		big = 0.0;
		for (j = 1; j <= n; j++)
			if ((temp = fabs(a(i, j)) > big))
			{
				big = temp;
			};
		if (big == 0.0) perror("Singular matrix in routine LUDCMP");

		vv[i] = 1.0 / big;  // is it included in the if statement 
	}
	for (j = 1; j <= n; j++)
	{
		for (i = 1; i < j; i++)
		{
			sum = a(i, j);
			for (k = 1; k < i; k++) sum -= a(i, k) * a(k, j);
			a(i, j) = sum;
		}
		big = 0.0;
		for (i = j; i <= n; i++)
		{
			sum = a(i, j);
			for (k = 1; k < j; k++)
				sum -= a(i, k) * a(k, j);
			a(i, j) = sum;
			if ((dum = vv[i] * fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 1; k <= n; k++)
			{
				dum = a(imax, k);
				a(imax, k) = a(j, k);
				a(j, k) = dum;
			}
			d = -d;
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a(j, j) == 0.0) a(j, k) = TINY;
		if (j != n)
		{
			dum = 1.0 / (a(j, j));
			for (i = j + 1; i <= n; i++)
			{
				a(i, j) *= dum;
			}
		}
	}
#undef vv
}

float gsim::obtain_SK_kriging_variance(cond_point& Cond_Pts, ipt point, int npoints, int ivar) {

	float variance,temp;
	int dx, dy, dz, i, acov;

	acov = lookup_table(ivar, ivar);
	variance = cova_table(0, 0, 0,acov);
	for (i = 1; i <= npoints; i++)
	{
		dx = Cond_Pts.ix[i] - point.x;
		dy = Cond_Pts.iy[i] - point.y;
		dz = Cond_Pts.iz[i] - point.z;
		acov = lookup_table(ivar, Cond_Pts.index[i]);

		
		variance -= (dx < 0 ? cova_table(-dx,-dy,-dz,acov) : 
			cova_table(dx, dy, dz,acov)) * Krig_Rhs[i];
	};

	return variance;
}

void gsim::solve_kriging_system(int nequations) {
	/* solution of the kriging system by means of the LU decomposition */

	ludcmp(Krig_Lhs, nequations, Krig_Indx, Krig_Nper);
	lubksb(Krig_Lhs, nequations, Krig_Indx, Krig_Rhs);
}

void gsim::build_kriging_system(cond_point& Cond_Pts, ipt point, int npoints, int ivar) {
	/* 
	The left hand side matrix and the right hand side vector of the kriging
	system are built here using the covariance look up tables previously
	intialized 
	*/

	int i, j, dx, dy, dz, acov;

	for (i = 1; i <= npoints; i++)
	{
		for (j = 1; j <= i; j++)
		{
			dx = Cond_Pts.ix[i] - Cond_Pts.ix[j];
			dy = Cond_Pts.iy[i] - Cond_Pts.iy[j];
			dz = Cond_Pts.iz[i] - Cond_Pts.iz[j];
			acov = lookup_table(Cond_Pts.index[i], Cond_Pts.index[j]);
			Krig_Lhs(i, j) = dx < 0 ? cova_table(-dx, -dy, -dz,acov) :
				cova_table(dx, dy, dz,acov);
			Krig_Lhs(j, i) = Krig_Lhs(i, j);

	
		}
		dx = Cond_Pts.ix[i] - point.x;
		dy = Cond_Pts.iy[i] - point.y;
		dz = Cond_Pts.iz[i] - point.z;
		acov = lookup_table(ivar, Cond_Pts.index[i]);
		Krig_Rhs[i] = dx < 0 ? cova_table(-dx, -dy, -dz,acov) :
			cova_table(dx, dy, dz,acov);
	}
}

inline void gsim::do_simulation(int isim) {

	/* This function performs the actual simulation */

	/*
	This version of this subroutine only considers simple co-kriging
	to infer the parameters of the conditional probability distribution
	In the future, ordinary co-kriging will be implemented
	*/

	ipt         point;
	int         j, k, npoints, flag, report, iv, ivar;
	int         nxy, nxyz, nz;
	long int    i;
	ivec   index;

	
	

	/* 
	   point     stores the coordinates of the point being simulated
	   npoints   is the number of conditioning points used by kriging
	   report    is an integer used to feed back the user with periodic
				 information about the point that is being simulated
	*/
	printf("Starting simulation\n\n");


	Sim_Nodes.x = to.x - from.x + 1;
	Sim_Nodes.y = to.y - from.y + 1;
	Sim_Nodes.z = to.z - from.z + 1;
	Total_Nodes = (long int)Sim_Nodes.x*(long int)Sim_Nodes.y*(long int)Sim_Nodes.z;

	nxyz = Sim_Nodes.x*Sim_Nodes.y*Sim_Nodes.z;
	 nxy = Sim_Nodes.x*Sim_Nodes.y;
	  nz = Sim_Nodes.z;

	index = genrandpath(nxy);  // use GSLIB to generate Random Numbers
	

	gsim::allocate_conditioning_data_matrices(nvar, max_data);
	/*
	Index_Number = initialize_random_path();
	report = (log((double)Total_Nodes) - 2.5);
	if (report < 0)	{		report = 0;	}
	report = pow(E, (double)report);
	*/

	ireport = max(1, min(int(Total_Nodes / 10), 10000));
	
	printf("Simulation number %d\n\n", isim);
	
	for (i = 1; i <= nxy; i++)	{
		//std::cout << "i = " << i << " out of =  " << nxy << std::endl;
		display_report(ireport, i, Total_Nodes);
		for (int ii = 1; ii <= nz; ii++) {

		//	if ((i%report) == 1){printf("Now simulating point number %ld\r", i);}

			/* find the coordinates of the point to simulate */
			Index_Number = index[i] + (ii - 1)*nxy;
			
			point = get_next_point();

			/* do nothing if it is a conditioning data with values for all the variables */

			if (node_flag(point.x, point.y, point.z) == 1) continue;

			/* find conditioning points within search neighborhood */

			npoints = find_conditioning_points(Cond_Pts, point);

			for (iv = 1; iv <= nvar_nsim; iv++)
			{
				ivar = simulation_var[iv];
				// if previously simulated
				if (node(point.x, point.y, point.z, ivar) != UNKNOWN)
				{
					continue;
				}

				if (npoints > 0)
				{
					build_kriging_system(Cond_Pts, point, npoints, ivar);
					solve_kriging_system(npoints);
				} // check Line 426

				Kriged_Mean = obtain_SK_estimate(npoints, ivar, expected_value);
				Kriged_Variance = obtain_SK_kriging_variance(Cond_Pts, point, npoints, ivar);
				Kriged_Std = std::sqrt((double)Kriged_Variance);

				// The sqrt root is the reason why the results are different  iter = 19

				node(point.x, point.y, point.z, ivar) = gasdev(seed)*Kriged_Std + Kriged_Mean;

				npoints++;
				Cond_Pts.ix[npoints] = point.x;
				Cond_Pts.iy[npoints] = point.y;
				Cond_Pts.iz[npoints] = point.z;
				Cond_Pts.index[npoints] = ivar;
				Cond_Pts.value[npoints] = node(point.x, point.y, point.z, ivar);
			}
		}
		node_flag(point.x, point.y, point.z) = 1;
	}
	
}

float gsim::obtain_SK_estimate(int npoints, int  ivar, fvec&expected_value){
	/* 
	The linear simple kriging estimate of the kriging estimate of the
	conditional probability is computed here 
	*/

	int i;
	float estimate;

	estimate = expected_value[ivar];
	for (i = 1; i <= npoints; i++)
	{
		estimate += Krig_Rhs[i] * (Cond_Pts.value[i] -
			expected_value[Cond_Pts.index[i]]);	
	}
	return estimate;
}

void gsim::allocate_conditioning_data_matrices(int nvar, int max_data) {

	int maxpt = nvar * max_data + 1;
	Cond_Pts.ix.zeros(maxpt);
	Cond_Pts.iy.zeros(maxpt);
	Cond_Pts.iz.zeros(maxpt);
	Cond_Pts.value.zeros(maxpt);
	Cond_Pts.index.zeros(maxpt);

	Krig_Lhs.zeros(maxpt + 1, maxpt + 1);
	Krig_Rhs.zeros(maxpt + 1);
	Krig_Indx.zeros(maxpt + 1);
	Krig_Vect.zeros(maxpt + 1);

}

long int gsim::initialize_random_path() {
	/* 
	  This function initialize the value of Index_Number and it also finds the
	   smallest power of 2 that is larger than Total_Nodes. This value is used
	   in a congruential equation that will generate the random path
	*/

	Mod = 1;
	do
	{
		Mod = Mod*2;
	} 
	while (Mod <= Total_Nodes);
	Index_Number = Total_Nodes / 2;

// Bratley et al., 1983):
// indeXt =(5 x indext _! + l)mod 2n

	return Index_Number;

}

ipt gsim::get_next_point(ipt from){

	/*
	   This function returns the coordinates of the next point to be simulated.
	   Each node in the simulation sub-area can be uniquely determined by an index
	   number between 1 and Total_Nodes according to the following system:

			index           coordinates realtive to (from.x, from.y, from.z)
			-----           ------------------------------------------------
			  0                 (0,0,0)
			  1                 (0,0,1)
			  .                    .                    z increments fastest
			  .                    .
			Sim_Nodes.z-1       (0,0,Sim_Nodes.z)
			Sim_Nodes.z         (0,1,0)
			Sim_Nodes.z+1       (0,1,1)
			  .                    .                    y increments next
			  .                    .
	Sim_Nodes.y*Sim_Nodes.z-1   (0,Sim_Nodes.y,Sim_Nodes.z)
	Sim_Nodes.y*Sim_Nodes.z     (1,0,0)
			  .                    .                    x increments slowest
			  .                    .
		   Total_Nodes          (Sim_Nodes.x,Sim_Nodes.y,Sim_Nodes.z)
	*/


	ipt             point;
	long int        dummy, nyz;


	/* calculate coordinates relative to sub-area being simulated */

	    nyz = (long int)Sim_Nodes.y*(long int)Sim_Nodes.z;
	point.x = Index_Number / nyz;
	  dummy = Index_Number - point.x*nyz;
	point.y = dummy / Sim_Nodes.z       ;
	point.z = dummy - point.y*Sim_Nodes.z;

	/* calculate true grid coordinates */

	point.x += from.x;
	point.y += from.y;
	point.z += from.z;

	/* advance index to be used for next point */

	do
	{
		Index_Number = (Index_Number * 5 + 1) % Mod;
	} while (Index_Number >= Total_Nodes);

	return point;
}

ipt gsim::get_next_point() { // modified from GSLIB


	

	ipt             point;
	long int       nx, nxy;


	/* calculate coordinates relative to sub-area being simulated */
	
	nxy = Sim_Nodes.x*(long int)Sim_Nodes.y;
	nx = Sim_Nodes.x;
	// nyz = (long int)Sim_Nodes.y*(long int)Sim_Nodes.z;
	point.z = int((Index_Number - 1)/nxy ) + 1;
	point.y = int((Index_Number - (point.z-1)*nxy - 1) / nx) + 1;
	point.x = Index_Number - (point.z-1)*nxy - (point.y-1)*nx;
	
	point.x -= 1;
	point.y -= 1;
	point.z -= 1;
	/* calculate true grid coordinates */

	point.x += from.x;
	point.y += from.y;
	point.z += from.z;

	return point;
}

int gsim::find_conditioning_points(cond_point& Cond_Pts,ipt point) {

	/*
	For the node with coordinates stored in 'point', this function will find
	the nodes within the search neighborhood that have already being simulated
	or that are initial data. To do this the function makes use of the template
	previously calculated in which the relative coordinates of all possible
	points have being stored in order of distance starting with the closest
	point. There is a maximum number of points per octant that can be used,
	so the maximum number points to be used is 8*max_per_octant. The function
	returns the number of conditioning points found.
	*/

	int     i, j, jx, jy, jz, ioct, points_per_octant[8], point_counter, iv;

	/* initialize octant_counters and point counter */

	if (max_per_octant > 0)
		for (i = 0; i < 8; i++)
			points_per_octant[i] = 0;

	point_counter = 0;


	/* retrieve the points from the template in order of decreasing distance */

	for (i = 1; i <= ntemp && point_counter < max_data; i++)
	{
		j = tmp.index[i];
		/* obtain absolute coordinates */
		jx = tmp.irx[j] + point.x;
		jy = tmp.iry[j] + point.y;
		jz = tmp.irz[j] + point.z;

		/* if the point is outside the grid or has not been simulated yet discard */

		if (jx<1 || jy<1 || jz<1 ||
			jx>n_nodes.x || jy>n_nodes.y || jz>n_nodes.z ||
			node_flag(jx, jy, jz) == -1) continue;

		/* 
		if there are already more than max_per_octant points
		   in the octant to which the point belongs, discard
		   otherwise move the point to the Cond_Pts structure
		   and increase octant counter and point counter
		 */

		if (max_per_octant > 0)
		{
			ioct = tmp.octant[j];
			if (points_per_octant[ioct] >= max_per_octant)
			{
				continue;
			}
			else points_per_octant[ioct]++;
		}

		for (iv = 1; iv <= nvar; iv++)
		{

			if (node(jx, jy, jz,iv) == UNKNOWN) continue;
			point_counter++;

			Cond_Pts.ix[point_counter] = jx;
			Cond_Pts.iy[point_counter] = jy;
			Cond_Pts.iz[point_counter] = jz;
			Cond_Pts.value[point_counter] = node(jx, jy, jz,iv);
			Cond_Pts.index[point_counter] = iv;
		}


		/* check if there is any octant that is not complete */

		if (max_per_octant > 0)
		{
			for (ioct = 0; ioct < 8; ioct++)
			{
				if (points_per_octant[ioct] <= max_per_octant) break;
			}

			if (ioct < 8) continue;
			else break;
		}
	}

	return point_counter;
}




