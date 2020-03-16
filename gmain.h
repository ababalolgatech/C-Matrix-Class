
#include <Armadillo>
#include<iostream>
#include<fstream>
#include "math.h"
#include<string>
#include"Config.h"
#include"ConfigVario.h"
#include"GSLIB.h"
#include"SignalProcessing.h"



class gmain : public Geostat{
public:

	std::string		fgeom, fvario, fdata, fout, fdebug;
	fvec            expected_value, value, value_org; // these are vectors
	pt              delta, origin, *coord; // coord has to be cos we are definitign an array of coordinates
	ipt             from, to, n_nodes, srch_dim, dim;
	int             max_per_octant, dbg, ndata, nvar, nvar_nsim,
					ntemp, cov, nsim, i, max_data, ncov,acov;
	ivec			simulation_var, index;
	Search			srch;
	variogram_model* variogram;
	geotemplate     tmp;
	fmat4d          node, cova_table;  // note cova_table(x,y,z,ncov)
	smat3d          node_flag;
	imat2d	        lookup_table;
	long int        seed;
	float			ltail,utail,ltpar,utpar, zmin, zmax;
	float			tmin, tmax;
	int	            cosine_flag, nscore_trans, ireport;


	


	 gmain() {};
	~gmain() {};



	void read_in_geometry(std::string fp);
	int read_in_variograms(std::string fp,variogram_model*& variogram);
	void input_structure( std::string , variogram_model*& var_model, int i);
	int	findout_rotation(fmat2d& cosines);
	int read_in_conditioning_data(std::string fp, pt*& coord);
	void allocate_nodes(int nvar, ipt n_nodes);
	void allocate_covariance_tables();
	int allocate_template();
	void initialize_covariance_tables(variogram_model*& variogram);
	void a_few_checks();
	void initialize_nodes(pt*& coord, int isim);
	void write_out_results(FILE* fout, int isim);
	void initialize_template(int acov);
	float lookup_cova(int dx, int dy, int dz, int ind);
	void hpsort(int n, fvec&ra, ivec&index);
	pt rotate(pt vect, fmat2d& dir_cos);
	float exponential(float h, float sill, float a);
	float spherical(float h, float sill, float a);
	float gaussian(float h, float sill, float a);
	float linear(float h, float slope, float power);
	float modulus(pt vect, pt anis);
	float cova(pt point1, pt point2, variogram_model& var_model);
	float gamma0(pt point1, pt point2, variogram_model&  var_model);
	void setrot(float ang1, float ang2, float ang3,fmat2d& cosines);
	void backtr(float simval, int nt, fvec& value, fvec& value_trans, float& Backtr);
	
	//void  cosgs_error(char *err);
};




void gmain::backtr(float simval, int nt, fvec& value, fvec& value_trans, float& Backtr) {

	/*	               Back Transform Univariate Data from Normal Scores
					   *************************************************

			  This subroutine backtransforms a standard normal deviate from a
			  specified back transform table and option for the tails of the
			  distribution.Call once with "first" set to true then set to false
			 unless one of the options for the tail changes.

			  INPUT VARIABLES :

				simval             normal score value to be back transformed
				nt               number of values in the back transform table
				value(nt)        original data values that were transformed
				value_trans(nt)          the corresponding transformed values
				zmin, zmax        limits possibly used for linear or power model
				ltail            option to handle values less than value_trans(1) :
				ltpar            parameter required for option ltail
				utail            option to handle values greater than value_trans(nt) :
				utpar            parameter required for option utail

			 Parameters

			Value in the lower tail ? 1 = linear, 2 = power, (3 and 4 are invalid) :  */

	float  cdflo, cdfhi, cdfbt, cpow, lambda;
	int j;



	// values inthe lower tail
	if (simval <= value_trans[0])
	{
		Backtr = value[1];
		cdflo = gcum(value_trans[0]);
		cdfbt = gcum(simval);
		if (ltail == 1)
		{
			Backtr = powint(0.0, cdflo, zmin, value[1], cdfbt, 1.0);
		}
		else if (ltail == 2) {
			cpow = 1.0 / ltpar;
			Backtr = powint(0.0, cdflo, zmin, value[1], cdfbt, cpow);
		}
	}
	// values in the upper tail 
	else if (simval >= value_trans[nt])
	{
		Backtr = value(nt);
		cdfhi = gcum(value_trans[nt]);
		cdfbt = gcum(simval);
		if (ltail == 1)
		{
			Backtr = powint(cdfhi, 1.0, value[nt], zmax, cdfbt, 1.0);
		}
		else if (ltail == 2)
		{
			cpow = 1.0 / ltpar;
			Backtr = powint(cdfhi, 1.0, zmax, value[nt], cdfbt, cpow);
		}
		else if (utail == 4)
		{
			lambda = float(((pow(value[nt], utpar))*(1.0 - gcum(value_trans[nt]))));
			Backtr = float((pow((lambda /
				(1.0 - gcum(simval))), (1.0 / utpar))));
		}

	}
	else
	{
		//Value within the transformation table:
		j = locate_c(value_trans, nt, simval);
		j = max(min((nt - 1), j), 1);
		Backtr = powint(value_trans(j), value_trans(j + 1), value(j), value(j + 1), simval, 1.0);
	}
}
//-------------------------------------------------------------------
void gmain::setrot(float ang1, float ang2, float ang3,fmat2d& cosines) {
	/*
	Sets up an Anisotropic Rotation Matrix
	**************************************

	Sets up the matrix to transform cartesian coordinates to coordinates
	accounting for angles and anisotropy (see GSLIB manual for a detailed
	definition):


	INPUT PARAMETERS:

	ang1             Azimuth angle for principal direction
	ang2             Dip angle for principal direction
	ang3             Third rotation angle
	anis1            First anisotropy ratio
	anis2            Second anisotropy ratio
	MAXROT           maximum number of rotation matrices dimensioned
	rotmat           rotation matrices
	*/
	float DEG2RAD, alpha, beta, theta;
	float sina, sinb, cosb, cosa, sint, cost;
	float afac1, afac2, neq;
	
	/*
	Converts the input angles to three angles which make more
	mathematical sense :

	alpha   angle between the major axis of anisotropy and the
	E - W axis.Note : Counter clockwise is positive.
	beta    angle between major axis and the horizontal plane.
	(The dip of the ellipsoid measured positive down)
	theta   Angle of rotation of minor axis about the major axis
	of the ellipsoid.
	*/
	DEG2RAD = 3.141592654 / 180.0;
	
		if (ang1 >= 0.0 & ang1 < 270) 
		{
			alpha = (90 - ang1) * DEG2RAD;
		}
		else 
		{
			alpha = (450 - ang1) * DEG2RAD;
		}
		beta = -1.0 * ang2 * DEG2RAD;
		theta = ang3* DEG2RAD;

		//afac1 = 1.0 / float(max(anis1, TINY));
		//afac2 = 1.0 / float(max(anis2, TINY));

		// Get the required sines and cosines:
		sina = float(std::sin(alpha));
		sinb = float(std::sin(beta));
		sint = float(std::sin(theta));
		cosa = float(std::cos(alpha));
		cosb = float(std::cos(beta));
		cost = float(std::cos(theta));

		//  Construct the rotation matrix in the required memory:



		cosines(1, 1) = float(cosb * cosa);   // note its ijk 
		cosines(1, 2) = float(cosb * sina);
		cosines(1, 3) = float(-sinb);

		cosines(2, 1) = float(-cost * sina + sint * sinb*cosa);
		cosines(2, 2) = float(cost * cosa + sint * sinb * sina);
		cosines(2, 3) = float(sint * cosb);

		cosines(3, 1) = float(sint * sina + cost * sinb * cosa);
		cosines(3, 2) = float(-sint * cosa + cost * sinb * sina);
		cosines(3, 3) = float(cost * cosb);
	

}
//-------------------------------------------------------------------
float gmain::lookup_cova(int dx, int dy, int dz, int ind){
	/*
	This function returns the value stored in the covariance table 'cova_table'
	for a distance vector (measured in grid nodes) given by, dx, dy and dz.
	Since the covariance tables were computed for positive values of dx only,
	a radial symmetry is applied if dx<0 
	 */

	if (dx < 0) 
	{
		dx = -dx;
		dy = -dy;
		dz = -dz;
	}
	return cova_table(dx,dy,dz,ind);
}
//-------------------------------------------------------------------
void gmain::initialize_template(int acov) {

	/* 
   This subroutine finds the relative position of all the nodes within the
   Search nighborhood with respect to the central node. The distance to
   the center and the octant is also obtained. Finally the nodes in the
   template are sorted from closest to farthest to the center. The octants
   are numbered from 0 to 7 
   */

	int     ix, iy, iz, i, counter;
	pt      h;
	float   norm_h, oct_rot[3][3];

	/* 
	rotation to dealing the cartesian axes to the anisotropy
	axes to make the octant Search
	*/

	/*
	oct_rot[0][0]=oct_rot[0][1]=oct_rot[0][2]=oct_rot[1][2]=oct_rot[2][2]=.58;
	oct_rot[1][0]=oct_rot[2][1]= -0.79;
	oct_rot[1][1]=0.21; oct_rot[2][0]= -0.21;
	*/

	printf("Initializing template\n\n");

	counter = 0;
	for (ix = -dim.x; ix <= dim.x; ix++)
		for (iy = -dim.y; iy <= dim.y; iy++)
			for (iz = -dim.z; iz <= dim.z; iz++) 
			{
				h.x = ix * delta.x;
				h.y = iy * delta.y;
				h.z = iz * delta.z;
				h = rotate(h, srch.cosines);
				norm_h = h.x/srch.radius.x*h.x / srch.radius.x +
					     h.y/srch.radius.y*h.y / srch.radius.y +
					     h.z/srch.radius.z*h.z / srch.radius.z;
				if (norm_h < 1) 
				{
					/* h=rotate(h,oct_rot); */
					counter++;
					tmp.irx[counter] = ix;
					tmp.iry[counter] = iy;
					tmp.irz[counter] = iz;
					if (iz < 0)
					{
						tmp.octant[counter] = 4;
						if (ix <= 0 && iy > 0) tmp.octant[counter] = 1;
						if (ix > 0 && iy >= 0) tmp.octant[counter] = 2;
						if (ix < 0 && iy <= 0) tmp.octant[counter] = 3;
					}
					else
					{
						tmp.octant[counter] = 0;
						if (ix <= 0 && iy > 0) tmp.octant[counter] = 5;
						if (ix > 0 && iy >= 0) tmp.octant[counter] = 6;
						if (ix < 0 && iy <= 0) tmp.octant[counter] = 7;
					}
					/*            
					tmp.octant[counter]=
				     (h.x<0 ? 4 : 0) + (h.y<0 ? 2 : 0) + (h.z<0 ? 1 : 0); */
					/*                    tmp.distance[counter]=ix*ix+iy*iy+iz*iz; */
					/* now a structural distance is used instead of the euclidean distance used
					   before. */

					tmp.distance[counter] = -lookup_cova(ix, iy, iz,acov)
						+ 0.0001*norm_h;
				}
			}

	hpsort(counter, tmp.distance, tmp.index);

}

//-------------------------------------------------------------------
float rand0(long& idum) {

/*
   returns a uniform random deviate between 0.0 and 1.0.
   set idum to a negative value to initialize or reinitialize
   the sequence. This random number generator has period effectively
   infinite. Its principal limitation is that it returns one of
   only 714025 possible values, equally spaced as a comb in the
   interval [0,1). From page 212 of Numerical Recipes in C, Press
   et al., Cambridge University Press, 1988 
   */


	static long iy, ir[98];  /* the number 98, of elements in the shuffle
							   table is uninmportant */
	static int iff = 0;
	int j;
	//void error();

	if (idum < 0 || iff == 0)
	{
		iff = 1;
		if ((idum = (IC - idum) % M) < 0)
		{
			idum = -(idum);
		}
		for (j = 1; j <= 97; j++)/* initialize the shuffle table */
		{           
			 idum = (IA*idum + IC) % M;
			ir[j] = idum;
		}
		idum = (IA*idum + IC) % M;
		iy = idum;
	}               /* here is where it starts except on initialization */
	j = 1 + 97.0*iy / M;
	if (j > 97 || j < 1)
	{
		perror("RAND: This cannot happen");
	}
	   iy = ir[j];
	 idum = (IA*idum + IC) % M;
	ir[j] = idum;
	return (float) iy/M;
	// not sure why its different from : float(iy/M)
}
//-------------------------------------------------------------------
void gmain ::write_out_results(FILE* fout, int isim)

{

	/* This is an example of a subroutine to output the simulation results. In this
	   case each of the output records contains as many columns as variables
	   have been simulated   */

	int i, ix, iy, iz, iv;

	/* if this is the first simulation write header to output file */

	if (isim == 1) {
		fprintf(fout, "Output from GCOSIM3D, %d simulations of the subarea\n", nsim);
		fprintf(fout, "%d\n", nvar_nsim);
		for (i = 1; i <= nvar_nsim; i++)
			fprintf(fout, "Variable number %d\n", simulation_var[i]);
	}


	/* output by layers, then by rows, then by columns */
	/* all variables in the same output record */

	for (iz = from.z; iz <= to.z; iz++) 
	{
		for (iy = from.y; iy <= to.y; iy++) 
		{
			for (ix = from.x; ix <= to.x; ix++) 
			{
				for (i = 1; i <= nvar_nsim; i++)
				{
					iv = simulation_var[i];
//					fprintf(fout, "%7.3f ", node(iv,ix,iy,iz));
				}
//				fprintf(fout, "\n");
			}
		}
	}
}
//-------------------------------------------------------------------
void gmain::initialize_nodes(pt*& coord, int isim){
	/* 
	This subroutine assigns each conditioning data to the closest node in
	the simulation grid. It also print warning messages if two points are 
	assigned to the same node or if a point is outside of the simulation area
	*/
	int ix, iy, iz, i;   /* ix,iy,iz are the grid indices assigned to the datum */
	int iv;              /* the variable index */

	printf("Initializing nodes\n\n");

	for (ix = 1; ix <= n_nodes.x; ix++)
	{
		for (iy = 1; iy <= n_nodes.y; iy++) 
		{
			for (iz = 1; iz <= n_nodes.z; iz++) 
			{
				node_flag(ix,iy,iz) = -1;
				for (int iv = 1; iv <= nvar; iv++)
				{
					node(ix, iy, iz,iv) = UNKNOWN;
				}
			}
		}
	}


	for (i = 1; i <= ndata; i++)
	{
		ix = NINT((coord[i].x - origin.x) / delta.x) + 1;
		iy = NINT((coord[i].y - origin.y) / delta.y) + 1;
		iz = NINT((coord[i].z - origin.z) / delta.z) + 1;
		iv = index[i];

		if (   ix >= 1 && ix <= n_nodes.x
			&& iy >= 1 && iy <= n_nodes.y
			&& iz >= 1 && iz <= n_nodes.z)
		{    // point is within simulation area 
			if (node(ix, iy, iz,iv) == UNKNOWN)
			{
				node(ix, iy, iz, iv) = value[i];
			   node_flag(ix, iy, iz) = 0;
			}
		}
	}

}
//-------------------------------------------------------------------
inline void gmain::a_few_checks(){

	/*
	This function ends the initialization tasks with a few consistency checks
   in the input data. If any of this test is failed the program is aborted
   */

   int i;

   if (ndata < 0)
   {
	   cosgs_error("Number of data less than zero");
   }
   if (ntemp <= 0)
   {
	   cosgs_error("Search neighborhood too small, no data points in template");
   }
   if (nvar <= 0)
   {
	   cosgs_error("Number of variables less than or equal to zero");
   }
   if (nvar_nsim > nvar)
   {
	   cosgs_error("Number of simulated variables larger than number of variables\n");
   }
   if (n_nodes.x <= 0 || n_nodes.y <= 0 || n_nodes.z <= 0)
   {
	   cosgs_error("Error in number of nodes, no nodes in one direction");
   }
   if (from.x > to.x || from.y > to.y || from.z > to.z)
   {
	   cosgs_error("Incosistency between \'from\' and \'to\'");
   }
   if (from.x <= 0 || from.x > n_nodes.x || from.y<0 || from.y>n_nodes.y ||
	   from.z <= 0 || from.z > n_nodes.z)
   {
	   std::cout << "Note : origin corresponds to node(1,1,1)" << std::endl;
	   cosgs_error("\'from\' out of limits");
   }
   if (to.x <= 0 || to.x > n_nodes.x || to.y <= 0 || to.y > n_nodes.y || to.z <= 0 ||
	   to.z > n_nodes.z)
   {
	   std::cout << "Note : to corresponds to node(n_nodes.x,n_nodes.y,n_nodes.z)" << std::endl;
	   cosgs_error("\'to\' out of limits");
   }
   if (delta.x == 0 || delta.y == 0 || delta.z == 0)
   {
	   cosgs_error("Size of grid block equal to zero");
   }
}
//-------------------------------------------------------------------
void gmain::initialize_covariance_tables(variogram_model*& variogram){
	
	/* 
   This subroutine initializes the covariance tables. There are (nvar*(nvar+1)/2)
   covariance tables. The size of the smallest parallelepiped that
   includes the search neighborhood is passed in the 'ipt' srch_dim. The
   covariances are computed only for positive values of dx, since an
   ergodic covariance is assumed 
   */

	pt  point1, point2;
	int ix, iy, iz, i, j, counter = 0;

	printf("Initializing covariance tables. Table:   \n");

	point1.x = 0;
	point1.y = 0;
	point1.z = 0;

	/* 
	   covariance values are computed for vectors from grid node (0,0,0)
	   to grid node (ix,iy,iz) as (ix,iy,iz) sweeps the search nieghborhood 
	*/

	for (i = 1; i <= ncov; i++) 
	{
		printf("-------->   %d .. \n", i);
		for (ix = 0; ix <= 2 * srch_dim.x; ix++)		
			for (iy = -2 * srch_dim.y; iy <= 2 * srch_dim.y; iy++)	
				for (iz = -2*srch_dim.z; iz <= 2 * srch_dim.z; iz++) 
				{
					point2.x = ix * delta.x;
					point2.y = iy * delta.y;
					point2.z = iz * delta.z;
 	  cova_table(ix,iy,iz,i) = cova(point1, point2, variogram[i]);
				}
			
		
	}



	for (i = 1; i <= ncov; i++) 
	{
		for (j = 1; j <= i; j++)
		{
			counter++;
			lookup_table(i,j) = lookup_table(j,i) = counter;
		}
	}
};

//-------------------------------------------------------------------
int gmain::allocate_template() {
	
	/*
	This subroutine allocates the members of the 'tmp' structure. The 'tmp'
   structure contains all the information associated with the search template,
   that is, relative coordinates of all the points within a search
   neighborhood and distance to the point around which the search is taking
   place, as well as the octant to which each point belongs. There is also
   an 'index' member to retrieve the template points in order, from the
   closest to the farthest 
   */


	int         ix, iy, iz, i, ntemp;
	pt          h;
	float       norm_h;

	printf("Allocating template\n\n");

	ntemp = 0;
	for (ix = -dim.x; ix <= dim.x; ix++) {
		for (iy = -dim.y; iy <= dim.y; iy++) {
			for (iz = -dim.z; iz <= dim.z; iz++)
			{
				h.x = ix * delta.x;
				h.y = iy * delta.y;
				h.z = iz * delta.z;
				  h = rotate(h, srch.cosines);
				h.x /= srch.radius.x;
				h.y /= srch.radius.y;
				h.z /= srch.radius.z;
				norm_h = h.x*h.x + h.y*h.y + h.z*h.z;
				ntemp += (norm_h < 1) ? 1 : 0;
			}
		}
	}

	tmp.irx.zeros(1,ntemp);
	tmp.iry.zeros(1,ntemp);
	tmp.irz.zeros(1,ntemp);
	tmp.octant.zeros(1,ntemp);
	tmp.index.zeros(1,ntemp);
	tmp.distance.zeros(1,ntemp);

	//fprintf(fperr, "\n\nNumber of points in template -------> %d\n\n", ntemp);
	return ntemp;
}
//-------------------------------------------------------------------
void gmain::allocate_covariance_tables(){

	/* 
   This subroutine returns the pointer to a tensor of dimension 4
   that will contain the covariance tables. The dimensions of the covariance
   tables is determined by the smallest parallelepiped that can contain
   the search neighborhood 'srch'. This dimensions are computed according
   to the rotation flag (refer to the 'isim3d.h' for description of the
   structure 'search'. An ergodic covariance is assumed, therfore
   cova(dx,dy,dz)= cova(-dx,-dy,-dz) so that it needs to be computed only
   for positive values of dx 
   */

   /* optimize the dimensions of the covariance tables */

	delta.x = fabs((double)delta.x);
	delta.y = fabs((double)delta.y);
	delta.z = fabs((double)delta.z);


	switch (srch.rotation)
	{
	case 0: dim.x = NINT(srch.radius.x / delta.x);             /* no rotation */
		dim.y = NINT(srch.radius.y / delta.y);
		dim.z = NINT(srch.radius.z / delta.z);
		break;
	case 1: dim.x = NINT(srch.radius.x / delta.x);       /* rotation around x */
		dim.y = NINT(MAX(srch.radius.y, srch.radius.z) / delta.y);
		dim.z = NINT(MAX(srch.radius.y, srch.radius.z) / delta.z);
		break;
	case 2: dim.x = NINT(MAX(srch.radius.x, srch.radius.z) / delta.x);
		dim.y = NINT(srch.radius.y / delta.y);       /* rotation around y */
		dim.z = NINT(MAX(srch.radius.x, srch.radius.z) / delta.z);
		break;
	case 3: dim.x = NINT(MAX(srch.radius.x, srch.radius.y) / delta.x);
		dim.y = NINT(MAX(srch.radius.x, srch.radius.y) / delta.y);
		dim.z = NINT(srch.radius.z / delta.z);       /* rotation around z */
		break;
	case 4: dim.x = NINT(MAX(srch.radius.x, MAX(srch.radius.y, srch.radius.z))
		/ delta.x);
		dim.y = NINT(MAX(srch.radius.x, MAX(srch.radius.y, srch.radius.z))
			/ delta.y);
		dim.z = NINT(MAX(srch.radius.x, MAX(srch.radius.y, srch.radius.z))
			/ delta.z);
		break;
	default: perror("rotation flag out of limits");

	}
	
	srch_dim.x = dim.x;
	srch_dim.y = dim.y;
	srch_dim.z = dim.z;
	// Create copy constructor later


	cova_table.zeros(0, 2*dim.x, -2*dim.y, 2*dim.y, -2*dim.z, 2*dim.z,1, nvar);
	/*
	 note the diference variable is placed in the fourth dimension for convience
	 tensor(0,2*dim->x,-2*dim->y,2*dim->y,-2*dim->z,2*dim->z);
	 from index -2*dim to 2dim
	*/

		lookup_table.zeros(ncov, ncov);

};
//-------------------------------------------------------------------
 void gmain::allocate_nodes(int nvar, ipt n_nodes) 
 { 

	     node.zeros(n_nodes.x, n_nodes.y, n_nodes.z,nvar);
		 // note the difference.. variables are placed in the fourth dimension
	node_flag.zeros(n_nodes.x, n_nodes.y, n_nodes.z);
}
//-------------------------------------------------------------------
int gmain::read_in_conditioning_data(std::string fp, pt*& coorde) {

	arma::mat temp;  // to read in the data 
	int   i, ndata, iwt;
	fvec weight;


	printf("Reading hard conditinioning data\n\n");
	temp.load(fp); // loading the ascii file with armadillo
	ndata = temp.n_rows;
	value_org.zeros(ndata);
	value.zeros(ndata);
	index.zeros(ndata);
	coord = new pt[ndata];
	coord--;

	for (int i = 1; i <= ndata; i++)
	{
		coord[i].x = temp(i-1, 0);
		coord[i].y = temp(i-1, 1);
		coord[i].z = temp(i-1, 2);
   value_org[i] = temp(i-1, 4);
		 index[i] = temp(i-1, 3);
	}

	if (nscore_trans = 1) 
	{
		value = value_org;
	}
	else
	{
		weight.ones(ndata);
		iwt = 0;
		int ierror;
		nscore_transform(value_org, tmin, tmax, iwt,weight, value,ierror);
	}

	return ndata;
};
//-------------------------------------------------------------------
void gmain::input_structure(std::string fp, variogram_model*& var_model, int ind) {

	Config_vario cfg_vec(fp);
	ConfigFile cfg(fp);
	// auto nugget_tmp = cfg_vec.return_value("nugget");

	var_model[ind].nugget = cfg.getValueOfKey<float>("nugget");
	var_model[ind].cmax = cfg.getValueOfKey<float>("cmax");
	var_model[ind].num_struct = cfg.getValueOfKey<float>("num_struct");

	auto type_tmp = cfg_vec.return_value("type");
	auto sill_tmp = cfg_vec.return_value("sill");
	auto anisx_tmp = cfg_vec.return_value("anis_x");
	auto anisy_tmp = cfg_vec.return_value("anis_y");
	auto anisz_tmp = cfg_vec.return_value("anis_z");
	auto c1 = cfg_vec.return_value("c1");
	auto c2 = cfg_vec.return_value("c2");
	auto c3 = cfg_vec.return_value("c3");


	for (int i = 0; i < var_model[ind].num_struct; i++)
	{
		var_model[ind].type[i+1] = int(type_tmp[i]);
		var_model[ind].sill[i+1] = sill_tmp[i];
		var_model[ind].anis[i+1].x = anisx_tmp[i];
		var_model[ind].anis[i+1].y = anisy_tmp[i];
		var_model[ind].anis[i+1].z = anisz_tmp[i];
	}



	for (i = 1; i <= var_model[ind].num_struct; i++)
	{
		     var_model[ind].a[i] = var_model[ind].anis[i].x;
		var_model[ind].anis[i].x = var_model[ind].a[i] / var_model[ind].anis[i].x;
		var_model[ind].anis[i].y = var_model[ind].a[i] / var_model[ind].anis[i].y;
		var_model[ind].anis[i].z = var_model[ind].a[i] / var_model[ind].anis[i].z;
	}


	if (cosine_flag == 0)
	{
		var_model[ind].dir_cosines(1, 1) = c1[0];
		var_model[ind].dir_cosines(1, 2) = c1[1];
		var_model[ind].dir_cosines(1, 3) = c1[2];

		var_model[ind].dir_cosines(2, 1) = c2[0];
		var_model[ind].dir_cosines(2, 2) = c2[1];
		var_model[ind].dir_cosines(2, 3) = c2[2];

		var_model[ind].dir_cosines(3, 1) = c3[0];
		var_model[ind].dir_cosines(3, 2) = c3[1];
		var_model[ind].dir_cosines(3, 3) = c3[2];
	}
	else
	{
    setrot(var_model[ind].ang1, var_model[ind].ang2, var_model[ind].ang3, var_model[ind].dir_cosines);
	}

	//delete cfg;
}
//-------------------------------------------------------------
int gmain::read_in_variograms(std::string fp, variogram_model*& variogram) {
	int i, nvar;

	printf("Reading Variogram\n\n");

	ConfigFile cfg(fp);
	Config_vario cfg_vec(fp);
	
	nvar_nsim = cfg.getValueOfKey<int>("nvar_nsim");
	    nvar = cfg.getValueOfKey<int>("nvar");
	
		auto tmp_expected_val = cfg_vec.return_value("expected_value");
		auto tmp_simulation_var = cfg_vec.return_value("simulation_var");

		variogram = new variogram_model[(nvar*(nvar + 1) / 2)];
		variogram--;  // first element becomes 1

		expected_value.zeros(1,nvar);
		simulation_var.zeros(1,nvar_nsim);  


		// Note : the cfg_vec output a vector that starts at 0 while my in-built 

			for (int i = 0; i < nvar_nsim; i++)
			{
				expected_value[i + 1] = tmp_expected_val[i];
				simulation_var[i + 1] = tmp_simulation_var[i];
			}


	for (i = 1; i <= (nvar*(nvar + 1)) / 2; i++)
	{
		
		input_structure(fp, variogram,i);

	}

	//std::cout << "lots of bugs to fix" << std::endl;
	// https://youtu.be/AHxCQKQLpqY

	return nvar;
}
//-------------------------------------------------------------------
int gmain::findout_rotation(fmat2d& cosines) {
	int rotation;

	/* no rotation --- identity matrix */


	if (
		fabs((double)cosines(1, 1) - 1.0) < TINY &&
		fabs((double)cosines(2, 2) - 1.0) < TINY &&
		fabs((double)cosines(3, 3) - 1.0) < TINY &&
		fabs((double)cosines(1, 2)) < TINY &&
		fabs((double)cosines(1, 3)) < TINY &&
		fabs((double)cosines(2, 1)) < TINY &&
		fabs((double)cosines(2, 3)) < TINY &&
		fabs((double)cosines(3, 1)) < TINY &&
		fabs((double)cosines(3, 2)) < TINY)
	{
		rotation = 0;
	}
	/* rotation about the x-axis */

	else if (fabs((double)cosines(1, 1) - 1.0) < TINY &&
		fabs((double)cosines(1, 2)) < TINY &&
		fabs((double)cosines(1, 3)) < TINY &&
		fabs((double)cosines(2, 1)) < TINY &&
		fabs((double)cosines(3, 1)) < TINY)
	{
		rotation = 1;
	}

	/* rotation about the y-axis */

	else if (fabs((double)cosines(2, 2) - 1.0) < TINY &&
		fabs((double)cosines(2, 1)) < TINY &&
		fabs((double)cosines(2, 3)) < TINY &&
		fabs((double)cosines(1, 2)) < TINY &&
		fabs((double)cosines(3, 2)) < TINY)
	{
		rotation = 2;
	}

	/* rotation about the z-axis */

	else if (fabs((double)cosines(3, 3) - 1.0) < TINY &&
		fabs((double)cosines(3, 1)) < TINY &&
		fabs((double)cosines(3, 2)) < TINY &&
		fabs((double)cosines(1, 3)) < TINY &&
		fabs((double)cosines(2, 3)) < TINY)
	{
		rotation = 3;
	}

	else rotation = 4;

	return (rotation);
}
//-------------------------------------------------------------------
void gmain::read_in_geometry(std::string fp) {

	printf("Reading Geometry\n\n");

	ConfigFile cfg(fp);
	seed = cfg.getValueOfKey<float>("seed");
	nsim = cfg.getValueOfKey<float>("nsim");
	// Grid definition
	delta.x = cfg.getValueOfKey<float>("xsiz");
	delta.y = cfg.getValueOfKey<float>("ysiz");
	delta.z = cfg.getValueOfKey<float>("zsiz");
	origin.x = cfg.getValueOfKey<float>("xmn");
	origin.y = cfg.getValueOfKey<float>("ymn");
	origin.z = cfg.getValueOfKey<float>("zmn");
	n_nodes.x = cfg.getValueOfKey<float>("nx");
	n_nodes.y = cfg.getValueOfKey<float>("ny");
	n_nodes.z = cfg.getValueOfKey<float>("nz");
	from.x = cfg.getValueOfKey<float>("from_x");
	from.y = cfg.getValueOfKey<float>("from_y");
	from.z = cfg.getValueOfKey<float>("from_z");
	to.x = cfg.getValueOfKey<float>("to_x");
	to.y = cfg.getValueOfKey<float>("to_y");
	to.z = cfg.getValueOfKey<float>("to_z");


	// additional parameters
	cosine_flag = cfg.getValueOfKey<int>("cosine_flag");
	nscore_trans = cfg.getValueOfKey<int>("nscore_trans");

	if (nscore_trans == 1)
	{
		zmin = cfg.getValueOfKey<float>("zmin");
		zmax = cfg.getValueOfKey<float>("zmax");
		ltail = cfg.getValueOfKey<float>("ltail");
		ltpar = cfg.getValueOfKey<float>("ltpar");
		utail = cfg.getValueOfKey<float>("utail");
		utpar = cfg.getValueOfKey<float>("utpar");
		tmin = cfg.getValueOfKey<float>("tmin");
		tmax = cfg.getValueOfKey<float>("tmax");		
	}


	// search radii
	if (cosine_flag == 0) 
	{

		srch.radius.x = cfg.getValueOfKey<float>("radius_x");
		srch.radius.y = cfg.getValueOfKey<float>("radius_y");
		srch.radius.z = cfg.getValueOfKey<float>("radius_z");


		srch.cosines(1, 1) = cfg.getValueOfKey<float>("C11");
		srch.cosines(1, 2) = cfg.getValueOfKey<float>("C12");
		srch.cosines(1, 3) = cfg.getValueOfKey<float>("C13");

		srch.cosines(2, 1) = cfg.getValueOfKey<float>("C21");
		srch.cosines(2, 2) = cfg.getValueOfKey<float>("C22");
		srch.cosines(2, 3) = cfg.getValueOfKey<float>("C23");

		srch.cosines(3, 1) = cfg.getValueOfKey<float>("C31");
		srch.cosines(3, 2) = cfg.getValueOfKey<float>("C32");
		srch.cosines(3, 3) = cfg.getValueOfKey<float>("C33");
	}
	else 
	{
		setrot(srch.sang1, srch.sang2, srch.sang3, srch.cosines);
	}

	 srch.rotation = findout_rotation(srch.cosines); // debug
	max_per_octant = cfg.getValueOfKey<float>("max_per_octant");
	      max_data = cfg.getValueOfKey<float>("max_data");
	           dbg = cfg.getValueOfKey<float>("dbg");

	// work on writing to debug files later



}
//-------------------------------------------------------------------
float gmain::cova(pt point1, pt point2, variogram_model& var_model) {
	/* 
	This function returns the covariance between points 'point1' and 'point2'
	according to the variogram model var_model. The covariance is obtained from
	the variogram value is substracted 
	*/

	
	float cov = var_model.nugget;
	
	
	for (i = 1; i <= var_model.num_struct; i++) {
		switch (var_model.type[i]) 
		{
		case 1:
		case 2:
		case 3:
			cov += var_model.sill[i];
			break;
		case 4: cov += var_model.cmax;
		}
	}
	
	cov = cov - gamma0(point1, point2, var_model);
	return cov;
}
//-------------------------------------------------------------------
float gmain::gamma0(pt point1, pt point2, variogram_model& var_model) {
	/*
	   This function returns the variogram value between points 'point1' & 'point2'
	   according to the variogram model var_model. The four, most common 
	   variogram models are considered.
	*/

	float  value, h;
	int    i;
	pt     delta;

	/*	obtain separation vector between point1 and point2 and store it in point1	*/

	delta.x = point1.x - point2.x;
	delta.y = point1.y - point2.y;
	delta.z = point1.z - point2.z;

	delta = rotate(delta, var_model.dir_cosines);

	h = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
	if (h <= 1e-10)
	{
		return 0.;
	}

	value = var_model.nugget;                           // add nugget effect 




	for (i = 1; i <= var_model.num_struct; i++) // loop through all structures 
	{
		// first obtain the moudulus of the rotated vector after applying
		//  the anisotropy corrections 

		h = modulus(delta, var_model.anis[i]);

		if (h <= 1.e-20F) continue; // for very small distances the variogram is 0 


		switch (var_model.type[i])
		{
		case 1: value += spherical(h, var_model.sill[i], var_model.a[i]);
			break;

		case 2: value += exponential(h, var_model.sill[i], var_model.a[i]);
			break;

		case 3: value += gaussian(h, var_model.sill[i], var_model.a[i]);
			break;

		case 4: value += linear(h, var_model.sill[i], var_model.a[i]);
			break;

		default: perror("Undefined variogram type");
		}
	}



	return value;

	// first obtain the moudulus of the rotated vector after applying
	 //  the anisotropy corrections 


}
//-------------------------------------------------------------------
float gmain::modulus(pt vect, pt anis) {
	/* apply anisotropy ratios and compute the modulus of a vector */

	float h;

	h = vect.x*vect.x*anis.x*anis.x +
		vect.y*vect.y*anis.y*anis.y +
		vect.z*vect.z*anis.z*anis.z;
	h = sqrt((double)h);

	return h;
}
//-------------------------------------------------------------------
float gmain::linear(float h, float slope, float power) {
	/* linear isotropic variogram */

	//double pow;

	return slope * std::pow((double)h, (double)power);
}
//-------------------------------------------------------------------
float gmain::gaussian(float h, float sill, float a) {
	/* gaussian isotropic variogram

		double exp();  */

	h = h / a;
	return (sill*(1.0 - exp((double)-h * h)));
}
//-------------------------------------------------------------------
float gmain::spherical(float h, float sill, float a) {
	/* spherical isotropic variogram */

	h = h / a;
	return (h < 1.0 ? sill * h*(1.5 - 0.5*h*h) : sill);
}
//-------------------------------------------------------------------
float gmain::exponential(float h, float sill, float a) {
	/* exponential isotropic variogram

		double exp();  */

	h = h / a;
	return (sill*(1.0 - cos((double)h*3.1415927)));
}
//-------------------------------------------------------------------
pt gmain::rotate(pt vect, fmat2d& dir_cos) {
	/* 
	this function returns the vector 'vect' after rotation 
	according to the direction cosines in the 3 by 3 matrix dir_cos. 
	The rotated vector  returns in the same vector vect 
	*/

	pt  rotated;

	rotated.x = vect.x*dir_cos(1,1) +
				vect.y*dir_cos(1,2) +
				vect.z*dir_cos(1,3);

	rotated.y = vect.x*dir_cos(2,1) +
				vect.y*dir_cos(2,2) +
				vect.z*dir_cos(2,3);

	rotated.z = vect.x*dir_cos(3,1) +
				vect.y*dir_cos(3,2) +
				vect.z*dir_cos(3,3);

	return rotated;
}
//---------------------------------------------------------------------
void gmain::hpsort(int n, fvec&ra, ivec&indx){
	/* Replaces qcksort, that has failed in some very large arrays */

	int i, ir, j, l, indxt;
	float rra;
	int x, y;

	for (j = 1; j <= n; j++) 
	{
		indx[j] = j; 
	}

	if (n < 2) return;

	l = (n >> 1) + 1;
	ir = n;
	for (;;) {
		if (l > 1) 
		{
			y = --l;
			rra = ra[y];
			indxt = indx[y];
		}
		else 
		{
			rra = ra[ir];
			indxt = indx[ir];
			ra[ir] = ra[1];
			indx[ir] = indx[1];
			if (--ir == 1) 
			{
				ra[1] = rra;
			  indx[1] = indxt; break;
			}
		}
		i = l; j = l + l;
		while (j <= ir) 
		{
			if (j < ir && ra[j] < ra[j + 1]) j++;
			if (rra < ra[j]) 
			{
				ra[i] = ra[j];
			  indx[i] = indx[j];
				i = j;
				j <<= 1;
			}
			else j = ir + 1;
		}
		ra[i] = rra;
	  indx[i] = indxt;
	}
}