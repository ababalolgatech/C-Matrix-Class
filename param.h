
#include<iostream>
#include<string>
#include"Geo11.h"

class param
{
public:
	 param();
	~param();

	std::string		fgeom, fvario, fdata, fout, fdebug;
	float			expected_value, value, ****cova_table;
	pt              delta, origin, coord;
	ipt             from, to, n_nodes, srch_dim;
	int             max_per_octant, dbg, ndata, nvar, nvar_sim, simulation_var,
		            ntemp, ncov, acov, lookup_table, index, nsim, i, max_data ;
	search			srch;
	variogram_model variogram;
	float           node;
	short           ***node_flag;
	long int        seed;
	geotemplate     tmp;


	void read_in_geometry(FILE *fp, FILE *fperr, long int *seed, int *nsim,
		int *max_data, pt *delta, pt *origin, ipt *n_nodes, ipt *from,
		ipt *to, search *srch, int *max_per_octant, int *dbg);
	



};

