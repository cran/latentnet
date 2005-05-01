int Z_up(int *heads, int *tails, int n_edges, double **Z, double zdelta, 
	 double prior_mu, double prior_sd, int g, int k,
	 double **Znew, double *beta, int p, int dir,
	 double *llk, double ***X);
double loglike_y(int *heads, int *tails, int n_edges, int g, double **Z, 
		 int k, double *beta, int p, int dir, double ***X);
int procr_transform(double **Z, double **Zo, int g, int k, double **pZ);
int beta_up(int *heads, int *tails, int n_edges, double **Z, int g, int k, 
	    double prior_mu, double prior_sd, 
	    double *llk, double *beta, int p, int dir,
	    double ***X, double bdelta);
void neg_dist_Z(double **Z, int g, int dim, double **dZ);
