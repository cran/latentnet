int Z_up2(int *heads, int *tails, int n_edges, double **Z, double zdelta, 
	  double prior_mu, double prior_sd, int g, int k,
	  double **Znew, double *beta, int p, int dir, double ***X, 
	  double *mu,double *Sigma,int *Ki, int ng, double *llk,
	  int chisqprop, double thetaprop, double **Z_mle,
	  double **A, double **tZ, double **tZo, double ** Ahalf, 
	  double **AhalfInv, double **tptrans, double **eAvectors, 
	  double **eADvalues, double **teAvectors, double *avZ, 
	  double *avZo, double *eAvalues, double **kk_helper, 
	  double **gg_helper, double **gk_helper, 
	  double **kg_helper, double **kk2_helper);

double loglike2_y(int *heads, int *tails, int n_edges, int g, double **Z, 
		  int k, double *beta, int p, int dir, double ***X, double *mu,
		  double *Sigma, int *Ki, int ng);

int beta_up2(int *heads, int *tails, int n_edges, double **Z, int g, int k, 
	     double prior_mu, double prior_sd, 
	     double *llk, double *beta, int p, int dir, double ***X, 
	     double bdelta, double *mu,double *Sigma,int *Ki, int ng);

void neg_dist_Z2(double **Z, int g, int dim, double **dZ);

int procr_transform2(double **Z, double **Zo, int g, int k, double **pZ,
		     double **A, double **tZ, double **tZo, double ** Ahalf, 
		     double **AhalfInv, double **tptrans, double **eAvectors, 
		     double **eADvalues, double **teAvectors, double *avZ, 
		     double *avZo, double *eAvalues, double **kk_helper, 
		     double **gg_helper, double **gk_helper, 
		     double **kg_helper, double **kk2_helper);
