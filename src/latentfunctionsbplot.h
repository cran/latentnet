int Z_upplot(int *heads, int *tails, int n_edges, double **Z, double zdelta, 
	 double prior_mu, double prior_sd, int g, int k,
	 double **Znew, double *beta, int dir,
	 double *llk, int nsubsamp);
double loglike_yplot(int *heads, int *tails, int n_edges, int g, double **Z, 
		 int k, double *beta, int dir, int nsubsamp);
int beta_upplot(int *heads, int *tails, int n_edges, double **Z, int g, int k, 
	    double prior_mu, double prior_sd, 
	    double *llk, double *beta, int dir, double bdelta, int nsubsamp);

