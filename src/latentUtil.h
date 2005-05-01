double *dvector(int n);
int *ivector(int n);
double **dmatrix(int n,int m);
int **imatrix(int n,int m);
void free_dmatrix(double **A, int n);
void free_imatrix(int **A, int n);
void free_dvector(double *x);
void free_ivector(int *x);
void print_dvector (double *a, int length, FILE *stream );
void print_drvector(double *a, int n, FILE *stream);
void print_ivector (int *a, int length, FILE *stream );
void print_dmatrix (double **a, int nrow, int ncol, FILE* stream);
void print_imatrix (int **a, int nrow, int ncol, FILE* stream);
void init_dvector(double *x, int n, double value);
void init_dmatrix(double **A, int n, int m, double value);
double *cat_dmectors(double *x, int nx, double *y, int ny);
double *cat_dmector_scaler(double *x, int nx, double y, int end);
double *dvector_times_matrix(double *x, int n,double **A, int m, double *b);
void dscaler_times_matrix(double x, double **A, int n, int m, double **B);
void dmatrix_multiply(double **A,int na,int ma, double **B, int mb, 
		      double **C);
void imatrix_multiply(int **A,int na,int ma, int **B, int mb, int **C);
void dmatrix_addition(double **A, int n, int m, double **B);
void init_imatrix(int **A, int n, int m, int value);
void t(double **A, int n, int m, double **tA);
void copy_dmatrix(double **source,double **dest,int n,int m);
double mean(double *x, int n);
void svd_symm(double **A, int n, double **eAl, double **eAleftVectors, 
	     double **eArightVectors/*, FILE *fout*/);
/* these two thanks to Rahpael Gottardo*/
void qr_solve(double **x, int *n1, double ** y, double **coef, int singular);
int inverse(double **mat1, int *n ,double **res);
int sym_eigen(double **A, int n, int vectors, double *EValues, double **EVectors);
/*  FORTRAN function declarations below */
void   F77_CALL(dqrdc2)(double *xt, int *n1, int *n2, int *p, double *tol, 
			int *rank, double *qraux, int *pivot, double *work);
void   F77_CALL(dqrcf)(double *xt, int *n1, int *rank, double *qraux, 
		       double *yt, int *n2, double *coeft, int *info);
void   F77_NAME(rs)(int *n1, int *n2, double *vA, double *l, int *vectorsflag, 
		    double *vEVectors, double *fv1, double *fv2, int *err);





