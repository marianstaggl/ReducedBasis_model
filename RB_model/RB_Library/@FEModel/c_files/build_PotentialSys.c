/*==========================================================
 * c code to create a mex function for the assembly of the 
 * heat transfer equation due to a finite element approach
 *
 * Inputs -->
 * elems: array containing the connectivity of the elements
 * dbdx: derivations of basefunctions along x (needed for mapping)
 * dbdy: derivations of basefunctions along y (needed for mapping)
 * weights: weights of the mappend integration points
 *
 * Outputs -->
 * k_idx: sparse systems matrix
 * j_idx: residuals vector
 * t_val: residuals vector
 *
 *
 * The calling syntax is:
 *
 * [k_idx, j_idx, t_val] = build_HeatTransf(elems, dbdx, dbdy, weights)
 *
 *========================================================*/

#include "mex.h"


/* The computational routine */
void build_PotentialSys(int n_base, int n_elem, int n_weight, 
        double *elems, double *dbdx, double *dbdy, double *weights, 
        double *k_idx, double *j_idx, double *t_val)
{
    /* declare the running variables */
    int i; int j; int k; int m; int count;
    
    /* loop over all the elements */
    count = 0;
    for (i=0; i<n_elem; i++) {
        
        /* loop over the base functions */
        for (k=0; k<n_base; k++) {
            for (j=0; j<n_base; j++) {
                /* set the global indices */
                j_idx[count] = elems[i+j*n_elem];
                k_idx[count] = elems[i+k*n_elem];
                
                /* integrate the residuals function */
                for (m=0; m<n_weight; m++){
                    t_val[count] = t_val[count] + weights[i+m*n_elem]*
                          (dbdx[j+(n_base*m)+i*(n_base*n_weight)]*
                           dbdx[k+(n_base*m)+i*(n_base*n_weight)]+ 
                           dbdy[j+(n_base*m)+i*(n_base*n_weight)]*
                           dbdy[k+(n_base*m)+i*(n_base*n_weight)]);
                }
                
                /* count one up */
                count++;
            }
        }
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /* define the variables for the functions input */
    int n_base; int n_elem; int n_weight;
    double *elems; double *dbdx; double *dbdy; double *weights;
    
    /* define the variables for the functions output */
    double *k_idx; double *j_idx; double *t_val;
    
    /* link the real values of the inputs */
    elems = mxGetDoubles(prhs[1]);
    n_base = mxGetN(prhs[1]); 
    n_elem = mxGetM(prhs[1]); 
    dbdx = mxGetDoubles(prhs[2]);
    dbdy = mxGetDoubles(prhs[3]);
    weights = mxGetDoubles(prhs[4]);
    n_weight = mxGetN(prhs[4]);

    /* preallocate and link vectors for the sparse system assembly */
    plhs[0] = mxCreateDoubleMatrix(1, n_elem*n_base*n_base, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, n_elem*n_base*n_base, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, n_elem*n_base*n_base, mxREAL);
    
    /* run the function to build the matrix */
    k_idx = mxGetDoubles(plhs[0]);
    j_idx = mxGetDoubles(plhs[1]);
    t_val = mxGetDoubles(plhs[2]);
    
    build_PotentialSys(n_base, n_elem, n_weight, elems, dbdx, dbdy, 
            weights, k_idx, j_idx, t_val);
}
