/*
 * computePhenotypeModel_Opt_mex.c
 * Optimized MEX function for one PDE time-step with masking and OpenMP
 // mex -largeArrayDims COMPFLAGS="$COMPFLAGS -O3 -fopenmp" RunnersAndSetup/computePhenotypeModelTime_mex.c
 */
#include "mex.h"
#include <stdbool.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Corrected the usage string to match the expected number of inputs (20)
    if (nrhs != 20 || nlhs != 1) {
        mexErrMsgIdAndTxt("pde_step_masked_opt_mex:nargs", "Usage: u_new = pde_step_masked_opt_mex(u,C,A,A_Ic,chix,Dx,dx,dy,rho_p,rho_m,amax,da,dt,i_e,i_w,i_n,i_s,i_c,dims,num_t);");
    }

    /* Extract pointers and scalars */
    double *u_in     = mxGetPr(prhs[0]);
    double *C        = mxGetPr(prhs[1]);
    double *A        = mxGetPr(prhs[2]);
    double *A_Ic     = mxGetPr(prhs[3]);
    double chix      = mxGetScalar(prhs[4]);
    double Dx        = mxGetScalar(prhs[5]);
    double inv_dx    = mxGetScalar(prhs[6]);
    double inv_dy    = mxGetScalar(prhs[7]);
    double rho_p_hat = mxGetScalar(prhs[8]);
    double rho_m_hat = mxGetScalar(prhs[9]);
    double amax      = mxGetScalar(prhs[10]);
    double inv_da    = mxGetScalar(prhs[11]);
    double dt        = mxGetScalar(prhs[12]);

    mxLogical *i_e = mxGetLogicals(prhs[13]);
    mxLogical *i_w = mxGetLogicals(prhs[14]);
    mxLogical *i_n = mxGetLogicals(prhs[15]);
    mxLogical *i_s = mxGetLogicals(prhs[16]);
    mxLogical *i_c = mxGetLogicals(prhs[17]);

    int32_T *dims_int = (int32_T*)mxGetData(prhs[18]);
    
    // FIX 1: Correctly get the scalar value for num_T
    mwSize num_T = (mwSize)mxGetScalar(prhs[19]);

    mwSize nx = dims_int[0], ny = dims_int[1], na = dims_int[2];
    mwSize slice = nx * ny, total = slice * na;

    /* Allocate memory for ping-pong buffers */
    double *u = mxMalloc(total * sizeof(double));
    double *u_new = mxMalloc(total * sizeof(double));
    memcpy(u, u_in, total * sizeof(double));

    /* Hoist invariants */
    double half_chi = 0.5 * chix;
    double inv_dx_sqrd = inv_dx * inv_dx;
    double inv_dy_sqrd = inv_dy * inv_dy;

    /* Time-stepping loop */
    for (mwSize time_iter = 0; time_iter < num_T; ++time_iter) {
        #ifdef _OPENMP
        #pragma omp parallel for collapse(3)
        #endif
        for (mwSize k = 0; k < na; ++k) {
            for (mwSize j = 0; j < ny; ++j) {
                for (mwSize i = 0; i < nx; ++i) {
                    mwSize idx = k * slice + j * nx + i;
                    double un = 0.0;

                    if (i_c[idx]) {
                        double u_c = u[idx];
                        double C_c = C[idx];

                        /* X flux */
                        double flux_e = 0, flux_w = 0;
                        if (i < nx - 1 && i_e[idx]) {
                            double u_e = u[idx + 1];
                            double diff_e = Dx * (u_e - u_c);
                            double chem_e = -half_chi * (u_e + u_c) * (C[idx + 1] - C_c);
                            flux_e = (chem_e + diff_e);
                        }
                        if (i > 0 && i_w[idx]) {
                            double u_w = u[idx - 1];
                            double diff_w = Dx * (u_c - u_w);
                            double chem_w = -half_chi * (u_c + u_w) * (C_c - C[idx - 1]);
                            flux_w = (chem_w + diff_w);
                        }
                        double grad_x = (flux_e - flux_w) * inv_dx_sqrd;

                        /* Y flux */
                        double flux_n = 0, flux_s = 0;
                        if (j < ny - 1 && i_n[idx]) {
                            double u_n = u[idx + nx];
                            double diff_n = Dx * (u_n - u_c);
                            double chem_n = -half_chi * (u_n + u_c) * (C[idx + nx] - C_c);
                            flux_n = (chem_n + diff_n);
                        }
                        if (j > 0 && i_s[idx]) {
                            double u_s = u[idx - nx];
                            double diff_s = Dx * (u_c - u_s);
                            double chem_s = -half_chi * (u_c + u_s) * (C_c - C[idx - nx]);
                            flux_s = (chem_s + diff_s);
                        }
                        double grad_y = (flux_n - flux_s) * inv_dy_sqrd;

                        /* Phenotypic advection */
                        double flux_ap = 0, flux_am = 0;
                        double A_Ic_c = A_Ic[idx], A_c = A[idx];
                        if (k > 0) {
                            double amh = 0.5 * (A_c + A[idx - slice]);
                            double fmh = rho_p_hat * (amax - amh) * A_Ic_c - rho_m_hat * (amh) * (1 - A_Ic_c);
                            double um = (fmh <= 0 ? u_c : u[idx - slice]);
                            flux_am = fmh * um;
                        }
                        if (k < na - 1) {
                            double aph = 0.5 * (A[idx + slice] + A_c);
                            double fph = rho_p_hat * (amax - aph) * A_Ic_c - rho_m_hat * (aph) * (1 - A_Ic_c);
                            double up = (fph <= 0 ? u[idx + slice] : u_c);
                            flux_ap = fph * up;
                        }
                        double grad_a = -(flux_ap - flux_am) * inv_da;

                        un = u_c + dt * (grad_x + grad_y + grad_a);
                    }
                    u_new[idx] = un;
                }
            }
        }

        /* Ping-pong buffer swap */
        double *temp = u;
        u = u_new;
        u_new = temp;
    }

    // FIX 2: Create output array AFTER the loop and copy final result
    mwSize dims[3] = {nx, ny, na};
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double *out_ptr = mxGetPr(plhs[0]);
    memcpy(out_ptr, u, total * sizeof(double));

    // Free the unused buffer to prevent memory leaks
    mxFree(u_new);
}