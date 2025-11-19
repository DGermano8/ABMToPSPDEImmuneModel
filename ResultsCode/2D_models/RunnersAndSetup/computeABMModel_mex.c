#include "mex.h"
#include <math.h>    /* For round() */
#include <stdlib.h>  /* For rand(), srand() */
#include <time.h>    /* For time() */
#include <stdbool.h> /* For bool, true, false (C99) */

// mex -largeArrayDims COMPFLAGS="$COMPFLAGS -O3 -fopenmp" RunnersAndSetup/computeABMModel_mex.c

/* Helper macros for min/max */
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

/* Forward declaration of the core computational routine */
void computeABMModel_c(
    double* walker_positions, double* walker_activation, double* DCLingerTime,
    const double* C, const double* Ic, const double* BoundaryDC,
    mwSize num_walkers, double dx, double dy, mwSize Nx, mwSize Ny,
    double p_move, double C_chi, double activatedAge,
    double P_A, double P_D, double dt, double da, int NumberOfDCs,
    mwSize grid_rows
);

/*=================================================================
 * MEX GATEWAY FUNCTION
 * This is the entry point that MATLAB calls.
 *=================================================================*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* --- Argument Checking --- */
    /* Expects 6 arrays and 12 scalar parameters */
    if (nrhs != 19) {
        mexErrMsgIdAndTxt("ABM:compute_C:nrhs", "19 input arguments are required.");
    }
    if (nlhs != 4) {
        mexErrMsgIdAndTxt("ABM:compute_C:nlhs", "4 output arguments are required.");
    }
    
    /* --- Seed random number generator (only on first call) --- */
    static bool is_seeded = false;
    if (!is_seeded) {
        srand((unsigned int)time(NULL));
        is_seeded = true;
    }
    
    /* --- Get pointers to input arrays (inputs 0-5) --- */
    double* C_in = mxGetPr(prhs[2]);
    double* Ic_in = mxGetPr(prhs[4]);
    double* BoundaryDC_in = mxGetPr(prhs[5]);

    /* --- Get scalar parameters (inputs 6-17) --- */
    mwSize num_walkers = (mwSize)mxGetScalar(prhs[6]);
    double dx = mxGetScalar(prhs[7]);
    double dy = mxGetScalar(prhs[8]);
    mwSize Nx = (mwSize)mxGetScalar(prhs[9]);
    mwSize Ny = (mwSize)mxGetScalar(prhs[10]);
    double p_move = mxGetScalar(prhs[11]);
    double C_chi = mxGetScalar(prhs[12]);
    double activatedAge = mxGetScalar(prhs[13]);
    double P_A = mxGetScalar(prhs[14]);
    double P_D = mxGetScalar(prhs[15]);
    double dt = mxGetScalar(prhs[16]);
    double da = mxGetScalar(prhs[17]);
    int NumberOfDCs = (int)mxGetScalar(prhs[18]);
    
    /* --- Basic Input Validation --- */
    if (mxGetM(prhs[0]) != num_walkers || mxGetN(prhs[0]) != 2) {
        mexErrMsgIdAndTxt("ABM:compute_C:badSize", "'walker_positions' must be a num_walkers-by-2 matrix.");
    }
    if (mxGetM(prhs[4]) != Nx || mxGetN(prhs[4]) != Ny) {
        mexErrMsgIdAndTxt("ABM:compute_C:badSize", "'Ic' and other grid matrices must be Nx-by-Ny.");
    }
    
    /* --- Create output arrays by duplicating inputs --- */
    plhs[0] = mxDuplicateArray(prhs[0]);
    plhs[1] = mxDuplicateArray(prhs[1]);
    plhs[2] = mxDuplicateArray(prhs[2]); /* C is passed through */
    plhs[3] = mxDuplicateArray(prhs[3]);
    
    /* --- Get pointers to output data --- */
    double* walker_pos_out = mxGetPr(plhs[0]);
    double* walker_act_out = mxGetPr(plhs[1]);
    double* DCLT_out = mxGetPr(plhs[3]);
    
    mwSize grid_rows = mxGetM(prhs[4]); /* Should be equal to Nx */
    // mwSize DCLT_size = mxGetNumberOfElements(prhs[3]);
    
    /* --- Call the core computational routine --- */
    computeABMModel_c(
        walker_pos_out, walker_act_out, DCLT_out,
        C_in, Ic_in, BoundaryDC_in,
        num_walkers, dx, dy, Nx, Ny,
        p_move, C_chi, activatedAge, P_A, P_D, dt, da, NumberOfDCs,
        grid_rows
    );
}

/*=================================================================
 * CORE COMPUTATIONAL ROUTINE
 * Contains the main simulation logic.
 *=================================================================*/
void computeABMModel_c(
    double* walker_positions, double* walker_activation, double* DCLingerTime,
    const double* C, const double* Ic, const double* BoundaryDC,
    mwSize num_walkers, double dx, double dy, mwSize Nx, mwSize Ny,
    double p_move, double C_chi, double activatedAge,
    double P_A, double P_D, double dt, double da, int NumberOfDCs,
    mwSize grid_rows
) {
    // DCLingerTime array is size 
    /* Loop over each walker */
    for (mwSize ii = 0; ii < num_walkers; ++ii) {
        /* Get current position and convert to 1-based grid indices */
        double wx = walker_positions[ii];
        double wy = walker_positions[ii + num_walkers]; /* Column-major access */
        long id_x = (long)round(wx / dx) + 1;
        long id_y = (long)round(wy / dy) + 1;
        
        id_x = MAX(1L, MIN(id_x, (long)Nx));
        id_y = MAX(1L, MIN(id_y, (long)Ny));
        
        /* --- 1. Random Motion --- */
        {
            long id_x_p1 = MIN((long)Nx, id_x + 1);
            long id_x_m1 = MAX(1L, id_x - 1);
            long id_y_p1 = MIN((long)Ny, id_y + 1);
            long id_y_m1 = MAX(1L, id_y - 1);

            /* Linear 0-based index for column-major arrays: (row-1) + num_rows * (col-1) */
            long left_idx  = (id_y - 1) + grid_rows * (id_x_m1 - 1);
            long right_idx = (id_y - 1) + grid_rows * (id_x_p1 - 1);
            long down_idx  = (id_y_m1 - 1) + grid_rows * (id_x - 1);
            long up_idx    = (id_y_p1 - 1) + grid_rows * (id_x - 1);
            
            double prob_left  = MAX(0.0, (1.0 - Ic[left_idx]) * (p_move / 4.0));
            double prob_right = MAX(0.0, (1.0 - Ic[right_idx]) * (p_move / 4.0));
            double prob_down  = MAX(0.0, (1.0 - Ic[down_idx]) * (p_move / 4.0));
            double prob_up    = MAX(0.0, (1.0 - Ic[up_idx]) * (p_move / 4.0));
            double prob_stay  = 1.0 - p_move;
            
            double total_prob = prob_left + prob_right + prob_down + prob_up + prob_stay;
            if (total_prob > 1e-9) { /* Normalize probabilities */
                prob_left /= total_prob;
                prob_right /= total_prob;
                prob_down /= total_prob;
                prob_up /= total_prob;
            }

            double p1 = prob_left;
            double p2 = p1 + prob_right;
            double p3 = p2 + prob_down;
            double p4 = p3 + prob_up;
            
            double Rad = (double)rand() / (double)RAND_MAX;

            if (Rad < p1)         id_x = id_x_m1;
            else if (Rad < p2)    id_x = id_x_p1;
            else if (Rad < p3)    id_y = id_y_m1;
            else if (Rad < p4)    id_y = id_y_p1;
        }

        /* --- 2. Chemotaxis Motion --- */
        {
            long id_x_p1 = MIN((long)Nx, id_x + 1);
            long id_x_m1 = MAX(1L, id_x - 1);
            long id_y_p1 = MIN((long)Ny, id_y + 1);
            long id_y_m1 = MAX(1L, id_y - 1);
            
            long here_idx  = (id_y - 1) + grid_rows * (id_x - 1);
            long left_idx  = (id_y - 1) + grid_rows * (id_x_m1 - 1);
            long right_idx = (id_y - 1) + grid_rows * (id_x_p1 - 1);
            long down_idx  = (id_y_m1 - 1) + grid_rows * (id_x - 1);
            long up_idx    = (id_y_p1 - 1) + grid_rows * (id_x - 1);
            
            double C_Here = C[here_idx];
            double prob_left  = MAX(0.0, (1.0 - Ic[left_idx]) * (C[left_idx] - C_Here) * (C_chi / 4.0));
            double prob_right = MAX(0.0, (1.0 - Ic[right_idx]) * (C[right_idx] - C_Here) * (C_chi / 4.0));
            double prob_down  = MAX(0.0, (1.0 - Ic[down_idx]) * (C[down_idx] - C_Here) * (C_chi / 4.0));
            double prob_up    = MAX(0.0, (1.0 - Ic[up_idx]) * (C[up_idx] - C_Here) * (C_chi / 4.0));
            double prob_stay  = MAX(0.0, 1.0 - (prob_left + prob_right + prob_down + prob_up));
            
            double total_prob = prob_left + prob_right + prob_down + prob_up + prob_stay;
            if (total_prob > 1e-9) { /* Normalize as in original code */
                prob_left /= total_prob;
                prob_right /= total_prob;
                prob_down /= total_prob;
                prob_up /= total_prob;
            }

            double p1 = prob_left;
            double p2 = p1 + prob_right;
            double p3 = p2 + prob_down;
            double p4 = p3 + prob_up;

            double Rad = (double)rand() / (double)RAND_MAX;

            if (Rad < p1)         id_x = id_x_m1;
            else if (Rad < p2)    id_x = id_x_p1;
            else if (Rad < p3)    id_y = id_y_m1;
            else if (Rad < p4)    id_y = id_y_p1;
        }
        
        /* --- 3. Update Position --- */
        walker_positions[ii] = dx * (id_x - 1);
        walker_positions[ii + num_walkers] = dy * (id_y - 1);

        /* --- 4. Update Activation --- */
        {
            bool states_can_change = !(walker_activation[ii] >= (activatedAge - 1e-9) );
            states_can_change = true;
            
            long here_idx = (id_y - 1) + grid_rows * (id_x - 1);
            int i_DC_id = BoundaryDC[here_idx];

            // bool near_boundary_id = (i_DC_id > 0.0 && states_can_change);
            bool near_boundary = (i_DC_id > 0.0);
            
            double prob_up = (states_can_change && near_boundary) ? (P_A * (1.0 - (walker_activation[ii]+da)/ activatedAge)) : 0.0;
            double prob_down = (states_can_change && !near_boundary) ? (P_D * ((walker_activation[ii]-da) / activatedAge)) : 0.0;
            
            if (near_boundary) {
                /* Compute the correct linear index for column-major order */
                mwSize dcl_index = (i_DC_id - 1) * num_walkers + (ii);
                if (dcl_index >= 0 && dcl_index < (num_walkers * NumberOfDCs)) {
                    DCLingerTime[dcl_index] += dt;
                }
            }
            
            double Rad = (double)rand() / (double)RAND_MAX;
            if (Rad < prob_up) {
                walker_activation[ii] += da;
            }
            else if (Rad > prob_up && Rad < prob_up + prob_down) {
                walker_activation[ii] -= da;
            }
            
            if (fabs(walker_activation[ii]) < da && Rad < prob_up)
            {
                for (int dc = 0; dc < NumberOfDCs; ++dc) {
                    mwSize dcl_index = dc * num_walkers + (ii);
                    DCLingerTime[dcl_index] = 0.0;
                }
            }

        }
    }
}