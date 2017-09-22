/* Mex file for 1) Faster NLM, and 2) Computing the corresponding weights
 * Inputs - Image patches V of size [(m*n) x (2*n_rad+1)^2] , Image of size of [m x n x ndims], h, p_rad, n_rad
 * Outputs - Output filtered image, Weight matrix data (W_rows, W_cols, W_vals, W_wmax_vals, sum_W)
*/

/* Reference file for Mex creation - 
 * edit([matlabroot '/extern/examples/mex/arrayProduct.c']);
   edit([matlabroot '/extern/examples/mex/mexfunction.c']);
 */

#include "mex.h"
#include "math.h"
#ifndef min
#define min(a, b)        ((a) < (b) ? (a): (b))
#endif
#ifndef max
#define max(a, b)        ((a) > (b) ? (a): (b))
#endif

void get_V1(double *V1, double *V, int row, int N_total, int p_size)
{   
    int i;
    for (i = 0; i < p_size; i++) 
    {
        V1[i] = V[i*N_total + row];
        /*mexPrintf("%f\t", V1[i]);*/
    }
}
void get_V2(double *V2, double *V, int row, int N_total, int p_size)
{
    int i;
    for (i = 0; i < p_size; i++) 
    {
        V2[i] = V[i*N_total + row];
        /*mexPrintf("%f\t", V2[i]);*/
    }
}
double get_d(double *V1, double *V2, int p_size)
{
    int i, j;
    double sum, v1, v2;
    double d;
    sum = 0;
    for(i=0;i<1;i++)
    {
        for(j=0;j<p_size;j++)
        {
            v1 = V1[j + i];
            v2 = V2[j + i];
            sum = sum + ((v1-v2)*(v1-v2));
        }
    }
    d = sum;
    return d;
}
/* The actual function */
void NLM_fast(double *J, double *W_rows, double *W_cols, double *W_vals, double *W_wmax_vals, double *sum_W, double *V, double *I, double *I_padded, int rows, int cols, int total_dims, double h, int f, int t)
{
  int m, n, ndims, i, j, i1, j1;
  int r, s, rmin, rmax, smin, smax;
  int pixel_no, neigh_no;
  int counter, ndims_i;
  int p_size;
  double d, w, I_val, *V1, *V2;
  double wmax, sweight;
  double average[3];
  
  m = rows; n = cols; ndims = total_dims;

  p_size = (2*f+1)*(2*f+1)*(ndims);
  V1 = mxGetPr(mxCreateDoubleMatrix(1, p_size, mxREAL));
  V2 = mxGetPr(mxCreateDoubleMatrix(1, p_size, mxREAL));
 
  mexPrintf("Parameters h=%.2f, p_rad=%d, n_rad=%d, img_size=[%d,%d,%d]\n", h, f, t, m, n, ndims);
   
  counter = -1;
  for(i=0;i<m;i++)
  {
      for(j=0;j<n;j++)
      {
          i1 = i + f;
          j1 = j + f;
          wmax = 0;
          average[0]=0;average[1]=0;average[2]=0;
          sweight = 0;
          rmin = max(i1-t,f);      /*rmin = max(i1-t,f+1);   - Corrected! */
          rmax = min(i1+t,m+f-1);  /*rmax = min(i1+t,m+f);   - Corrected! */
          smin = max(j1-t,f);      /*smin = max(j1-t,f+1);   - Corrected! */
          smax = min(j1+t,n+f-1);  /*smax = min(j1+t,n+f);   - Corrected! */
          pixel_no = i + j*m;      /*pixel_no = i + (j-1)*m  - Corrected! */
          /*mexPrintf("m=%d, n=%d, p_size=%d, average=[%f, %f, %f]\n", m, n, p_size, average[0], average[1], average[2]);*/
          for(r=rmin;r<=rmax;r++)
          {
              for(s=smin;s<=smax;s++)
              {
                  /*mexPrintf("r=%d, s=%d, f=%d, m=%d, n=%d, p_size=%d\n", r, s, f, m, n, p_size);*/
                  if( (r==i1) && (s==j1) )
                  {
                      ; /*do nothing*/
                  }
                  else
                  {
                    /*mexPrintf("r=%d, s=%d, f=%d, m=%d, n=%d, p_size=%d\n", r, s, f, m, n, p_size);*/
                    counter = counter + 1;
                    neigh_no = (r-f) + (s-f)*m; /*neigh_no = (r-f) + (s-f-1)*m;  - Corrected!*/
                    get_V1(V1, V, pixel_no, m*n, p_size);
                    get_V2(V2, V, neigh_no, m*n, p_size);
                    d = get_d(V1, V2, p_size); 
                    w = exp(-d/pow(h,2));
                    /*mexPrintf("For pixel_no=%d, neigh_no=%d, d=%f\n", pixel_no, neigh_no, d);*/
                    if(w>wmax) {wmax = w;}
                    sweight = sweight + w;
                    for(ndims_i=0; ndims_i<ndims; ndims_i++)
                    {
                        I_val = I_padded[s*(m+2*f) + r + ndims_i*(m+2*f)*(n+2*f)];
                        average[ndims_i] = average[ndims_i] + w*I_val;
                    }
                    /*mexPrintf("Ival = %.2f, average=[%f, %f, %f]\n", I_val, average[0], average[1], average[2]);*/
                    W_rows[counter] = pixel_no+1; /*To use it in MATLAB, icrement this value!*/
                    W_cols[counter] = neigh_no+1; /*To use it in MATLAB, icrement this value!*/
                    W_vals[counter] = w;
                  }
              }
          }
          for(ndims_i=0; ndims_i<ndims; ndims_i++)
          {
              I_val = I_padded[j1*(m+2*f) + i1 + ndims_i*(m+2*f)*(n+2*f)];
              average[ndims_i] = average[ndims_i] + wmax*I_val;
          }   
          sweight = sweight + wmax;
          W_wmax_vals[pixel_no] = wmax;
          sum_W[pixel_no] = sweight; 
          for(ndims_i=0; ndims_i<ndims; ndims_i++)
          {
              if(sweight > 0) 
              {
                  J[j*m + i + ndims_i*m*n] = average[ndims_i]/sweight;
              }
              else 
              { 
                  J[j*m + i + ndims_i*m*n] = I[j*m + i + ndims_i*m*n];
              }  
          }
      }
  }
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
/* variable declarations here */
double *V;
double *I;
double *I_padded;
double *J;
double *W_rows, *W_cols, *W_vals, *W_wmax_vals, *sum_W;
double h;
int  p_rad, n_rad;
double d, *V1, *V2;
int rows, cols, total_dims, isRGB; 
const mwSize *Idims_pr;


/* check for the correct no. of inputs */
if(nrhs != 6) 
{
    mexErrMsgIdAndTxt("computeNLMweights:nrhs", "Six inputs required.");
}

/* Check whether the input image is RGB */
Idims_pr = mxGetDimensions(prhs[1]);
total_dims = Idims_pr[2];
if(total_dims != 3) {total_dims = 1;}
rows = Idims_pr[0];
cols = Idims_pr[1];
    

/* Verify that V is double */
if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) 
{
    mexErrMsgIdAndTxt("computeNLMweights:notDouble", "Input matrix V must be of type double.");
}
/* Verify that it is of the correct size */
if ( !(mxGetM(prhs[0]) == (rows*cols)) )
{
    mexPrintf("\n m_V=%d, m=%d, n=%d", mxGetM(prhs[0]), rows, cols); 
    mexErrMsgIdAndTxt("computeNLMweights:VsizeIncorrect", "Input matrix V first dimesnion must be of size (m*n).");      
}
if ( !(mxGetN(prhs[0]) == ((2*mxGetScalar(prhs[4]) + 1)*(2*mxGetScalar(prhs[4]) + 1)*total_dims) ) )
{
    mexPrintf("total_dims=%d", total_dims);
    mexErrMsgIdAndTxt("computeNLMweights:VsizeIncorrect", "Input matrix V second dimesnion must be of size (2*p_rad+1)^2*(total_dims).");
}
/* Verify that I is double */
if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) 
{
    mexErrMsgIdAndTxt("computeNLMweights:notDouble", "Input matrix I must be of type double.");
}
/* Verify that all other arguments are scalar */
if( !mxIsDouble(prhs[3]) || !mxIsDouble(prhs[4]) || !mxIsDouble(prhs[5]) ||
     mxIsComplex(prhs[3]) || mxIsComplex(prhs[4]) || mxIsComplex(prhs[5]) ||
     (mxGetNumberOfElements(prhs[3]) != 1) || (mxGetNumberOfElements(prhs[4]) != 1) || (mxGetNumberOfElements(prhs[5]) != 1) ) 
{
    mexErrMsgIdAndTxt("computeNLMweights:notScalar", "NLM parameters must be scalar.");
}
        
/* Get the inputs */
V = mxGetPr(prhs[0]);
I = mxGetPr(prhs[1]);
I_padded = mxGetPr(prhs[2]);
h = mxGetScalar(prhs[3]);
p_rad = (int) mxGetScalar(prhs[4]);
n_rad = (int) mxGetScalar(prhs[5]);

/* Prepare the outputs */
if(total_dims == 1) 
{
    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
}
else if(total_dims == 3)
{
    plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[1]),mxGetDimensions(prhs[1]), mxDOUBLE_CLASS, mxREAL);
}
plhs[1] = mxCreateDoubleMatrix(1, rows*cols*pow(2*n_rad+1,2), mxREAL);
plhs[2] = mxCreateDoubleMatrix(1, rows*cols*pow(2*n_rad+1,2), mxREAL);
plhs[3] = mxCreateDoubleMatrix(1, rows*cols*pow(2*n_rad+1,2), mxREAL);
plhs[4] = mxCreateDoubleMatrix(1, rows*cols, mxREAL);
plhs[5] = mxCreateDoubleMatrix(1, rows*cols, mxREAL);
J = mxGetPr(plhs[0]);
W_rows = mxGetPr(plhs[1]);
W_cols = mxGetPr(plhs[2]);
W_vals = mxGetPr(plhs[3]);
W_wmax_vals = mxGetPr(plhs[4]);
sum_W       = mxGetPr(plhs[5]);

/* code here */
NLM_fast(J, W_rows, W_cols, W_vals, W_wmax_vals, sum_W, V, I, I_padded, rows, cols, total_dims, h, p_rad, n_rad);
}
