 /* Calculating the correlation between two spike trains or two vectors
 * MEX file
 * 
 * Raghav 2014
 * 
 * inputs: Data1, Data2 : The two vectors for which correlation has to be calculated
 *
 * outputs:outputs
 *         outputs   : This is the correlation
 *          
 * version 1.0
 * this is the matlab version of the same code
 * Correlation is (Data1 - mean(Data1)) * (Data2 - mean(Data2))'/(norm(Data1 - mean(Data1)) * norm(Data2 - mean(Data2)))
 -----------------------------------*/

#include "mex.h"
#include <math.h>
#include <matrix.h>

void mexFunction(
  int nOUT, mxArray *outputs[],
  int nINP, const mxArray *inputs[])
{
  
      int i, Data1Len, Data2Len;
      double *Data1, *Data2, *OutputValues;
      double Data1MEAN, Data2MEAN, VectorProduct, Data1Norm, Data2Norm;
      
      if (nINP < 2) 
        mexErrMsgTxt("1 input required");
      
      if (nOUT > 1) 
        mexErrMsgTxt("Too many output arguments");
      
      Data1 = (double *)mxGetPr(inputs[0]);
      Data2 = (double *)mxGetPr(inputs[1]);
      
      outputs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
      OutputValues = mxGetPr(outputs[0]);
      
      
      Data1Len = mxGetM(inputs[0]) * mxGetN(inputs[0]);
      Data2Len = mxGetM(inputs[0]) * mxGetN(inputs[0]);
      
      if (Data1Len != Data2Len)
          mexErrMsgTxt("Both input vectors should be of same length");
      
      Data1MEAN = 0;
      Data2MEAN = 0;
      
      // First calculate means for both sets of data
      /*for (i = 0; i < Data1Len; i++)
      {
          Data1MEAN = Data1MEAN + Data1[i];
          Data2MEAN = Data2MEAN + Data2[i];
      }
      
      Data1MEAN = Data1MEAN/Data1Len;
      Data2MEAN = Data2MEAN/Data2Len;

      // Now mean subtract both data sets
      for (i = 0; i < Data1Len; i++)
      {
          Data1[i] = Data1[i] - Data1MEAN;
          Data2[i] = Data2[i] - Data2MEAN;
      }
      */
      
      // Now multiply two vectors - this is numerator of correlation
      VectorProduct - 0;
      for (i = 0; i < Data1Len; i++)
          VectorProduct = VectorProduct + Data1[i]*Data2[i];
      
      // Now to calculate norm of the two vectors
      Data1Norm = 0;
      Data2Norm = 0;
      
      for (i = 0; i < Data1Len; i++)
      {
          Data1Norm = Data1Norm + (Data1[i]*Data1[i]);
          Data2Norm = Data2Norm + (Data2[i]*Data2[i]);
      }
      
      // printf("%g\t%g\t%g\n", VectorProduct, Data1Norm, Data2Norm);
      Data1Norm = sqrt(Data1Norm);
      Data2Norm = sqrt(Data2Norm);
      
      OutputValues[0] = VectorProduct/(Data1Norm * Data2Norm);
}
