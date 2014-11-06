 /* Calculating the mean and standard deviation for template match
 * MEX file
 * 
 * Raghav 2014
 * 
 * inputs: Data      :   This is a vector of the data for which variance has to be calculated
 *
 * outputs:outputs
 *         outputs   : This is the variance
 *          
 * version 1.0
 * this is the matlab version of the same code
 *   
 -----------------------------------*/

#include "mex.h"
#include <math.h>
#include <matrix.h>

void mexFunction(
  int nOUT, mxArray *outputs[],
  int nINP, const mxArray *inputs[])
{
  
      int i, DataLen;
      double *Data, *OutputValues;
      double VAR, MEAN;
      
      if (nINP < 1) 
        mexErrMsgTxt("1 input required");
      
      if (nOUT > 1) 
        mexErrMsgTxt("Too many output arguments");
      
      Data = (double *)mxGetPr(inputs[0]);
      
      outputs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
      OutputValues = mxGetPr(outputs[0]);
      
      DataLen = mxGetM(inputs[0]) * mxGetN(inputs[0]);
      
      MEAN = 0;
      
      for (i = 0; i < DataLen; i++)
          MEAN = MEAN + Data[i];
      MEAN = MEAN/DataLen;
     
      VAR = 0;
      for (i = 0; i < DataLen; i++)
          VAR = VAR + (Data[i] - MEAN)*(Data[i] - MEAN);
      VAR = VAR/(DataLen - 1);
      
      OutputValues[0] = VAR;
}
