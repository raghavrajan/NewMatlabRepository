 /* Calculating the template match
 * MEX file
 * 
 * Raghav 2009
 * 
 * inputs: m x n array 
 *          
 *       
 * outputs:outputs
 *          
 *         outputs: this is an 1 X (n1 - n + 1) array with indices indicating the matches for each location
 *          
 * version 1.0
 -----------------------------------*/

#include "mex.h"
#include <math.h>
#include <matrix.h>

void mexFunction(
  int nOUT, mxArray *outputs[],
  int nINP, const mxArray *inputs[])
{
  
      double *Data, *Output;
      int DataRows, DataCols, i, j, Cols, Rows, ArrayCounter, StartPoint;
      double Temp[100000];
      
      if (nINP < 3) 
        mexErrMsgTxt("3 inputs required");
      
      if (nOUT > 1) 
        mexErrMsgTxt("Too many output arguments");
      
      Data = (double *) mxGetPr(inputs[0]);
      Rows = mxGetScalar(inputs[1]);
      Cols = mxGetScalar(inputs[2]);
      StartPoint = mxGetScalar(inputs[3]);
      DataRows = mxGetM(inputs[0]);
      DataCols = mxGetN(inputs[0]);
      
      ArrayCounter = 0;

      for (i = 0; i < Cols*Rows; i++)
      {
              Temp[ArrayCounter] = Data[ArrayCounter + StartPoint*DataRows];
              printf("i is %i, Temp is %g\n", i, Temp[ArrayCounter]);
              ArrayCounter++;
      }
      
      /*      for (i = 0; i < (0 + Cols); i++)
      {
          for (j = 0; j < Rows; j++)
          {
              Temp[ArrayCounter] = Data[i*DataRows + j];
              ArrayCounter++;
              printf("i is %i, j is %i, counter is %i, and Temp is %g\n", i, j, ((i-0)*DataRows + j), Temp[(i-0)*DataRows + j]);
          }
          
      }*/
              
      outputs[0] = mxCreateDoubleMatrix(Rows, Cols, mxREAL);
      Output = mxGetPr(outputs[0]);

      for(i = 0; i < Rows*Cols; i++)
      {
          Output[i] = Temp[i];
      }
}
