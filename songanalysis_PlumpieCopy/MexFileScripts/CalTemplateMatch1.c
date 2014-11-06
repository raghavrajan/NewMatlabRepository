 /* Calculating the template match
 * MEX file
 * 
 * Raghav 2009
 * 
 * inputs: Template :   This is an m x n array of the spectrogram of the template
 *         Data     :   This is an m x n1 array of the spectrogram of the file in which matches to template have to be found 
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
  
      double *Data, *Template, *WinMean, *WinSTD, *OutputMatchValues;
      int TempLen, TempRows, TempCols, DataRows, DataCols, i, j, k, ArrayCounter;
      double Match, Sum;
      double Window[1000000];
      
      if (nINP < 4) 
        mexErrMsgTxt("4 inputs required");
      
      if (nOUT > 1) 
        mexErrMsgTxt("Too many output arguments");
      
      Template = (double *)mxGetPr(inputs[0]);
      Data = (double *) mxGetPr(inputs[1]);
      WinMean = (double *) mxGetPr(inputs[2]);
      WinSTD = (double *) mxGetPr(inputs[3]);
      
      TempLen = mxGetM(inputs[0]) * mxGetN(inputs[0]);
      TempRows = mxGetM(inputs[0]);
      TempCols = mxGetN(inputs[0]);
      
      DataRows = mxGetM(inputs[1]);
      DataCols = mxGetN(inputs[1]);
      
      /*printf("Data # of rows and cols is %i, %i, and template # of rows and cols is %i, %i\n", DataRows, DataCols, TempRows, TempCols);
      printf("First Data value is %g\n", Data[0]);*/
      
      outputs[0] = mxCreateDoubleMatrix((DataCols - TempCols + 1), 1, mxREAL);
      OutputMatchValues = mxGetPr(outputs[0]);
      
      for(i = 0; i < (DataCols - TempCols + 1); i++)
      {
          Match = 0;
          for(j = 0; j < (TempCols*TempRows); j++)
          {
              Match += (((Data[j + i*DataRows] - WinMean[i])/WinSTD[i]) - Template[j]) * (((Data[j + i*DataRows] - WinMean[i])/WinSTD[i]) - Template[j]);
          }
          OutputMatchValues[i] = 1/Match;
      }
      

      
/*      for(i = 0; i < (DataCols - TempCols - 10); i++)
      {
          OutputMatchValues[i] = 0;
      }

/*      Now set index values for range 1-100 with index indicating number of decibels from maximum */
/*      for(i = 0; i < (DataCols - TempCols - 10); i++)
      {
          OutputMatchValues[i] = MatchValue[i];
      }*/
}
