 /* Calculating the mean and standard deviation for template match
 * MEX file
 * 
 * Raghav 2009
 * 
 * inputs: Data      :   This is a vector of the data which has to be windowed and mean and std for each window has to be calculated
 *         WindowLen :   This is a scalar with the length of each window of data to be considered
 *         WinNo     :   This is a scalar with the total number of windows to be considered
 *         StepSize  :   This is a scalar with the step size between the start of two consecutive windows
 * outputs:outputs
 *         outputs: this is an array with mean values in the first half and corresponding std values in the next half
 *          
 * version 1.0
 * this is the matlab version of the same code
 *   WinMean = zeros((size(TempS,2) - size(WMotif,2) + 1), 1);
                WinSTD = zeros((size(TempS,2) - size(WMotif,2) + 1), 1);
                
                for ColNo = 1:(size(TempS,2) - size(WMotif, 2) + 1),
                    StartIndex = ((ColNo - 1)*size(WMotif,1)) + 1;
                    WinIndices = StartIndex:1:(StartIndex + size(WMotif,1)*size(WMotif,2) - 1);
                    WinMean(ColNo) = mean(TempS(WinIndices));
                    WinSTD(ColNo) = std(TempS(WinIndices));
%                        Data = (Data - mean(mean(Data)))./std(reshape(Data, 1, (size(Data,1) * size(Data,2))));
%                       Match1(ColNo,1) = 1/sum(sum((abs(Data - WMotif))));
                end
 -----------------------------------*/

#include "mex.h"
#include <math.h>
#include <matrix.h>

void mexFunction(
  int nOUT, mxArray *outputs[],
  int nINP, const mxArray *inputs[])
{
  
      double *Data, *OutputValues;
      int WinLen, WinNo, StepSize, StartIndex, i, j, DataLen;
      double Mean, STD;
      
      if (nINP < 4) 
        mexErrMsgTxt("4 inputs required");
      
      if (nOUT > 1) 
        mexErrMsgTxt("Too many output arguments");
      
      Data = (double *)mxGetPr(inputs[0]);
      WinLen = mxGetScalar(inputs[1]);
      WinNo = mxGetScalar(inputs[2]);
      StepSize = mxGetScalar(inputs[3]);
      
      outputs[0] = mxCreateDoubleMatrix(2*WinNo, 1, mxREAL);
      OutputValues = mxGetPr(outputs[0]);
      
      DataLen = mxGetM(inputs[0]) * mxGetN(inputs[0]);
      
      /*printf("Data length is %i and Window length is %i and number of windows is %i and step size is %i\n", DataLen, WinLen, WinNo, StepSize);*/
      
      for (i = 0; i < WinNo; i++)
      {
          Mean = 0;
          STD = 0;
          StartIndex = i*StepSize;
          
          for (j = StartIndex; j <= (StartIndex + WinLen - 1); j++)
          {
              Mean = Mean + Data[j];
          }
          Mean = Mean/WinLen;
          
          for (j = StartIndex; j <= (StartIndex + WinLen - 1); j++)
          {
              STD = STD + (Data[j] - Mean)*(Data[j] - Mean);
          }
          STD = sqrt(STD/(WinLen - 1));
                
          OutputValues[i] = Mean;
          OutputValues[WinNo + i] = STD;
      }
}
