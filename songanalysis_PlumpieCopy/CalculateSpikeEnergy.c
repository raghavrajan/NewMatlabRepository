 /* find_spikes
 * Finding spikes in multi-unit data
 * MEX file
 * 
 * Raghav 2006
 * 
 * The program finds spikes using the following algorithm.
 * The algorithm follows one used by DeAngelis and Newsome (J.Neurosci. Feb 15, 1999)
 *
 * If the voltage crosses an upper threshold and if it crosses a lower threshold within a 
 * specific time-window following the upper threshold crossing, it is considered a spike.
 *
 * Similarly, if the voltage crosses a lower threshold and if it crosses an upper threshold within a 
 * specific time-window following the upper threshold crossing, it is considered a spike.
 *
 * The point of the first threshold crossing is considered as the time of the spike and the program
 * jumps one window_size to avoid picking up the same spike again.
 *
 * inputs: data, upper_threshold, lower_threshold, window_size, sampling_rate
 *       
 *        data: this is the spike amplitudes for the time period over which spikes have to be found
 *        upper_threshold: generally this is (mean + (multiplying_factor) * (standard deviation))
 *        lower_threshold: generally this is (mean - (multiplying_factor) * (standard deviation))
 *        window_size: this is the size of the time window in sec - generally this is 0.001sec or 1msec
 *                     the window starts 0.15ms after the first threshold crossing and extends upto the time specified by window_size      
 *        start_time: this is the time of the first data point
 *                  
 * outputs:outputs
 *          
 *         outputs: this is an array of the indices of the voltage data that correspond to spikes
 *                  using the algorithm described above. The indices are converted to spike times by 
 *                  the original analysis program that calls this one.
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
  
      double *data, *output_data;
      int length_data, window_size, index, i, energy_data_index;
      double energy_data[1000000];
 
      if (nINP < 1) 
        mexErrMsgTxt("One inputs required");
      
      if (nOUT > 2) 
        mexErrMsgTxt("Too many output arguments");
      
      window_size = mxGetScalar(inputs[1]);
      
      printf("Window size is %i\n",window_size);
      
      data = (double *)mxGetPr(inputs[0]);
      
      length_data = mxGetM(inputs[0]) * mxGetN(inputs[0]);
      
      printf("Length of data is %i and the first and last data points are %g,%g\n",length_data,data[0],data[length_data - 1]);
      
      index = window_size;
      energy_data_index = 0;
      do
      {
          energy_data[energy_data_index] = 0;
          for(i = index; i < (index + window_size); i++)
          {
              energy_data[energy_data_index] = energy_data[energy_data_index] + (data[index] * data[index]);
          }
          energy_data[energy_data_index] = sqrt(energy_data[energy_data_index])/window_size;
          energy_data_index = energy_data_index + 1;
          index = index + 1;
      } while(index < (length_data - window_size)); 
      
      outputs[0] = mxCreateDoubleMatrix(energy_data_index,1,mxREAL);
      output_data = mxGetPr(outputs[0]);
      
      for (i = 0;i < energy_data_index;i++)
      {
        output_data[i] = energy_data[i];
      }  
     
}
  
  
		 
