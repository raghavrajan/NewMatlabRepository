 /* find_spikes_negative_first
 * Finding spikes in multi-unit data
 * MEX file
 * 
 * Raghav 2006
 * 
 * The program finds spikes using the following algorithm.
 * The algorithm follows one used by DeAngelis and Newsome (J.Neurosci. Feb 15, 1999)
 *
 * If the voltage crosses a lower threshold and if it crosses a upper threshold within a 
 * specific time-window following the lower threshold crossing, it is considered a spike.
 *
 * Then it stores 8 points before the minimum and 23 points after the minimum, which makes a total of 32
 * points for the waveform.
 *
 * The minimum value between the two threshold crossings is considered as the time of the spike and the program
 * jumps 32 samples to avoid picking up the same spike again.
 *
 * inputs: data, upper_threshold, lower_threshold, window_size, start_time
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
  
      double *data, *output_spikes, *plot_output, *output_waveforms;
      int length_data, window_size, index, spike_index, i, j, k, waveform_index, max_index;
      double lower_threshold, upper_threshold, start_time, time, max_value;
      double spike_times[100000];
      double spike_waveforms[320000];
      
      FILE *fp = fopen("test","w");
      
      if (nINP < 5) 
        mexErrMsgTxt("Five inputs required");
      
      if (nOUT > 2) 
        mexErrMsgTxt("Too many input arguments");
      
      lower_threshold = mxGetScalar(inputs[1]);
      upper_threshold = mxGetScalar(inputs[2]);
      window_size = mxGetScalar(inputs[3]);
      start_time = mxGetScalar(inputs[4]);
      
      printf("Lower threshold is %g, Upper threshold is %g and Window size is %i\n",lower_threshold,upper_threshold,window_size);
      
      data = (double *)mxGetPr(inputs[0]);
      
      length_data = mxGetM(inputs[0]) * mxGetN(inputs[0]);
      
      printf("Length of data is %i and the first and last data points are %g,%g\n",length_data,data[0],data[length_data - 1]);
      
      index = 32;
      spike_index = 0;
      waveform_index = 0;

      do
      {
        if (data[index] < lower_threshold)
        {
            for(i = index; i < (index + window_size); i++)
            {
                if (data[i] > upper_threshold)
                {
                    min_index = index;
                    min_value = data[i];
                    
                    for(k = index; k < i; k++)
                    {
                        if (data[k] < min_value)
                        {
                            min_index = k;
                            min_value = data[k];
                        }
                    }
                    index = min_index;
                    spike_times[spike_index] = index;
                    for(j = (index-8); j < (index + 23); j++)
                    {
                        time = start_time + (j)/32000.0;
                        fprintf(fp,"%g\t%g\n",time,data[j]);
                        spike_waveforms[waveform_index] = data[j];
                        waveform_index = waveform_index + 1;
                    }
                    index = index + 32;
                    spike_index = spike_index + 1;
                    break;
                }
            }
        }
      index = index + 1;
      } while(index < (length_data - 32)); 
      
      outputs[0] = mxCreateDoubleMatrix(spike_index,1,mxREAL);
      output_spikes = mxGetPr(outputs[0]);
      
      for (i = 0;i < spike_index;i++)
      {
        output_spikes[i] = spike_times[i];
      }  

      fclose(fp);
 
      
}
  
  
		 
