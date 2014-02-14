 /* FindSpikes
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
      int length_data, window_size, index, spike_index, i, j, k, waveform_index, max_index, skip_window_size,increment_index;
      double lower_threshold, upper_threshold, start_time, time, max_value, sampling_rate;
      double spike_times[100000];
      double spike_waveforms[3200000];
      
      if (nINP < 7) 
        mexErrMsgTxt("Seven inputs required");
      
      if (nOUT > 2) 
        mexErrMsgTxt("Too many output arguments");
      
      lower_threshold = mxGetScalar(inputs[1]);
      upper_threshold = mxGetScalar(inputs[2]);
      window_size = mxGetScalar(inputs[3]);
      skip_window_size = mxGetScalar(inputs[4]);
      start_time = mxGetScalar(inputs[5]);
      sampling_rate = mxGetScalar(inputs[6]);
      
      /* printf("Lower threshold is %g, Upper threshold is %g and Window size is %i\n",lower_threshold,upper_threshold,window_size); */
      
      data = (double *)mxGetPr(inputs[0]);
      
      length_data = mxGetM(inputs[0]) * mxGetN(inputs[0]);
      
      printf("Length of data is %i and the first and last data points are %g,%g\n",length_data,data[0],data[length_data - 1]);
      
      index = skip_window_size;
      spike_index = 0;
      waveform_index = 0;

      do
      {
        if (data[index] > upper_threshold)
        {
            for(i = index; i < (index + window_size); i++)
            {
                if (data[i] < lower_threshold)
                {
                    max_index = index;
                    max_value = data[i];
                    increment_index = i;
                    
                    for(k = index; k < i; k++)
                    {
                        if (data[k] > max_value)
                        {
                            max_index = k;
                            max_value = data[k];
                        }
                    }
                    /*index = max_index;

                    index = index + skip_window_size;*/
                    index = increment_index;
                    
                    spike_index = spike_index + 1;
                    break;
                }
                
            }
        }
/*        printf("Index is %i\n",index);*/
      index = index + 1;
      } while(index < (length_data - skip_window_size)); 
      
      printf("Found %i spikes\n",spike_index);  

      outputs[0] = mxCreateDoubleMatrix(spike_index,1,mxREAL);
      output_spikes = mxGetPr(outputs[0]);
      
      outputs[1] = mxCreateDoubleMatrix(32, spike_index,mxREAL);
      output_waveforms = mxGetPr(outputs[1]);
      
      index = skip_window_size;
      spike_index = 0;
      waveform_index = 0;
      
      do
      {
        if (data[index] > upper_threshold)
        {
            for(i = index; i < (index + window_size); i++)
            {
                if (data[i] < lower_threshold)
                {
                    max_index = index;
                    max_value = data[i];
                    increment_index = i;
                    
                    for(k = index; k < i; k++)
                    {
                        if (data[k] > max_value)
                        {
                            max_index = k;
                            max_value = data[k];
                        }
                    }
                    index = max_index;
                    time = start_time + index/sampling_rate;
                    output_spikes[spike_index] = time;
                    
                    for(i = (index - window_size/4 - 1); i < (index + (3 * window_size/4) - 1); i++)
                    {
                        output_waveforms[waveform_index] = data[i];
                        waveform_index = waveform_index + 1;
                    }
                    
                    /*index = index + skip_window_size;*/
                    index = increment_index;
                    spike_index = spike_index + 1;
                    break;
                }
            }
        }
      index = index + 1;
      } while(index < (length_data - skip_window_size)); 
      
}
  
