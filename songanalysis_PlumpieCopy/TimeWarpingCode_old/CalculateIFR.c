 /* CalculateIFR
 * Given a spike train, calculate the Instantaneous Firing Rate
 * MEX file
 * 
 * Raghav 2008
 * 
 * This program finds the instantaneous firing rate defined at each time point as follows:
 * IFR(t) = 1/(SpikeTime(t1) - SpikeTime(t2)) where t1 > t > t2
 *
 * inputs: spiketimes
 *       
 *        spiketimes: this is the spike times for each spike train - each row is one trial or one rendition of song
 *                  
 * outputs:outputs
 *          
 *         outputs: this is the instantaneous firing rates for each spike train.
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
  
      double *spiketimes, *output_ifr;
      int length_data, window_size, index, spike_index, i, j, k, waveform_index, max_index, skip_window_size;
      double IFR[100000];
      
      if (nINP < 1) 
        mexErrMsgTxt("One input required");
      
      if (nOUT > 1) 
        mexErrMsgTxt("Too many output arguments");
      
      spiketimes = (double *)mxGetPr(inputs[0]);
      
      length_data = mxGetM(inputs[0]) * mxGetN(inputs[0]);
      
      /* printf("Length of data is %i and the first and last data points are %g,%g\n",length_data,data[0],data[length_data - 1]); */
      
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
                    
                    for(k = index; k < i; k++)
                    {
                        if (data[k] > max_value)
                        {
                            max_index = k;
                            max_value = data[k];
                        }
                    }
                    index = max_index;

                    index = index + skip_window_size;
                    spike_index = spike_index + 1;
                    break;
                }
            }
        }
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
                    
                    index = index + skip_window_size;
                    spike_index = spike_index + 1;
                    break;
                }
            }
        }
      index = index + 1;
      } while(index < (length_data - skip_window_size)); 
      
}
  
