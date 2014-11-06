 /* tetrode_find_spikes
 * Finding spikes in multi-unit data from a tetrode
 * MEX file
 * 
 * Raghav 2006
 * 
 * The program finds spikes using the following algorithm.
 * The algorithm follows one used by DeAngelis and Newsome (J.Neurosci. Feb 15, 1999)
 *
 * If the voltage crosses an upper threshold on any one of the channels of the tetrode,
 * 32 points are stored on each channel - 8 before threshold crossing and 24 after threshold crossing
 *
 * The point of the first threshold crossing is considered as the time of the spike and the program
 * jumps one window_size to avoid picking up the same spike again.
 *
 * inputs: data, upper_threshold, lower_threshold, window_size, start_time
 *       
 *        data: this is the spike amplitudes for the time period over which spikes have to be found - 4 columns each from one tetrode channel
 *        upper_threshold: generally this is (mean + (multiplying_factor) * (standard deviation))
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
  
      double *data_ch1, *data_ch2, *data_ch3, *data_ch4, *output_spike_times, *out_waveforms_ch1, *out_waveforms_ch2, *out_waveforms_ch3, *out_waveforms_ch4, *upper_threshold;
      int length_data, window_size, index, spike_index, i, j, k, waveform_index, max_index,length_thresholds;
      double window_to_skip, start_time, time, max_value;
      double spike_times[100000];
      double spike_waveforms[320000];
      
/*      FILE *fp = fopen("TempSpikeTimes.txt","w");
      FILE *fp1 = fopen("TempSpikeWaveformsCh1.txt","w");
      FILE *fp2 = fopen("TempSpikeWaveformsCh2.txt","w");
      FILE *fp3 = fopen("TempSpikeWaveformsCh3.txt","w");
      FILE *fp4 = fopen("TempSpikeWaveformsCh4.txt","w"); */
      
      if (nINP < 7) 
        mexErrMsgTxt("Eight inputs required");
      
      if (nOUT > 5) 
        mexErrMsgTxt("Too many output arguments");
      
      upper_threshold = mxGetPr(inputs[4]);
      length_thresholds = mxGetM(inputs[4]) * mxGetN(inputs[4]);
      
      window_size = mxGetScalar(inputs[5]);
      start_time = mxGetScalar(inputs[6]);
      window_to_skip = mxGetScalar(inputs[7]);
      
      printf("Start time is %g, window_to_skip is %g and window_size is %i\n",start_time,window_to_skip,window_size);
      printf("The thresholds for the 4 channels are %g, %g, %g and %g\n",upper_threshold[0],upper_threshold[1],upper_threshold[2],upper_threshold[3]);
      
      data_ch1 = (double *)mxGetPr(inputs[0]);
      data_ch2 = (double *)mxGetPr(inputs[1]);
      data_ch3 = (double *)mxGetPr(inputs[2]);
      data_ch4 = (double *)mxGetPr(inputs[3]);            
      
      length_data = mxGetM(inputs[0]) * mxGetN(inputs[0]);
      
      printf("Length of data is %i \n",length_data);
      printf("First points for all 4 tetrode channels are %g, %g, %g and %g\n",data_ch1[0], data_ch2[0],data_ch3[0],data_ch4[0]);
      
      index = 32;
      spike_index = 0;
      waveform_index = 0;

      do
      {
        if ((data_ch1[index] > upper_threshold[0]) || (data_ch2[index] > upper_threshold[1]) || (data_ch3[index] > upper_threshold[2]) || (data_ch4[index] > upper_threshold[3]))
        {
            spike_index = spike_index + 1;
            
            if ((data_ch1[index] > upper_threshold[0]))
            {
                max_index = index;
                max_value = data_ch1[index];
                for(j = index; j < (index + 2*window_size); j++)
                {
                    if (data_ch1[j] > max_value)
                    {
                        max_index = j;
                        max_value = data_ch1[j];
                    }
                }
                index = max_index;
            }
            
            if ((data_ch2[index] > upper_threshold[1]))
            {
                max_index = index;
                max_value = data_ch2[index];
                for(j = index; j < (index + 2*window_size); j++)
                {
                    if (data_ch2[j] > max_value)
                    {
                        max_index = j;
                        max_value = data_ch2[j];
                    }
                }
                index = max_index;
            }
            
            if ((data_ch3[index] > upper_threshold[2]))
            {
                max_index = index;
                max_value = data_ch3[index];
                for(j = index; j < (index + 2*window_size); j++)
                {
                    if (data_ch3[j] > max_value)
                    {
                        max_index = j;
                        max_value = data_ch3[j];
                    }
                }
                index = max_index;
            }
            
            if ((data_ch4[index] > upper_threshold[2]))
            {
                max_index = index;
                max_value = data_ch4[index];
                for(j = index; j < (index + 2*window_size); j++)
                {
                    if (data_ch4[j] > max_value)
                    {
                        max_index = j;
                        max_value = data_ch4[j];
                    }
                }
                index = max_index;
            }
            index = index + window_to_skip - 1;
        }
        index = index + 1;
      } while(index < (length_data - 32));

      printf("Found %i spikes\n",spike_index);  

      outputs[0] = mxCreateDoubleMatrix(spike_index,1,mxREAL);
      output_spike_times = mxGetPr(outputs[0]);
      
      outputs[1] = mxCreateDoubleMatrix(32, spike_index,mxREAL);
      out_waveforms_ch1 = mxGetPr(outputs[1]);
      
      outputs[2] = mxCreateDoubleMatrix(32, spike_index,mxREAL);
      out_waveforms_ch2 = mxGetPr(outputs[2]);
            
      outputs[3] = mxCreateDoubleMatrix(32, spike_index,mxREAL);
      out_waveforms_ch3 = mxGetPr(outputs[3]);
      
      outputs[4] = mxCreateDoubleMatrix(32, spike_index,mxREAL);
      out_waveforms_ch4 = mxGetPr(outputs[4]);
      
      spike_index = 0;
      waveform_index = 0;
      
      index = 32;
      do
      {
        if ((data_ch1[index] > upper_threshold[0]) || (data_ch2[index] > upper_threshold[1]) || (data_ch3[index] > upper_threshold[2]) || (data_ch4[index] > upper_threshold[3]))
        {
            
            if ((data_ch1[index] > upper_threshold[0]))
            {
                max_index = index;
                max_value = data_ch1[index];
                for(j = index; j < (index + 2*window_size); j++)
                {
                    if (data_ch1[j] > max_value)
                    {
                        max_index = j;
                        max_value = data_ch1[j];
                    }
                }
                index = max_index;
            }
            
            if ((data_ch2[index] > upper_threshold[1]))
            {
                max_index = index;
                max_value = data_ch2[index];
                for(j = index; j < (index + 2*window_size); j++)
                {
                    if (data_ch2[j] > max_value)
                    {
                        max_index = j;
                        max_value = data_ch2[j];
                    }
                }
                index = max_index;
            }
            
            if ((data_ch3[index] > upper_threshold[2]))
            {
                max_index = index;
                max_value = data_ch3[index];
                for(j = index; j < (index + 2*window_size); j++)
                {
                    if (data_ch3[j] > max_value)
                    {
                        max_index = j;
                        max_value = data_ch3[j];
                    }
                }
                index = max_index;
            }
            
            if ((data_ch4[index] > upper_threshold[2]))
            {
                max_index = index;
                max_value = data_ch4[index];
                for(j = index; j < (index + 2*window_size); j++)
                {
                    if (data_ch4[j] > max_value)
                    {
                        max_index = j;
                        max_value = data_ch4[j];
                    }
                }
                index = max_index;
            }
                    
            time = start_time + (index)/32000.0;
            output_spike_times[spike_index] = time;
            spike_index = spike_index + 1;
/*            fprintf(fp,"%g\n",time); */
            
            for(i = (index - window_size/4 - 1); i < (index + (3 * window_size/4) - 1); i++)
            {
                out_waveforms_ch1[waveform_index] = data_ch1[i];
                out_waveforms_ch2[waveform_index] = data_ch2[i];
                out_waveforms_ch3[waveform_index] = data_ch3[i];
                out_waveforms_ch4[waveform_index] = data_ch4[i];
                waveform_index = waveform_index + 1;
/*                fprintf(fp1,"%g\t",data_ch1[i]);
                fprintf(fp2,"%g\t",data_ch2[i]);
                fprintf(fp3,"%g\t",data_ch3[i]);
                fprintf(fp4,"%g\t",data_ch4[i]); */
            }
/*            fprintf(fp1,"\n");
            fprintf(fp2,"\n");
            fprintf(fp3,"\n");
            fprintf(fp4,"\n"); */
            
            index = index + window_to_skip - 1;
        }
        index = index + 1;
      } while(index < (length_data - 32));

/*      fclose(fp);
      fclose(fp1);
      fclose(fp2);
      fclose(fp3);
      fclose(fp4); */
 
      
}
  
  
		 
