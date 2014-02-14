function [spike_set] = getssmodelinfo(ssfile, ad_freq, ss_invert, ss_gain_external)
%GETSSMODELINFO get the model info for a SpikeSort inference.
%
%   Usage: [spike_set] = getssmodelinfo(ssfile, ad_freq, ss_invert, ss_gain_external)
%
%   Description:
%     Getssmodelinfo parses the model (.mdl) file associated with a
%     SpikeSort spike classification, returning a spike_set structure
%     which includes the number of models
%     and the waveforms and proportions for each model. Waveform
%     amplitudes should be in the units of Volts at the electrode.
%
%     Output spike_set structure fields: 
%       model                Cell array of n_model model structures.
%       n_models             Number of neural models.
%
%     Each model structure contains the fields:
%       x_data               Time axis (ms) for waveform model
%                              segment (len X 1).
%       y_data               Model Waveform (in Volts as seen at
%                              the electrode) (len X 1).
%       len                  Number of points in the model waveform.
%       len_ms               Duration of model waveform in msec.
%       pre_len              Number of points preceding the model
%                              waveform peak.
%       post_len             len = pre_len + post_len + 1.
%       pre_len_ms           pre_len expressed in msec.
%       class                Model number assigned to this model.
%       prop                 Approximate proportion of data in this
%                              model.
%   Input variables: 
%     ssfile:              SpikeSort .mdl file.
%     ad_freq:             Sampling rate of the waveform data.
%     ss_invert:           Waveform inversion setting (+/- 1) from SpikeSort.
%     ss_gain_external:    External gain used in SpikeSort.
%
%   Output variables:
%     spike_set:           As spike set structure with the fields
%                            described above.
%  
%   Uses: ss_include_globals
%
%   Author(s): Brian D. Wright <bdwright@phy.ucsf.edu>

ss_include_globals;

% Parse out comment lines from mdl files of the form:
% # Model 6, number of model points = 80, number of spikes = 43332
fid = fopen(ssfile,'rt');
if fid == -1
  error(['Can not open file ', ssfile])
end

imdl = -1;
while 1
  if imdl ~= -1
    commentline = fgetl(fid);
  end
  commentline = fgetl(fid);
  if commentline == -1
    break
  end
  stridx = findstr(commentline,'Model');
  commaidx = findstr(commentline,',');
  imdl = str2double(commentline(stridx+6:commaidx(1)-1));
  stridx = findstr(commentline,'=');
  npts = str2double(commentline(stridx(1)+2:commaidx(2)-1));
  nspikes(imdl+1) = str2double(commentline(stridx(2)+2:length(commentline)));
  [ssdata,nread]  = fscanf(fid,'%f',[2,npts]);
  ss_times = ssdata(1,:)';
  ss_amps(:,imdl+1) = ssdata(2,:)';
end

fclose(fid);
%ss_amps = ss_amps - ones(npts,1)*mean(ss_amps); 

n_models = imdl+1
maxmodelno = imdl

spike_set.n_models = n_models;

for i=1:n_models
  spike_set.model{i}.prop = nspikes(i)./sum(nspikes);
  spike_set.model{i}.y_data = ss_invert*ss_amps(:,i)*ss_gain_external/SS_TICKS_PER_VOLT;
  spike_set.model{i}.x_data = ss_times;
  spike_set.model{i}.pre_len_ms = -ss_times(1);
  spike_set.model{i}.pre_len = round((-ss_times(1))*ad_freq/1000.0);
  spike_set.model{i}.len = npts;
  spike_set.model{i}.len_ms = npts*1000.0/ad_freq;
  spike_set.model{i}.post_len = npts - spike_set.model{i}.pre_len - 1;  
  spike_set.model{i}.class = i-1;
end
