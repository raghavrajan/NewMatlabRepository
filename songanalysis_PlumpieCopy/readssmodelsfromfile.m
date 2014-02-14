function [spiketimes,spikeclass,nummodelspikes,nspikes] = readssmodelsfromfile(spikefile,maxmodelno)
% READSSMODELS: get SpikeSort times and class assignments.
%
% readssmodels reads in the information in a given SpikeSort
% classification (.spk) file and outputs the model ID and time
% associated with each spike. Outliers are not included.
%
% Usage: [spiketimes, spikeclass,nummodelspikes,nspikes] = readssmodels(spikefile,maxmodelno)
% Input variables: 
%     spikefile            SpikeSort .spk file.
%     maxmodelno           Largest model number to use, e.g. to
%                             exclude outlier models.
%
% Output variables: 
%     spiketimes           Spike times (in ms) array.
%     spikeclass           Model numbers giving each spike's
%                             classification (starting at 0).
%     nummodelspikes       Number of spikes belonging to each model.
%     nspikes              Total number of spikes.


if ~isempty(spikefile)
  fid = fopen(spikefile,'rt');
  [ssdata,nread]  = fscanf(fid,'%d %f',[2,inf]);
  spikeclass = ssdata(1,:)';
  spiketimes = ssdata(2,:)';
  nspikes = nread/2;
else 
  error('You must supply a spike data file.')
end
fclose(fid); 
clear ssdata;

for modelno=0:maxmodelno
  modelind = find(spikeclass==modelno);
%  modeltimes = sptimes(modelind);
  nummodelspikes(modelno+1) = length(modelind);
  disp(['Number of spikes for model ', num2str(modelno), ' is ', ...
	num2str(nummodelspikes(modelno+1)), '.'])
end  

