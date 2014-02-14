% Generalize later
ad_freq = 32000;

% This enables a generic batch processing mode.
% Interactive queries will only happen if necessary variables
% are not present.
if ~exist('do_batch','var')
  do_batch = 0
end

%%
%% Get a list of files from the user and extract file properties.
%%

if ~do_batch | ~exist('use_batch_file','var')
  use_batch_file = input(['Do you want to (0) process a single file or (1) process a batch file? [1]: ']);
end
if isempty(use_batch_file)
  use_batch_file = 1
end

% What data directory and file(s) to use?

if use_batch_file
  % Get batch file information and open it.
  if ~do_batch | ~exist('batch_info','var')
    batch_info = [];
  end
  batch_info = batch_select(batch_info);
  % Read in the file list
  batch_info = batch_read(batch_info)

  for ifile=1:batch_info.nfiles
    datapath = fullfile(batch_info.basedir,batch_info.filenames{ifile});
    % [datafiles{ifile}.dir,dname,dext,dver] = fileparts(datapath)
    % datafiles{ifile}.name = [dname, dext, dver];
    % Use full path whenever possible!
    datafiles{ifile}.name = datapath;
    datafiles{ifile}.type = 'auto';
  end
else
  if ~do_batch | ~exist('datafile','var')
    datafile = input(['Enter name of data file (typically rec or dcp file): '], 's');
  end
  while 1
    if ~exist(datafile,'file')
      datafile = input('Reenter name of data file: ', 's');
    else
      break;
    end
  end
  % [datafiles{1}.dir,dname,dext,dver] = fileparts(datafile);
  % datafiles{1}.name = [dname, dext, dver];
  datafiles{1}.name = datafile;
  datafiles{1}.type = 'auto';
end

% Get spike classification files: model id. and times
% Guess the correct directory from the first data file path.
[datadir,dname,dext,dver] = fileparts(datafiles{1}.name);
if ~do_batch | ~exist('spkdir','var')
  spkdir = input(['Enter name of directory containing SpikeSort output files [', datadir, ']: '], 's');
  if isempty(spkdir)
    spkdir = datadir
  end
end
while 1
  if ~exist(spkdir,'dir')
    spkdir = input('Reenter name of directory containing SpikeSort output files: ', 's');
  else
    break
  end
end

% Get syllable segmentation files:
% Guess the correct directory from the first data file path.
[datadir,dname,dext,dver] = fileparts(datafiles{1}.name);
if ~do_batch | ~exist('notedir','var')
  notedir = input(['Enter name of directory containing uisonganal output files [', datadir, ']: '], 's');
  if isempty(notedir)
    notedir = datadir
  end
end
while 1
  if ~exist(notedir,'dir')
    notedir = input('Reenter name of directory containing uisonganal output files: ', 's');
  else
    break
  end
end

% Load in .spk and .not.mat files
if use_batch_file
  for ifile=1:batch_info.nfiles
    [fdir,fname,fext,fver] = fileparts(batch_info.filenames{ifile});
    spkfiles{ifile}.name = fullfile(spkdir, [fname, '.spk']);
    mdlfiles{ifile}.name = fullfile(spkdir, [fname, '.mdl']);
    spkfiles{ifile}.type = 'ss';
    notefiles{ifile}.name = fullfile(notedir, [fname, fext, '.not.mat']);
  end
else
    [fdir,fname,fext,fver] = fileparts(datafiles{1}.name);
    spkfiles{1}.name = fullfile(spkdir, [fname, '.spk']);
    mdlfiles{1}.name = fullfile(spkdir, [fname, '.mdl']);
    spkfiles{1}.type = 'ss';
    notefiles{1}.name = fullfile(notedir, [fname, fext, '.not.mat']);
end

%%
%%  Get some extra user settings.
%%
% Do we invert the model
if ~do_batch | ~exist('ss_invert','var')
  yesno = input('Did you invert the waveforms in SpikeSort (Y/n): ', 's');  
end
if isempty(yesno) | strncmpi(yesno,'y',1)
  ss_invert = -1
else
  ss_invert = 1
end

% External gain?
if ~do_batch | ~exist('ss_gain_external','var')
  ss_gain_external = input(['Enter the external gain on the neural channel [1000]: ']);  
end
if isempty(ss_gain_external)
  ss_gain_external = 1000.0
end

% Get the relevant information for pulling out the spikes for a 
% particular syllable/sequence of syllables

syl_seq=input('What syllable/sequence of syllables? ','s');
last_syl=syl_seq(length(syl_seq));
syl_seq_len = length(syl_seq)

selmdl = input('Enter the spike model to analyze: ');
nselmdl = length(selmdl);

latency_ms=input('What is the latency between the song and the neural activity?  (msec)   ');
% range in which the neural data will be scanned:  
% syllable/sequence onset time-latency_ms to 
% (syllable/sequence offset_time-syllable/sequence onset_time) - latency_ms

% Get spike times and syllable onsets and offsets for each motif

iseq = 0;

for ifile=1:batch_info.nfiles

  % Set up the spike_set.model structures. We assume all files are
  % using the same models!
  [spike_set] = getssmodelinfo_new(mdlfiles{ifile}.name, ad_freq, ss_invert, ss_gain_external);

  [spiketimes,spikeclass,nummodelspikes,nspikes]= ...
      readssmodelsfromfile(spkfiles{ifile}.name,spike_set.n_models-1);

  % Will define "trials" based on motifs 
  % Reset this for each file.
  mdlidx = [];
  for i=1:nselmdl
    mdlidx = [mdlidx; find(spikeclass==selmdl(i))];
  end
  mdlidx = sort(mdlidx)
  mdltimes = spiketimes(mdlidx);
  
  % Get the notefile and onset times of the syllable/sequence
  if exist([notefiles{ifile}.name])
    load(notefiles{ifile}.name); %Fs, onsets, offsets and labels are defined
  else 
    disp(['Cannot find',notefiles{ifile}.name])
    break
  end
  
  % count the number of matches in this file
  matches=findstr(syl_seq,labels);    %# of matches for the entire sequence
  matches_last=findstr(last_syl,labels);  %# of matches for the last syllable
  
  %if length(matches)~=length(matches_last)
  %  disp('Some syllables occur outside of the sequence.')
  %  break
  %end
  num_seq_matches_file=length(matches)

  for i=1:num_seq_matches_file
    % iseq is a global "motif"/sequence counter
    iseq = iseq + 1
    for isyl = 1:syl_seq_len
      seq(iseq).syl(isyl).start = onsets(matches(i)+isyl-1);
      seq(iseq).syl(isyl).end = offsets(matches(i)+isyl-1);
      seq(iseq).syl(isyl).dur = seq(iseq).syl(isyl).end - ...
          seq(iseq).syl(isyl).start;
      seq(iseq).syl(isyl).label = syl_seq(isyl);
    end
    for isyl = 1:syl_seq_len-1
      seq(iseq).interval(isyl).start = offsets(matches(i)+isyl-1);
      seq(iseq).interval(isyl).end = onsets(matches(i)+isyl);
      seq(iseq).interval(isyl).dur = seq(iseq).interval(isyl).end - ...
          seq(iseq).interval(isyl).start;
    end
    
    % "Motif"/sequence level info.
    seq(iseq).start = onsets(matches(i));
    seq(iseq).end = offsets(matches(i)+syl_seq_len-1);
    seq(iseq).dur = seq(iseq).end - seq(iseq).start;
    seq(iseq).nsyl = syl_seq_len;
    seq(iseq).labels = syl_seq;
    seq(iseq).file = notefiles{ifile}.name;
    
    for imdl=0:spike_set.n_models-1
      seqspktimes = ...
          mdltimes(find(mdltimes >= seq(iseq).start - latency_ms ...
                        & mdltimes < seq(iseq).end - latency_ms ));
      seq(iseq).spike_times{1} = seqspktimes - seq(iseq).start ...
          + latency_ms;
      for isyl=1:syl_seq_len
        seq(iseq).syl(isyl).spike_times{1} = ...
            seqspktimes(find(seqspktimes >= seq(iseq).syl(isyl).start - latency_ms ...
                            & seqspktimes < seq(iseq).syl(isyl).end ...
                            - latency_ms)) - seq(iseq).start + latency_ms;    
      end
      for isyl=1:syl_seq_len-1
        seq(iseq).interval(isyl).spike_times{1} = ...
            seqspktimes(find(seqspktimes >= seq(iseq).interval(isyl).start - latency_ms ...
                            & seqspktimes < seq(iseq).interval(isyl).end ...
                            - latency_ms)) - seq(iseq).start + latency_ms;    
      end
    end
  
  end % loop over motif matches

end % loop over files

nseqs = iseq;

% This sequence will be the reference
iref = 1;

% OK now do a piecewise linear warp.

for iseq=1:nseqs
  if iseq ~= iref
    seq(iseq).spike_times_warp{1}(:) = [];
    
    for k=1:2*syl_seq_len-1
      isyl = ceil(k/2)
      if rem(k,2)
        is_syllable = 1;
      else
        is_syllable = 0;
      end
      if is_syllable == 1
        r = seq(iref).syl(isyl).dur/seq(iseq).syl(isyl).dur
        offset = seq(iseq).syl(isyl).start - seq(iseq).start;
        origin = seq(iref).syl(isyl).start - seq(iref).start;
      else
        r = seq(iref).interval(isyl).dur/seq(iseq).interval(isyl).dur
        offset = seq(iseq).interval(isyl).start - seq(iseq).start;
        origin = seq(iref).interval(isyl).start - seq(iref).start;
      end        

      if is_syllable
        seq(iseq).syl(isyl).spike_times_warp{1}(:) = [];
        seq(iseq).syl(isyl).spike_times_warp{1} = ...
            r*(seq(iseq).syl(isyl).spike_times{1}(:) - ... 
               offset) + origin;
        % Concatenate everything together for final warp result.
        seq(iseq).spike_times_warp{1} = [ ...
            seq(iseq).spike_times_warp{1}(:); seq(iseq).syl(isyl).spike_times_warp{1}(:)];     
      else
        seq(iseq).interval(isyl).spike_times_warp{1}(:) = [];
        seq(iseq).interval(isyl).spike_times_warp{1} = ...
            r*(seq(iseq).interval(isyl).spike_times{1}(:) - ... 
               offset) + origin;
        % Concatenate everything together for final warp result.
        seq(iseq).spike_times_warp{1} = [ ...
            seq(iseq).spike_times_warp{1}(:); seq(iseq).interval(isyl).spike_times_warp{1}(:)];     
      end
    end
    
  else
    for isyl = 1:syl_seq_len
      seq(iseq).syl(isyl).spike_times_warp{1} = ...
          seq(iseq).syl(isyl).spike_times{1};
    end
    for isyl = 1:syl_seq_len-1
      seq(iseq).interval(isyl).spike_times_warp{1} = ...
          seq(iseq).interval(isyl).spike_times{1};
    end
    seq(iseq).spike_times_warp{1} = seq(iseq).spike_times{1};
    
  end

end % loop over "motifs"


% Now plot unwarped and warped rasters.

figure;
title(['Batch file: ', batch_info.batchfilename, ' Model(s): ', num2str(selmdl)])
subplot(2,1,1)
for iseq = 1:nseqs
  n_events = length(seq(iseq).spike_times{1}(:));
  t_rast = ones(2,1)*seq(iseq).spike_times{1}(:)';
  y_rast = [iseq-0.3; iseq+0.3]*ones(1,n_events);
  line(t_rast,y_rast,'Color','k','LineWidth',1);
  if iseq == 1
    axislims = [0.0 max([seq.dur] + latency_ms) 0 nseqs+1] 
    axis(axislims);
    axis ij;
    %xlabel('t (msec)')
    ylabel('Trial number')
    title(['Original data'])
  end
end

subplot(2,1,2)
for iseq = 1:nseqs
  n_events = length(seq(iseq).spike_times_warp{1}(:));
  t_rast = ones(2,1)*seq(iseq).spike_times_warp{1}(:)';
  y_rast = [iseq-0.3; iseq+0.3]*ones(1,n_events);
  line(t_rast,y_rast,'Color','k','LineWidth',1);
  if iseq == 1
    axislims = [0.0 max([seq.dur]) + latency_ms 0 nseqs+1] 
    axis(axislims);
    axis ij;
    xlabel('t (msec)')
    ylabel('Trial number')
    title(['Warped data'])
  end
end
%hold on
%for isyl=1:seq(iref).nsyl
%  t_off = seq(iref).syl(isyl).start - seq(iref).syl(1).start + latency_ms
%  plot([t_off t_off], [0 nseqs+1], 'g-')
%  t_off = seq(iref).syl(isyl).end - seq(iref).syl(1).start + latency_ms
%  plot([t_off t_off], [0 nseqs+1], 'r-')
%end
%hold off

if 0
figure;
title(['Syllables, Batch file: ', batch_info.batchfilename])
for iseq = 1:nseqs
  y_off = iseq+0.5;
  for isyl = 1:seq(iseq).nsyl
    line([seq(iseq).syl(isyl).start - seq(iseq).syl(1).start ...
          seq(iseq).syl(isyl).end - seq(iseq).syl(1).start], [y_off ...
                        y_off], 'Color', 'r','LineWidth',1 );
  end
  if iseq == 1
    axislims = [0.0 max([seq.dur]) 0 nseqs+1] 
    axis(axislims);
    xlabel('t (sec)')
    ylabel('Trial number')
  end
end
end

nseqs_plot = nseqs
figure;
title(['Batch file: ', batch_info.batchfilename, ' Model(s): ', num2str(selmdl)])
subplot(2,1,1)
for iseq = 1:nseqs_plot
  %seq_off = nseqs_plot-iseq+1;
  seq_off = iseq;
  n_events = length(seq(iseq).spike_times{1}(:));
  t_rast = ones(2,1)*seq(iseq).spike_times{1}(:)';
  y_rast = [seq_off-0.3; seq_off+0.3]*ones(1,n_events);
  line(t_rast,y_rast,'Color','k','LineWidth',1);
  if iseq == 1
    axislims = [0.0 max([seq.dur])+latency_ms 0 nseqs_plot+1] 
    axis(axislims);
    axis ij;
    %xlabel('t (sec)')
    ylabel('Trial number')
    title(['Original data'])
  end
  hold on
  y_off = seq_off+0.4;
  for isyl = 1:seq(iseq).nsyl
    line([seq(iseq).syl(isyl).start - seq(iseq).syl(1).start + latency_ms ...
          seq(iseq).syl(isyl).end - seq(iseq).syl(1).start + latency_ms], [y_off ...
                        y_off], 'Color', 'r','LineWidth',0.5 );
  end
  hold off
end

subplot(2,1,2)
for iseq = 1:nseqs_plot
  %seq_off = nseqs_plot-iseq+1;
  seq_off = iseq;
  n_events = length(seq(iseq).spike_times_warp{1}(:));
  t_rast = ones(2,1)*seq(iseq).spike_times_warp{1}(:)';
  y_rast = [seq_off-0.3; seq_off+0.3]*ones(1,n_events);
  line(t_rast,y_rast,'Color','k','LineWidth',1);
  if iseq == 1
    axislims = [0.0 max([seq.dur])+latency_ms 0 nseqs_plot+1] 
    axis(axislims);
    axis ij
    %xlabel('t (sec)')
    ylabel('Trial number')
    title(['Warped data'])
  end
end
hold on
for isyl=1:seq(iref).nsyl
  t_off = seq(iref).syl(isyl).start - seq(iref).syl(1).start + latency_ms
  plot([t_off t_off], [0 nseqs_plot+1], 'g-')
  t_off = seq(iref).syl(isyl).end - seq(iref).syl(1).start + latency_ms
  plot([t_off t_off], [0 nseqs_plot+1], 'r-')
end
xlabel('t (msec)')
hold off

% Make trial raster histograms for mutual information calculations

% Make a 1 msec resolution raster.
dt_hires = 1.0
nrasterbins = fix(seq(iref).dur/dt_hires); % All sequences warped to iref duration
raster_binedges = 0.0:dt_hires:dt_hires*nrasterbins;
imdl = 1;
for iseq=1:nseqs
  
  % Raster is the spike train histogram, binned to resolution
  % dt, for each trial.

  % Note histc adds one more bin for when the value is equal to
  % the last bin edge.
  [raster{imdl}(:,iseq),raster_bins] = histc(seq(iseq).spike_times_warp{imdl}(:),raster_binedges);
end


[bdir,bname,bext,bver] = fileparts(batch_info.batchfilename);
%savefile = [bdir, bname, '_pwlwarp_save.mat']
raster_savefile = [bdir, bname,  '_rasters.mat']

save(raster_savefile, 'seq', 'nseqs', 'iref', 'dt_hires', 'raster', ...
    'raster_binedges', 'selmdl', 'latency_ms', ...
    'batch_info', 'ss_invert', 'ss_gain_external', 'syl_seq');

%save(raster_savefile, 'seq', 'nseqs', 'selmdl', 'syl_seq', 'raster', ...
%    'raster_binedges', 'dt_hires');


