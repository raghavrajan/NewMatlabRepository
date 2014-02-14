function [] = make_time_warped_raster_pst(SpikeFile,NoteFile,bin_size,motif,RawDataFile)

SpikeTimes = load(SpikeFile);
% RawSpikeTimes = load(SpikeFile);
% SpikeTimes = RawSpikeTimes(find(RawSpikeTimes(:,1) < 3),2);
% SpikeTimes = SpikeTimes/1000;
Notes = load(NoteFile);

Notes.onsets = Notes.onsets/1000;
Notes.offsets = Notes.offsets/1000;

for i = 1:length(motif),
    syllable(:,i) = find(Notes.labels == motif(i));
    syllable_lengths(:,i) = Notes.offsets(syllable(:,i)) - Notes.onsets(syllable(:,i));
end

for i = 1:(length(motif)-1),
    gap_lengths(:,i) = Notes.onsets(syllable(:,i+1)) - Notes.offsets(syllable(:,i));
end

motif_durations = Notes.offsets(syllable(:,end)) - Notes.onsets(syllable(:,1));
[sorted_motif_durations,sorted_indices] = sort(motif_durations);

if (mod(length(motif_durations),2) == 0)
    median_motif = sorted_indices(length(motif_durations)/2);
else
    median_motif = sorted_indices((length(motif_durations) + 1)/2);
end

median_syllable_lengths = syllable_lengths(median_motif,:);
median_gap_lengths = gap_lengths(median_motif,:);

figure;
Raster = [];
trial_pst = [];
for i = 1:size(motif_durations,1),
    TrialSpikeTimes = SpikeTimes(find((SpikeTimes > (Notes.onsets(syllable(i,1)) - 0.2)) & (SpikeTimes < Notes.onsets(syllable(i,1)))));
    TrialSpikeTimes = TrialSpikeTimes - Notes.onsets(syllable(i,1));
    y_value = ones(length(TrialSpikeTimes),1) * i/5;
    Raster = [Raster; [TrialSpikeTimes y_value]];

    syll_count = 0;
    gap_count = 0;
    for j = 1:(2*length(motif) - 1),
        if (mod(j,2) == 1)
            syll_count = syll_count + 1;
            TrialSpikeTimes = SpikeTimes(find((SpikeTimes > Notes.onsets(syllable(i,syll_count))) & (SpikeTimes < Notes.offsets(syllable(i,syll_count)))));
            TrialSpikeTimes = (TrialSpikeTimes - Notes.onsets(syllable(i,syll_count)))/syllable_lengths(i,syll_count);
            TrialSpikeTimes = (TrialSpikeTimes * median_syllable_lengths(1,syll_count)) + Notes.onsets(syllable(median_motif,syll_count)) - Notes.onsets(syllable(median_motif,1));
            y_value = ones(length(TrialSpikeTimes),1) * i/5;
            Raster = [Raster; [TrialSpikeTimes y_value]];
        else
            gap_count = gap_count + 1;
            TrialSpikeTimes = SpikeTimes(find((SpikeTimes > Notes.offsets(syllable(i,gap_count))) & (SpikeTimes < Notes.onsets(syllable(i,(gap_count + 1))))));
            TrialSpikeTimes = (TrialSpikeTimes - Notes.offsets(syllable(i,gap_count)))/gap_lengths(i,gap_count);
            TrialSpikeTimes = (TrialSpikeTimes * median_gap_lengths(1,gap_count)) + Notes.offsets(syllable(median_motif,gap_count)) - Notes.onsets(syllable(median_motif,1));
            y_value = ones(length(TrialSpikeTimes),1) * i/5;
            Raster = [Raster; [TrialSpikeTimes y_value]];
        end
        
%         switch j
%             case 1
%                 TrialSpikeTimes = SpikeTimes(find((SpikeTimes > Notes.onsets(syllable(i,1))) & (SpikeTimes < Notes.offsets(syllable(i,1)))));
%                 TrialSpikeTimes = (TrialSpikeTimes - Notes.onsets(syllable(i,1)))/syllable_lengths(i,1);
%                 TrialSpikeTimes = (TrialSpikeTimes * median_syllable_lengths(1,1)) + Notes.onsets(syllable(median_motif,1)) - Notes.onsets(syllable(median_motif,1));
%                 y_value = ones(length(TrialSpikeTimes),1) * i/5;
%                 Raster = [Raster; [TrialSpikeTimes y_value]];
%             case 2
%                 TrialSpikeTimes = SpikeTimes(find((SpikeTimes > Notes.offsets(syllable(i,1))) & (SpikeTimes < Notes.onsets(syllable(i,2)))));
%                 TrialSpikeTimes = (TrialSpikeTimes - Notes.offsets(syllable(i,1)))/gap_lengths(i,1);
%                 TrialSpikeTimes = (TrialSpikeTimes * median_gap_lengths(1,1)) + Notes.offsets(syllable(median_motif,1)) - Notes.onsets(syllable(median_motif,1));
%                 y_value = ones(length(TrialSpikeTimes),1) * i/5;
%                 Raster = [Raster; [TrialSpikeTimes y_value]];
%             case 3
%                 TrialSpikeTimes = SpikeTimes(find((SpikeTimes > Notes.onsets(syllable(i,2))) & (SpikeTimes < Notes.offsets(syllable(i,2)))));
%                 TrialSpikeTimes = (TrialSpikeTimes - Notes.onsets(syllable(i,2)))/syllable_lengths(i,2);
%                 TrialSpikeTimes = (TrialSpikeTimes * median_syllable_lengths(1,2)) + Notes.onsets(syllable(median_motif,2)) - Notes.onsets(syllable(median_motif,1));
%                 y_value = ones(length(TrialSpikeTimes),1) * i/5;
%                 Raster = [Raster; [TrialSpikeTimes y_value]];
%             case 4
%                 TrialSpikeTimes = SpikeTimes(find((SpikeTimes > Notes.offsets(syllable(i,2))) & (SpikeTimes < Notes.onsets(syllable(i,3)))));
%                 TrialSpikeTimes = (TrialSpikeTimes - Notes.offsets(syllable(i,2)))/gap_lengths(i,2);
%                 TrialSpikeTimes = (TrialSpikeTimes * median_gap_lengths(1,2)) + Notes.offsets(syllable(median_motif,2)) - Notes.onsets(syllable(median_motif,1));
%                 y_value = ones(length(TrialSpikeTimes),1) * i/5;
%                 Raster = [Raster; [TrialSpikeTimes y_value]];                
%             case 5
%                 TrialSpikeTimes = SpikeTimes(find((SpikeTimes > Notes.onsets(syllable(i,3))) & (SpikeTimes < Notes.offsets(syllable(i,3)))));
%                 TrialSpikeTimes = (TrialSpikeTimes - Notes.onsets(syllable(i,3)))/syllable_lengths(i,3);
%                 TrialSpikeTimes = (TrialSpikeTimes * median_syllable_lengths(1,3)) + Notes.onsets(syllable(median_motif,3)) - Notes.onsets(syllable(median_motif,1));
%                 y_value = ones(length(TrialSpikeTimes),1) * i/5;
%                 Raster = [Raster; [TrialSpikeTimes y_value]];                
%         end
    end
    edges = -0.2:bin_size:motif_durations(median_motif);
    temp_trial_pst = histc(Raster(:,1),edges);
    temp_trial_pst(end) = [];
    trial_pst = [trial_pst; temp_trial_pst'];
end
axes('position',[0.1 0.075 0.8 0.725]);
edges = -0.2:bin_size:motif_durations(median_motif);
pst = histc(Raster(:,1),edges);
pst(end) = [];
edges(end) = [];
pst = pst/size(motif_durations,1);
pst = pst/bin_size;
pst = pst/1000;
bar(edges,pst,'histc');
hold on;
Raster(:,2) = Raster(:,2) + max(pst) + 1;
plot(Raster(:,1),Raster(:,2),'w+');
marker_string = repmat('|',size(Raster,1),1);
text(Raster(:,1),Raster(:,2),marker_string,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',6);
axis([-0.2 motif_durations(median_motif) 0 (max(Raster(:,2)) + 0.5)]);
set(gca,'ytick',[0,max(pst)/2,max(pst)]);
ylabel('Firing rate (spikes / 100ms)','FontSize',14,'FontWeight','bold');
xlabel('Time (sec)','FontSize',14,'FontWeight','bold');

channel_string = strcat('obs0r');
pathname = pwd;
pathname = strcat(pathname,'/');
[RawData,Fs] = soundin_copy(pathname,RawDataFile,channel_string);
Times = 0:1/Fs:length(RawData)/Fs;
Times(end) = [];
Indices = (find((Times >= (Notes.onsets(syllable(median_motif,1)) - 0.2)) & (Times <= Notes.offsets(syllable(median_motif,length(motif))))));

axes('position',[0.1 0.825 0.8 0.15]);
%plot((Times(Indices) - Notes.onsets(syllable(median_motif,1))),RawData(Indices),'b');
filtsong = bandpass(RawData(Indices),Fs,300,10000);

nfft=round(Fs*8/1000);
nfft = 2^nextpow2(nfft);
spect_win = hanning(nfft);
noverlap = round(0.9*length(spect_win)); %number of overlapping points       
%now calculate spectrogram
[spect, freq, time] = specgram(filtsong, nfft, Fs, spect_win, noverlap);
idx_spect=scale_spect(spect);  %calculate index array for spectrogram
f_min = freq(1);
f_max = freq(length(freq));
freq_spect = [f_min, f_max];
time = time - 0.2;
t_min = time(1); %convert to ms
t_max = time(length(time)); %convert to ms
%adjust time axis for spectrogram offset (1/2 window duration in ms)
%t_min = t_min + 0.5*nfft*1000/Fs;  
%t_max = t_max + 0.5*nfft*1000/Fs;  
t_min = t_min;  
t_max = t_max;  
time_spect = [t_min, t_max];                

disp_idx_spect(idx_spect, time_spect, freq_spect, -55, ...
	0, 2, 'gray', 'classic');
axis([-0.2 t_max 300 10000]);
set(gca,'ytick',[1000 5000 10000]);
set(gca,'xtick',[]);
