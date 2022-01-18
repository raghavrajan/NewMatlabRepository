function [PSTHTime, MotifStartActivity, BoutStartActivity, SongStartActivity, ValidSongBouts] = CalculateGaussianSmoothedMotifOnsetActivity(SURecordingDetails, GaussWidth)

Files = SURecordingDetails.SongFileNames;
DataDir = SURecordingDetails.DataDirectory;
FileType = SURecordingDetails.FileType;
Motif = SURecordingDetails.CommonMotifs{1};
INs = SURecordingDetails.INLabels;
if (strfind(INs, 'NA'))
    INs = [];
end

ArtefactThreshold = 0.75; % in mv

[RawData, Fs] = GetData(DataDir, Files{1}, FileType, 0);

Width = GaussWidth; % in seconds - Gaussian width or std
GaussianLen = 3; % length of gaussian

XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));

% if (SURecordingDetails.Tetrode == 1)
%     SpikeChanNo = SURecordingDetails.SpikeChanNo + [0 1 2 3];
% else
    SpikeChanNo = SURecordingDetails.SpikeChanNo;
% end

ValidSongBouts = find((SURecordingDetails.Bouts(:,7) == 1) & (SURecordingDetails.Bouts(:,8) > 0) & (SURecordingDetails.Bouts(:,9) > 1));


SongStartActivity = [];
MotifStartActivity = [];
BoutStartActivity = [];
PreTime = 2; % in seconds
PostTime = 1; % in seconds

Index = 1;

% for i = ValidSongBouts(:)',
%     NoteInfo = load(fullfile(DataDir, 'ASSLNoteFiles', [Files{SURecordingDetails.Bouts(i,3)}, '.not.mat']));
%     Labels = NoteInfo.labels(SURecordingDetails.Bouts(i,1):SURecordingDetails.Bouts(i,2));
%     Onsets = NoteInfo.onsets(SURecordingDetails.Bouts(i,1):SURecordingDetails.Bouts(i,2));
%     Offsets = NoteInfo.offsets(SURecordingDetails.Bouts(i,1):SURecordingDetails.Bouts(i,2));
%     Motifs = strfind(Labels(:)', Motif);
%     
%     if (~isempty(Motifs))
%         MotifOnsetTime = Onsets(Motifs(1))/1000;
%         MotifOffsetTime = Offsets(Motifs(1))/1000; 
%         for j = 1:length(SpikeChanNo)',
%             [SpikeData, Fs] = GetData(DataDir, Files{SURecordingDetails.Bouts(i,3)}, FileType, SpikeChanNo(j));
%             MotifSpikeData = SpikeData(round(Fs*(MotifOnsetTime - PreTime)):round(Fs * (MotifOffsetTime + PreTime)));
%             % Remove data that is an outlier as an artefact
%             Threshold = [(prctile(MotifSpikeData, 25) - 3*iqr(MotifSpikeData)) (prctile(MotifSpikeData, 75) + 3*iqr(MotifSpikeData))]; 
%             MotifSpikeData(find((MotifSpikeData < Threshold(1)) | (MotifSpikeData > Threshold(2)))) = [];
%             RMS(i,j) = sqrt(mean(MotifSpikeData.^2));
%             MaxSpike(i,j) = max(abs(MotifSpikeData));
%         end
%     end
% end
% 
% [MaxVal, BestChan] = max(mean(MaxSpike));

SmoothedTraceFs = 1000; % 1000 Hz = 1ms resolution
PSTHTime = -PreTime:1/SmoothedTraceFs:PostTime;

TrialsToExclude = [];
for i = ValidSongBouts(:)',
    NoteInfo = load(fullfile(DataDir, 'ASSLNoteFiles', [Files{SURecordingDetails.Bouts(i,3)}, '.not.mat']));
    Labels = NoteInfo.labels(SURecordingDetails.Bouts(i,1):SURecordingDetails.Bouts(i,2));
    Onsets = NoteInfo.onsets(SURecordingDetails.Bouts(i,1):SURecordingDetails.Bouts(i,2));
    Motifs = strfind(Labels(:)', Motif);
    
    if (~isempty(Motifs))
        [SpikeData, Fs] = GetData(DataDir, Files{SURecordingDetails.Bouts(i,3)}, FileType, SpikeChanNo);
        SpikeData = SpikeData * 1000 / 10000; % to convert it to mV by first dividing by the gain (10000) and then multiplying by 1000
        [b, a] = butter(2, [300/(Fs/2) 6000/(Fs/2)]);
        SpikeData = filtfilt(b, a, SpikeData);
        
        SmoothedSpikeData = conv((SpikeData.^2), GaussWin, 'same');
        SmoothedSpikeData = spline((1:1:length(SpikeData))/Fs, SmoothedSpikeData, (1/Fs:1/SmoothedTraceFs:length(SpikeData)/Fs));
        
        [SongData, Fs] = GetData(DataDir, Files{SURecordingDetails.Bouts(i,3)}, FileType, 1);
        SmoothedSongData = ASSLCalculateLogAmplitudeAronovFee(SongData, Fs, (1:1:length(SongData))/Fs, 5, []);
        SmoothedSongData = spline((1:1:length(SpikeData))/Fs, SmoothedSongData, (1/Fs:1/SmoothedTraceFs:length(SongData)/Fs));
        
        % Now check that there are no syllables before the first syllable
        % of the motif and that there is > 2s of data before
        if (isempty(INs))
            if (Motifs(1) == 1)
                MotifOnsetTime = Onsets(Motifs(1))/1000;
                
                if (~isempty(find(abs(SpikeData(ceil((MotifOnsetTime - PreTime) * Fs):ceil((MotifOnsetTime + PostTime) * Fs))) >= ArtefactThreshold)))
                    TrialsToExclude(end+1) = i;
                    continue;
                end
        
                
                MotifStartActivity(end+1,:) = SmoothedSpikeData(ceil((MotifOnsetTime - PreTime) * SmoothedTraceFs):ceil((MotifOnsetTime - PreTime) * SmoothedTraceFs) + length(PSTHTime) - 1);
                SongStartActivity(end+1,:) = SmoothedSongData(ceil((MotifOnsetTime - PreTime) * SmoothedTraceFs):ceil((MotifOnsetTime - PreTime) * SmoothedTraceFs) + length(PSTHTime) - 1);
            end
        else
            TotalINs = 0;
            for IN_Index = 1:length(INs),
                TotalINs = TotalINs + length(find(Labels(1:Motifs(1)) == INs(IN_Index)));
            end
            if (TotalINs == (Motifs(1) - 1))
                MotifOnsetTime = Onsets(Motifs(1))/1000;
                BoutOnsetTime = Onsets(1)/1000;
                
                if ((~isempty(find(abs(SpikeData(ceil((MotifOnsetTime - PreTime) * Fs):ceil((MotifOnsetTime + PostTime) * Fs))) >= ArtefactThreshold))) || (~isempty(find(abs(SpikeData(ceil((BoutOnsetTime - PreTime) * Fs):ceil((BoutOnsetTime + PostTime) * Fs))) >= ArtefactThreshold))))
                    TrialsToExclude(end+1) = i;
                    continue;
                end
                
                MotifStartActivity(end+1,:) = SmoothedSpikeData(ceil((MotifOnsetTime - PreTime) * SmoothedTraceFs):ceil((MotifOnsetTime - PreTime) * SmoothedTraceFs) + length(PSTHTime) - 1);
                BoutStartActivity(end+1,:) = SmoothedSpikeData(ceil((BoutOnsetTime - PreTime) * SmoothedTraceFs):ceil((BoutOnsetTime - PreTime) * SmoothedTraceFs) + length(PSTHTime) - 1);
                SongStartActivity(end+1,:) = SmoothedSongData(ceil((MotifOnsetTime - PreTime) * SmoothedTraceFs):ceil((MotifOnsetTime - PreTime) * SmoothedTraceFs) + length(PSTHTime) - 1);
            end
        end
    end
    Index = Index + 1;
end

% Now to exclude trials that have transients based on trials that have any
% individual value greater than 75th percentile of all MotifStartActivity +
% 5 * iqr(MotifStartActivity)
% ALready done with a fixed artefact threshold above

% Threshold = prctile(MotifStartActivity(:), 75) + (11 * iqr(MotifStartActivity(:)));
% TrialsToExclude = [];
% for i = 1:size(MotifStartActivity, 1),
%     Artefacts = find(MotifStartActivity(i,:) >= Threshold);
%     if (~isempty(Artefacts))
%         TrialsToExclude(end+1) = i;
%     end
% end

ValidSongBouts = setdiff(ValidSongBouts, TrialsToExclude);
% MotifStartActivity(TrialsToExclude, :) = [];
% SongStartActivity(TrialsToExclude, :) = [];

disp(['Skipped ', num2str(length(TrialsToExclude)), ' bouts with artefacts : ', num2str(length(TrialsToExclude)), '/', num2str(length(TrialsToExclude) + length(ValidSongBouts)), ' bouts (', num2str(100 * length(TrialsToExclude)/(length(TrialsToExclude) + length(ValidSongBouts))), ' %)']);
PSTHTime = linspace(-PreTime, PostTime, size(MotifStartActivity, 2));

disp('Finished');