function [] = AnalyzeStimulusResponse(FileList, StimulusLabel, PreStimTime, PostStimTime)

% =========================================================================
% First load up all the note files and get information about syllables,
% their labels, onsets and offsets.
% Data is recorded in continuous fashion, so the onsets and offsets can be
% joined together.
% For some syllables, the first part is at the end of one file and the next
% part is at the beginning of the next file. The two parts have been
% labelled with the same capital letter label. So, we have to check if two
% consecutive labels are in capital letters and offset time of first
% syllable is same as onset time of next syllable. Such syllables will then
% be merged.
% =========================================================================

% First load up all notefiles.

Fid = fopen(FileList, 'r');
Files = textscan(Fid, '%s', 'DeLimiter', '\n');
Files = Files{1};
fclose(Fid);

Labels = [];
Onsets = [];
Offsets = [];

PresentDir = pwd;

StartTime = 0;
for i = 1:length(Files),
    % First load note file
    cd ASSLNoteFiles;
    Temp = load([Files{i}, '.not.mat']);
    Labels = [Labels; Temp.labels(:)];
    Onsets = [Onsets; (StartTime + Temp.onsets(:)/1000)]; % Divide by 1000 is to convert from ms to s
    Offsets = [Offsets; (StartTime + Temp.offsets(:)/1000)]; % Divide by 1000 is to convert from ms to s
    
    % Now get total duration of data file
    cd(PresentDir);
    [Data, Fs] = wavread(Files{i});
    DataFileDur = length(Data)/Fs;
    StartTime = StartTime + DataFileDur;
end

% =========================================================================
% Now to find syllables that are split across two files.
InterSyllGaps = Onsets(2:end) - Offsets(1:end-1);
SplitSylls = find(InterSyllGaps < 0.0001);
SyllsToBeChanged = [];
for i = 1:length(SplitSylls),
    if (Labels(SplitSylls(i)) == Labels(SplitSylls(i) + 1))
        SyllsToBeChanged = [SyllsToBeChanged; SplitSylls(i)];
    end
end

Labels(SyllsToBeChanged + 1) = [];
Onsets(SyllsToBeChanged + 1) = [];
Offsets(SyllsToBeChanged) = [];
% =========================================================================

% =========================================================================
% Now remove all the '0' labels as these are just noise
NoiseSylls = find(Labels == '0');
Labels(NoiseSylls) = [];
Onsets(NoiseSylls) = [];
Offsets(NoiseSylls) = [];
% =========================================================================

% =========================================================================
% Now, find all syllables in the PreTime and PostTime relative to the
% StimulusLabel provided as input

% Things to be analysed
% 1) Plot vocalizations in the pre-stim and post-stim period for each
% stimulus trial. Each trial is one row.
% 2) Plot vocalizations in the pre-stim and post-stim period for each
% stimulus trial. Each trial is one row and the duration of the response is
% also plotted.
% 3) No of post stim period vocalizations - No of pre stim period
% vocalizations
% 4) Latency of first vocalization
% 5) Response reliability - % of trials when there is a response in the
% first 1s

Stimuli = find(Labels == StimulusLabel);
StimulusOnsets = Onsets(Stimuli);
StimulusOffsets = Offsets(Stimuli);

Onsets(Stimuli) = [];
Offsets(Stimuli) = [];
Labels(Stimuli) = [];

% === Analysis 1 ==========================================================
figure;
hold on;

for i = 1:length(Stimuli),
    StimulusOnsetTime = StimulusOnsets(i);
    StimulusOffsetTime = StimulusOffsets(i);
    
    TrialVocalizations = find((Onsets >= (StimulusOnsetTime - PreStimTime)) & (Onsets <= (StimulusOnsetTime + PostStimTime)));
    plot(Onsets(TrialVocalizations) - StimulusOnsetTime, ones(size(TrialVocalizations))*i, 'ks', 'MarkerSize', 6);
end

axis([(-PreStimTime - 0.01) (PostStimTime + 0.01) 0 (length(Stimuli) + 1)]);
plot([0 0], [0 (length(Stimuli) + 1)], 'k--', 'LineWidth', 2);

% === Analysis 2 ==========================================================
figure;
hold on;

for i = 1:length(Stimuli),
    StimulusOnsetTime = StimulusOnsets(i);
    StimulusOffsetTime = StimulusOffsets(i);
    
    TrialVocalizations = find((Onsets >= (StimulusOnsetTime - PreStimTime)) & (Onsets <= (StimulusOnsetTime + PostStimTime)));
    for j = 1:length(TrialVocalizations),
        plot([(Onsets(TrialVocalizations(j)) - StimulusOnsetTime) (Offsets(TrialVocalizations(j)) - StimulusOnsetTime)], [i i], 'k', 'LineWidth', 2);
    end
end

axis([(-PreStimTime - 0.01) (PostStimTime + 0.01) 0 (length(Stimuli) + 1)]);
plot([0 0], [0 (length(Stimuli) + 1)], 'k--', 'LineWidth', 2);


    