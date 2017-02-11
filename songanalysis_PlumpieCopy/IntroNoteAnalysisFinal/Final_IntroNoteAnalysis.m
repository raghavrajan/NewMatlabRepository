function [] = Final_IntroNoteAnalysis(DataDir, RecFileDir, NoteFileList, NoteFileDir, FileType, ContinuousOrNot, OutputPattern, varargin)

PresentDir = pwd;

if (nargin > 9)
    SpikeFileDir = varargin{1};
    SpikeChanNo = varargin{2};
    ClusterNos = varargin{3};
end

%===================== General parameters used by program =================

InterBoutInterval = 2;

Fid = fopen(NoteFileList, 'r');
TempNoteFiles = textscan(Fid, '%s', 'DeLimiter', '\n');
TempNoteFiles = TempNoteFiles{1};
fclose(Fid);

for i = 1:length(TempNoteFiles),
    NoteFiles(i).name = [TempNoteFiles{i}, '.not.mat'];
end

%==========================================================================

[Notes, ANotes] = INALoadNoteFiles(DataDir, RecFileDir, NoteFiles, NoteFileDir, FileType, ContinuousOrNot);

InterBoutInt = cell(size(Notes.Onsets));
InterBoutInt(1:end) = {InterBoutInterval};

[BoutIndices, BoutLengths] = cellfun(@INAFindBouts, Notes.Onsets, Notes.Offsets, InterBoutInt, 'UniformOutput', 0);

Labels = [];
for i = 1:length(Notes.Labels),
    for j = 1:2:length(BoutIndices{i}),
        Labels = [Labels 'Q' Notes.Labels{i}(BoutIndices{i}(j):BoutIndices{i}(j+1)) 'q'];
    end
end
UniqueSylls = unique(Labels);

SyllIndex = find(UniqueSylls == 'q');

for i = 1:SyllIndex - 1,
    Matches = find(Labels == UniqueSylls(i));
    disp(['% of ', UniqueSylls(i), 's: ', num2str(length(Matches)/length(Labels) * 100)]);
    PostMatches = Labels(Matches+1);
    Sylls = unique(PostMatches);
    for j = 1:length(Sylls),
        disp([UniqueSylls(i), '--->', Sylls(j), ': prob. is ', num2str(length(find(PostMatches == Sylls(j)))/length(PostMatches))]);
    end
end

%Sylls = input('List the syllables that are part of the motif:', 's');
Sylls = MotifSylls;

Trans = zeros(length(Sylls));
for i = 1:length(Sylls)-1,
    Matches = find(Labels == Sylls(i));
    PostMatches = Labels(Matches+1);
    PMs = [];
    for j = 1:length(Sylls),
        PMs = [PMs PostMatches(find(PostMatches == Sylls(j)))];
    end
    for j = 1:length(Sylls),
        Trans(i,j) = (length(find(PMs == Sylls(j)))/length(PMs));
    end
end

% Labels = INAGenerateSyllSequence(Trans, Sylls, 10000);

%StartTrans = input('List the transition that represent the start of the motif: (for eg. type iaib if i to a and i to b represent motif starts):', 's');

if (nargin > 9)
    [SpikeData, ASpikeData] = INALoadSpikeData(ANotes.FileName, SpikeFileDir, ClusterNos, SpikeChanNo, ANotes.FileLength, ContinuousOrNot);
end

for i = 1:length(ANotes.Labels),
    TempMotifIndices = [];
    for j = 1:2:length(StartTrans),
        TempMotifIndices = [TempMotifIndices (strfind(ANotes.Labels{i}, StartTrans(j:j+1)) + 1)];
    end
    TempMotifIndices = sort(TempMotifIndices);
    MotifIndices{i} = TempMotifIndices;
end

ValidFiles = ~cellfun(@isempty, MotifIndices, 'UniformOutput', 1);

ANotes.Onsets(~ValidFiles) = [];
ANotes.Offsets(~ValidFiles) = [];
ANotes.Labels(~ValidFiles) = [];
ANotes.FileLength(~ValidFiles) = [];
ANotes.FileName(~ValidFiles) = [];
ANotes.FileType(~ValidFiles) = [];
ANotes.DataDir(~ValidFiles) = [];
ANotes.RecFileDir(~ValidFiles) = [];

if (nargin > 9)
    ASpikeData(~ValidFiles) = [];
end

if (isempty(strfind(ContinuousOrNot, 'continuous')))
    Notes.Onsets(~ValidFiles) = [];
    Notes.Offsets(~ValidFiles) = [];
    Notes.Labels(~ValidFiles) = [];
    Notes.FileLength(~ValidFiles) = [];
    Notes.FileName(~ValidFiles) = [];
    Notes.FileType(~ValidFiles) = [];
    Notes.DataDir(~ValidFiles) = [];
    Notes.RecFileDir(~ValidFiles) = [];
end

if (nargin > 9)
    SpikeData(~ValidFiles) = [];
end

MotifIndices(~ValidFiles) = [];

InterBoutInt = cell(size(ANotes.Onsets));
InterBoutInt(1:end) = {InterBoutInterval};

BoutIndices = [];
BoutLengths = [];

[BoutIndices, BoutLengths] = cellfun(@INAFindBouts, ANotes.Onsets, ANotes.Offsets, InterBoutInt, 'UniformOutput', 0);

[MotifIntroNoteIndices, NonMotifIntroNoteIndices] = cellfun(@INAFindIntroNotes, BoutIndices, MotifIndices, ANotes.Onsets, ANotes.Offsets, ANotes.Labels, InterBoutInt, ANotes.FileLength, 'UniformOutput', 0);

% for i = 1:length(Labels),
%     TempSimMotifIndices = [];
%     for j = 1:2:length(StartTrans),
%         TempSimMotifIndices = [TempSimMotifIndices (strfind(Labels{i}, StartTrans(j:j+1)) + 1)];
%     end
%     TempSimMotifIndices = sort(TempSimMotifIndices);
%     MotifSimIndices{i} = TempSimMotifIndices;
% end
% ValidFiles = ~cellfun(@isempty, MotifSimIndices, 'UniformOutput', 1);
% Labels(~ValidFiles) = [];
% MotifSimIndices(~ValidFiles) = [];
% 
% [SimMotifIntroNoteIndices] = cellfun(@INASimFindIntroNotes, Labels, MotifSimIndices, 'UniformOutput', 0);
% 

[MotifIntroNoteIndices, NonMotifIntroNoteIndices] = cellfun(@INACalcSAPFeatures, MotifIntroNoteIndices, NonMotifIntroNoteIndices, ANotes.Onsets, ANotes.Offsets, ANotes.FileName, ANotes.FileType, ANotes.DataDir, ANotes.RecFileDir, 'UniformOutput', 0);

% if (nargin > 9)
%     [SpontActivity] = cellfun(@INACalcSpontActivity, MotifIntroNoteIndices, ASpikeData);
%     [MotifIntroNoteIndices, NonMotifIntroNoteIndices] = cellfun(@INACalcEvokedActivity2, MotifIntroNoteIndices, NonMotifIntroNoteIndices, ASpikeData, ANotes.Labels, ANotes.FileLength, 'UniformOutput', 0);
% end

NoofIntroNotes = [];
VocalInterval = [];
Indices = [];
for i = 1:length(MotifIntroNoteIndices),
    if (~isempty(MotifIntroNoteIndices{i}))
        NoofIntroNotes = [NoofIntroNotes cell2mat(MotifIntroNoteIndices{i}.NoofINs)];
        TotalMotifs = length(cell2mat(MotifIntroNoteIndices{i}.NoofINs));
        Indices = [Indices [ones(1,TotalMotifs)*i; [1:1:TotalMotifs]]];
        VocalInterval = [VocalInterval cell2mat(MotifIntroNoteIndices{i}.VocalInterval)];
    end
end

VIIndices = find(VocalInterval <= InterBoutInterval);
MeanNoofIntroNotes(1) = mean(NoofIntroNotes(VIIndices));
StdNoofIntroNotes(1) = std(NoofIntroNotes(VIIndices));

VIIndices = find(VocalInterval > InterBoutInterval);
MeanNoofIntroNotes(2) = mean(NoofIntroNotes(VIIndices));
StdNoofIntroNotes(2) = std(NoofIntroNotes(VIIndices));

TempNoofIntroNotes = NoofIntroNotes(find(NoofIntroNotes > 0));
TempVocalInterval = VocalInterval(find(NoofIntroNotes > 0));

TVIIndices = find(TempVocalInterval <= InterBoutInterval);
TMeanNoofIntroNotes(1) = mean(TempNoofIntroNotes(TVIIndices));
TStdNoofIntroNotes(1) = std(TempNoofIntroNotes(TVIIndices));

TVIIndices = find(TempVocalInterval > InterBoutInterval);
TMeanNoofIntroNotes(2) = mean(TempNoofIntroNotes(TVIIndices));
TStdNoofIntroNotes(2) = std(TempNoofIntroNotes(TVIIndices));

% figure;
% subplot(2,1,1);
% Edges = 0:1:max(NoofIntroNotes);
% HistNoofIntroNotes = hist(NoofIntroNotes, Edges);
% HistBar = bar(Edges, HistNoofIntroNotes/sum(HistNoofIntroNotes) * 100);
% set(HistBar, 'FaceColor', 'w', 'EdgeColor', 'k');
% set(gcf, 'Color', 'w');
% xlabel('# of Intro notes', 'FontSize', 14, 'FontWeight', 'bold');
% ylabel('%', 'FontSize', 14, 'FontWeight', 'bold');
% axis auto;
% 
% subplot(2,1,2);
% plot(VocalInterval, NoofIntroNotes, 'k+');
% set(gca, 'XScale', 'log');
% set(gcf, 'Color', 'w');
% xlabel('Time since last vocalisation (sec)', 'FontSize', 14, 'FontWeight', 'bold');
% ylabel('# of Intro notes', 'FontSize', 14, 'FontWeight', 'bold');
% axis auto;
% temp2 = axis;
% temp2(3:4) = [-1 (Edges(end) + 1)];
% temp2(1) = 0.05;
% if (temp2(2) < 10)
%     temp2(2) = 10;
% end
% axis(temp2);
% 
% figure;
MaxINotes = max(NoofIntroNotes);
for i = 1:MaxINotes,
    Durs{i} = [];
    Intervals{i} = [];
    PG{i} = [];
    Amplitude{i} = [];
    Entropy{i} = [];
    Pitch{i} = [];
    FM{i} = [];
    AM{i} = [];
    MeanFreq{i} = [];
end

for i = 1:MaxINotes,
    INoteIndices = find(NoofIntroNotes >= i);
    for j = 1:length(INoteIndices),
        Durs{i} = [Durs{i}; MotifIntroNoteIndices{Indices(1, INoteIndices(j))}.Durs{Indices(2, INoteIndices(j))}(end - (i - 1))];
        Intervals{i} = [Intervals{i}; MotifIntroNoteIndices{Indices(1, INoteIndices(j))}.Intervals{Indices(2, INoteIndices(j))}(end - (i - 1))];
        PG{i} = [PG{i}; MotifIntroNoteIndices{Indices(1, INoteIndices(j))}.PG{Indices(2, INoteIndices(j))}(end - (i - 1))];
        Amplitude{i} = [Amplitude{i}; MotifIntroNoteIndices{Indices(1, INoteIndices(j))}.Amplitude{Indices(2, INoteIndices(j))}(end - (i - 1))];
        Entropy{i} = [Entropy{i}; MotifIntroNoteIndices{Indices(1, INoteIndices(j))}.Entropy{Indices(2, INoteIndices(j))}(end - (i - 1))];
        Pitch{i} = [Pitch{i}; MotifIntroNoteIndices{Indices(1, INoteIndices(j))}.Pitch{Indices(2, INoteIndices(j))}(end - (i - 1))];
        FM{i} = [FM{i}; MotifIntroNoteIndices{Indices(1, INoteIndices(j))}.FM{Indices(2, INoteIndices(j))}(end - (i - 1))];
        MeanFreq{i} = [MeanFreq{i}; MotifIntroNoteIndices{Indices(1, INoteIndices(j))}.MeanFreq{Indices(2, INoteIndices(j))}(end - (i - 1))];
        AM{i} = [AM{i}; MotifIntroNoteIndices{Indices(1, INoteIndices(j))}.AM{Indices(2, INoteIndices(j))}(end - (i - 1))];
    end
end

NMDurs = [];
NMIntervals = [];
NMPG = [];
NMAmplitude = [];
NMEntropy = [];
NMPitch = [];
NMFM = [];
NMAM = [];
NMMeanFreq = [];

for i = 1:length(NonMotifIntroNoteIndices),
    if (isfield(NonMotifIntroNoteIndices{i}, 'PG'))
        NMDurs = [NMDurs; NonMotifIntroNoteIndices{i}.Durs{1}];
        NMIntervals = [NMIntervals; NonMotifIntroNoteIndices{i}.Intervals{1}];
        NMPG = [NMPG; NonMotifIntroNoteIndices{i}.PG{1}'];
        NMAmplitude = [NMAmplitude; NonMotifIntroNoteIndices{i}.Amplitude{1}'];
        NMEntropy = [NMEntropy; NonMotifIntroNoteIndices{i}.Entropy{1}'];
        NMPitch = [NMPitch; NonMotifIntroNoteIndices{i}.Pitch{1}'];
        NMFM = [NMFM; NonMotifIntroNoteIndices{i}.FM{1}'];
        NMAM = [NMAM; NonMotifIntroNoteIndices{i}.AM{1}'];
        NMMeanFreq = [NMMeanFreq; NonMotifIntroNoteIndices{i}.MeanFreq{1}'];
    end
end

% FeatureFig = figure;
% set(gcf, 'Color', 'w');
% Features = [{'Durs'}, {'Intervals'}, {'Entropy'}, {'Amplitude'}, {'PG'}, {'AM'}, {'FM'}, {'MeanFreq'}];
% FeatureNames = [{'Duration (sec)'}, {'Interval to next syllable (sec)'}, {'Mean Wiener Entropy'}, {'Mean log amplitude'}, {'Mean Pitch goodness'}, {'Mean AM'}, {'Mean FM'}, {'Mean Frequency (hz)'}];
% FeatureVect = [];
% GroupVect = [];
% for i = 1:length(Features),
%     figure(FeatureFig);
%     subplot(4,2,i);
%     [TempFeatureVect, TempGroupVect] = eval(['INAPlotFeatures(', Features{i}, ',NM', Features{i}, ',''', FeatureNames{i}, ''')']);
%     FeatureVect = [FeatureVect TempFeatureVect];
%     GroupVect = [GroupVect TempGroupVect];
% end

cd(PresentDir);
if (exist('SpikeChanNo', 'var'))
    save([NoteFiles(1).name, '.Chan', num2str(SpikeChanNo), '.INoteAnalysisResults', OutputPattern, '.mat']);
else
    save([NoteFiles(1).name, '.INoteAnalysisResults', OutputPattern, '.mat']);
end

disp('Finished Analysis');
