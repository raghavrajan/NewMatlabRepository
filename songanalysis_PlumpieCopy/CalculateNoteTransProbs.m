function [BoutStats] = CalculateNoteTransProbs(RawDataDir, FileType, DirectoryName, NoteFiles, Motif, Motif2)

InterBoutInterval = 0.5; % in seconds

if ispc
    if (RawDataDir(end) ~= '\')
        RawDataDir(end+1) = '\';
    end
else
    if (RawDataDir(end) ~= '/')
        RawDataDir(end+1) = '/';
    end
end
      
cd(DirectoryName);

BoutNo = 1;

for NoteFileNo = 1:length(NoteFiles),
    Notes = load(NoteFiles(NoteFileNo).name);
    disp(NoteFiles(NoteFileNo).name);
    Notes.onsets = Notes.onsets/1000;
    Notes.offsets = Notes.offsets/1000;
    BoutIndices = find((Notes.onsets(2:end) - Notes.offsets(1:end-1)) > InterBoutInterval);
    Bouts = [];
    if (~isempty(BoutIndices))
        if (Notes.onsets(1) > InterBoutInterval)
            Bouts(1,:) = [1 BoutIndices(1)];
            for i = 1:(length(BoutIndices)-1),
                Bouts(i+1,:) = [(BoutIndices(i)+1) (BoutIndices(i+1))];
            end
        else
            for i = 1:(length(BoutIndices)-1),
                Bouts(i,:) = [(BoutIndices(i)+1) (BoutIndices(i+1))];
            end
        end
        Bouts(end+1,:) = [(BoutIndices(end) + 1) length(Notes.labels)];
    else
        if (Notes.onsets(1) > InterBoutInterval)
            Bouts = [1 length(Notes.labels)];
        end
    end
    for i = 1:size(Bouts,1),
        Motifs = union(strfind(Notes.labels(Bouts(i,1):Bouts(i,2)), Motif), strfind(Notes.labels(Bouts(i,1):Bouts(i,2)), Motif2));
        Motifs = Motifs + Bouts(i,1) - 1;
        if (isempty(Motifs))
            BoutStats.Motifs.NoofMotifs(BoutNo) = 0;
            continue;
        else
            IntroNotes = find(Notes.labels(Bouts(i,1):Motifs(1)) == 'i');
            IntroNotes = IntroNotes + Bouts(i,1) - 1;
            if (isempty(IntroNotes))
                BoutStats.IntroNotes.NoofIntroNotes(BoutNo) = 0;
            else
                BoutStats.IntroNotes.NoofIntroNotes(BoutNo) = length(IntroNotes);
            end
            BoutStats.Motifs.NoofMotifs(BoutNo) = length(Motifs);
        end
        BoutNo = BoutNo + 1;
    end
end

BoutStats.IntroNotes.Durations = ones(BoutNo-1, max(BoutStats.IntroNotes.NoofIntroNotes))*-100;
BoutStats.IntroNotes.Intervals = ones(BoutNo-1, max(BoutStats.IntroNotes.NoofIntroNotes))*-100;
BoutStats.IntroNotes.PG = ones(BoutNo-1, max(BoutStats.IntroNotes.NoofIntroNotes))*-100;
BoutStats.IntroNotes.Entropy = ones(BoutNo-1, max(BoutStats.IntroNotes.NoofIntroNotes))*100;
BoutStats.IntroNotes.Amplitude = ones(BoutNo-1, max(BoutStats.IntroNotes.NoofIntroNotes))*-100;
BoutStats.IntroNotes.FM = ones(BoutNo-1, max(BoutStats.IntroNotes.NoofIntroNotes))*-100;
BoutStats.IntroNotes.Pitch = ones(BoutNo-1, max(BoutStats.IntroNotes.NoofIntroNotes))*-100;
BoutStats.IntroNotes.MeanFreq = ones(BoutNo-1, max(BoutStats.IntroNotes.NoofIntroNotes))*-100;
BoutStats.IntroNotes.AM = ones(BoutNo-1, max(BoutStats.IntroNotes.NoofIntroNotes))*-100;

BoutNo = 1;

% BoutStats.Motifs.Durations = ones(BoutNo-1, max(BoutStats.Motifs.NoofMotifs))*-100;
% BoutStats.Motifs.Intervals = ones(BoutNo-1, max(BoutStats.Motifs.NoofMotifs))*-100;
% BoutStats.Motifs.PG = ones(BoutNo-1, max(BoutStats.Motifs.NoofMotifs))*-100;
% BoutStats.Motifs.Entropy = ones(BoutNo-1, max(BoutStats.Motifs.NoofMotifs))*100;
% BoutStats.Motifs.Amplitude = ones(BoutNo-1, max(BoutStats.Motifs.NoofMotifs))*-100;
% BoutStats.Motifs.FM = ones(BoutNo-1, max(BoutStats.Motifs.NoofMotifs))*-100;
% BoutStats.Motifs.Pitch = ones(BoutNo-1, max(BoutStats.Motifs.NoofMotifs))*-100;
% BoutStats.Motifs.MeanFreq = ones(BoutNo-1, max(BoutStats.Motifs.NoofMotifs))*-100;
% BoutStats.Motifs.AM = ones(BoutNo-1, max(BoutStats.Motifs.NoofMotifs))*-100;

for NoteFileNo = 1:length(NoteFiles),
    Index = strfind(NoteFiles(NoteFileNo).name, '.not.mat');
    SongFile = NoteFiles(NoteFileNo).name(1:Index-1);
    if (strfind(FileType, 'okrank'))
        [Song, Fs] = ReadOKrankData(RawDataDir, SongFile, 1);
    else
        if (strfind(FileType, 'wav'))
            cd(RawDataDir);
            [Song, Fs] = wavread(SongFile);
            cd(DirectoryName);
        else
            if (strfind(FileType, 'obs'));
                [Song, Fs] = soundin_copy(RawDataDir, SongFile, 'obs0r');
                Song = Song * 5/32768;
            end
        end
    end
    Time = (0:1:(length(Song)-1))/Fs;
    [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude , m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight] = deriv(Song, Fs);
    T = linspace(Time(1), Time(end), length(m_Entropy));
    
    Notes = load(NoteFiles(NoteFileNo).name);
    disp(NoteFiles(NoteFileNo).name);
    Notes.onsets = Notes.onsets/1000;
    Notes.offsets = Notes.offsets/1000;
    BoutIndices = find((Notes.onsets(2:end) - Notes.offsets(1:end-1)) > InterBoutInterval);
    Bouts = [];
    if (~isempty(BoutIndices))
        if (Notes.onsets(1) > InterBoutInterval)
            Bouts(1,:) = [1 BoutIndices(1)];
            for i = 1:(length(BoutIndices)-1),
                Bouts(i+1,:) = [(BoutIndices(i)+1) (BoutIndices(i+1))];
            end
        else
            for i = 1:(length(BoutIndices)-1),
                Bouts(i,:) = [(BoutIndices(i)+1) (BoutIndices(i+1))];
            end
        end
        Bouts(end+1,:) = [(BoutIndices(end) + 1) length(Notes.labels)];
    else
        if (Notes.onsets(1) > InterBoutInterval)
            Bouts = [1 length(Notes.labels)];
        end
    end
    for i = 1:size(Bouts,1),
        Motifs = union(strfind(Notes.labels(Bouts(i,1):Bouts(i,2)), Motif), strfind(Notes.labels(Bouts(i,1):Bouts(i,2)), Motif2));
        Motifs = Motifs + Bouts(i,1) - 1;
        if (isempty(Motifs))
            continue;
        else
            IntroNotes = find(Notes.labels(Bouts(i,1):Motifs(1)) == 'i');
            IntroNotes = IntroNotes + Bouts(i,1) - 1;
            if (~isempty(IntroNotes))
                IntroNotes = sort(IntroNotes, 'descend');
                BoutStats.IntroNotes.Durations(BoutNo, 1:length(IntroNotes)) = Notes.offsets(IntroNotes) - Notes.onsets(IntroNotes);
                BoutStats.IntroNotes.Intervals(BoutNo, 1:length(IntroNotes)) = Notes.onsets(IntroNotes + 1) - Notes.offsets(IntroNotes);
                for j = 1:length(IntroNotes),
                    StartIndex = find(T >= Notes.onsets(IntroNotes(j)), 1, 'first');
                    EndIndex = find(T >= Notes.offsets(IntroNotes(j)), 1, 'first');
                    
                    BoutStats.IntroNotes.PG(BoutNo, j) = mean(m_PitchGoodness(StartIndex:EndIndex));
                    BoutStats.IntroNotes.Entropy(BoutNo, j) = mean(m_Entropy(StartIndex:EndIndex));
                    BoutStats.IntroNotes.Amplitude(BoutNo, j) = mean(m_amplitude(StartIndex:EndIndex));
                    BoutStats.IntroNotes.FM(BoutNo, j) = mean(m_FM(StartIndex:EndIndex));
                    BoutStats.IntroNotes.Pitch(BoutNo, j) = mean(Pitch_chose(StartIndex:EndIndex));
                    BoutStats.IntroNotes.MeanFreq(BoutNo, j) = mean(m_Freq(StartIndex:EndIndex));
                    BoutStats.IntroNotes.AM(BoutNo, j) = mean(m_AM(StartIndex:EndIndex));
                end
            end
%             if (length(Motifs) > 1)
%                 for j = 1:length(Motifs)-1,
%                     BoutStats.Motifs{BoutNo}.MotifPos{j}.Durations(1,:) = Notes.offsets(Motifs(j):(Motifs(j+1)-1)) - Notes.onsets(Motifs(j):(Motifs(j+1)-1));
%                     BoutStats.Motifs{BoutNo}.MotifPos{j}.Intervals(1,:) = Notes.onsets((Motifs(j)+1):(Motifs(j+1)-1)) - Notes.offsets(Motifs(j):(Motifs(j+1)-1-1));
%                 end
%                 BoutStats.Motifs{BoutNo}.MotifPos{j+1}.Durations(1,:) = Notes.offsets(Motifs(end):(Motifs(end)+length(Motif)-1)) - Notes.onsets(Motifs(end):(Motifs(end)+length(Motif)-1));
%                 BoutStats.Motifs{BoutNo}.MotifPos{j+1}.Intervals(1,:) = Notes.onsets((Motifs(end)+1):(Motifs(end)+length(Motif)-1)) - Notes.offsets(Motifs(end):(Motifs(end)+length(Motif)-1-1));
%             else
%                 BoutStats.Motifs{BoutNo}.MotifPos{1}.Durations(1,:) = Notes.offsets(Motifs(end):(Motifs(end)+length(Motif)-1)) - Notes.onsets(Motifs(end):(Motifs(end)+length(Motif)-1));
%                 BoutStats.Motifs{BoutNo}.MotifPos{1}.Intervals(1,:) = Notes.onsets((Motifs(end)+1):(Motifs(end)+length(Motif)-1)) - Notes.offsets(Motifs(end):(Motifs(end)+length(Motif)-1-1));
%             end
        end
        BoutNo = BoutNo + 1;
    end
end

AllFeatures = [];
AllFeatureIndices = [];
%for i = 1:size(BoutStats.IntroNotes.Durations, 2),
for i = 1:5,
    Indices = find(BoutStats.IntroNotes.Durations(:,i) > 0);
    AllFeatures = [AllFeatures; [BoutStats.IntroNotes.Durations(Indices,i) BoutStats.IntroNotes.PG(Indices,i) BoutStats.IntroNotes.Entropy(Indices,i) BoutStats.IntroNotes.Amplitude(Indices,i) BoutStats.IntroNotes.FM(Indices,i) BoutStats.IntroNotes.Pitch(Indices,i) BoutStats.IntroNotes.MeanFreq(Indices,i) BoutStats.IntroNotes.AM(Indices,i)]];
    AllFeatureIndices = [AllFeatureIndices; ones(length(Indices), 1)*i];
end

AllFeatures = zscore(AllFeatures);
[coeff, score, latent] = princomp(AllFeatures);
disp(num2str(latent*100/sum(latent)));

BoutStats.IntroNotes.AllFeatures = AllFeatures;
BoutStats.IntroNotes.AllFeatureIndices = AllFeatureIndices;

% PC1-2, acoustic features figure
figure;
set(gcf, 'Color', 'w');
Cols = ['rbcmk'];
%for i = 1:size(BoutStats.IntroNotes.Durations, 2),
for i = 1:5,
    Indices = find(AllFeatureIndices == i);
    plot(score(Indices,1), score(Indices,2), [Cols(mod(i,length(Cols))+1), '+'], 'MarkerSize', 4);
    BoutStats.IntroNotes.PCMeans(i,:) = mean(score(Indices,1:6));
    hold on;
end
axis tight;
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
xlabel('PC 1', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('PC 2', 'FontSize', 14, 'FontWeight', 'bold');
title('Acoustic features of intro notes', 'FontSize', 14, 'FontWeight', 'bold');

% PC1-2, acoustic features figure
figure;
set(gcf, 'Color', 'w');
Cols = ['rbcmk'];
%for i = 1:size(BoutStats.IntroNotes.Durations, 2),
for i = 1:5,
    Indices = find(AllFeatureIndices == i);
    [x, y, z] = ellipsoid(mean(score(Indices,1)), mean(score(Indices,2)), mean(score(Indices,3)), std(score(Indices,1)), std(score(Indices,2)), std(score(Indices, 3)));
    Surface(i) = surf(x, y, z);
    set(Surface(i), 'FaceColor', Cols(mod(i,length(Cols))+1));
    hold on;
end
axis tight;
legend(num2str([1:1:size(BoutStats.IntroNotes.Durations, 2)]'))
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
xlabel('PC 1', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('PC 2', 'FontSize', 14, 'FontWeight', 'bold');
zlabel('PC 3', 'FontSize', 14, 'FontWeight', 'bold');
title('Acoustic features of intro notes - mean + 1sd ellipsoids', 'FontSize', 14, 'FontWeight', 'bold');

% PC1-2, acoustic features figure
figure;
set(gcf, 'Color', 'w');
Cols = ['rbcmk'];
%for i = 1:size(BoutStats.IntroNotes.Durations, 2),
for i = 1:5,
    Indices = find(AllFeatureIndices == i);
    [x, y, z] = ellipsoid(mean(score(Indices,1)), mean(score(Indices,2)), mean(score(Indices,3)), std(score(Indices,1)), std(score(Indices,2)), std(score(Indices, 3)));
    for j = 1:3,
        subplot(2,2,j);
        Surface(i) = surf(x, y, z);
        set(Surface(i), 'FaceColor', Cols(mod(i,length(Cols))+1));
        hold on;
        axis tight;
        set(gca, 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('PC 1', 'FontSize', 14, 'FontWeight', 'bold');
        ylabel('PC 2', 'FontSize', 14, 'FontWeight', 'bold');
        zlabel('PC 3', 'FontSize', 14, 'FontWeight', 'bold');
        if (j == 1)
            view(0,90);
        else
            if (j == 2)
                view(0,180);
            else
                view(90,180);
            end
        end
        legend(num2str([1:1:size(BoutStats.IntroNotes.Durations, 2)]'))
    end
end
title('Acoustic features of intro notes - mean + 1sd ellipsoids', 'FontSize', 14, 'FontWeight', 'bold');

% No of Intro notes figure
figure;
set(gcf, 'Color', 'w');
MaxINotes = max(BoutStats.IntroNotes.NoofIntroNotes);
NoofINotesHist = hist(BoutStats.IntroNotes.NoofIntroNotes, [0:1:MaxINotes])*100/(BoutNo-1);
NoofINotesBar = bar([0:1:MaxINotes], NoofINotesHist);
set(NoofINotesBar, 'EdgeColor', 'k', 'FaceColor', 'w');
axis([0 (MaxINotes+1) 0 1.2*max(NoofINotesHist)]);
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
xlabel('No of intro notes', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('%', 'FontSize', 14, 'FontWeight', 'bold');
title('No of intro notes per bout', 'FontSize', 14, 'FontWeight', 'bold');

% Duration of different intro notes figure
figure;
set(gcf, 'Color', 'w');
for i = 1:size(BoutStats.IntroNotes.Durations, 2);
    MeanDurBar(i) = bar(i, median(BoutStats.IntroNotes.Durations(find(BoutStats.IntroNotes.Durations(:,i) > 0), i)));
    set(MeanDurBar(i), 'EdgeColor', 'k', 'FaceColor', 'w');
    hold on;
    plot(ones(size(BoutStats.IntroNotes.Durations(find(BoutStats.IntroNotes.Durations(:,i) > 0), i)))*i, BoutStats.IntroNotes.Durations(find(BoutStats.IntroNotes.Durations(:,i) > 0), i), 'k+', 'MarkerSize', 4);
end
axis tight;
temp = axis;
temp(1:2) = [0 i+1];
axis(temp);
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Position of intro note (1 is the closest to the motif)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Duration (sec)', 'FontSize', 14, 'FontWeight', 'bold');
title('Duration of intro notes', 'FontSize', 14, 'FontWeight', 'bold');

% Inter syllable interval of different intro notes figure
figure;
set(gcf, 'Color', 'w');
for i = 1:size(BoutStats.IntroNotes.Intervals, 2);
    MeanDurBar(i) = bar(i, median(BoutStats.IntroNotes.Intervals(find(BoutStats.IntroNotes.Intervals(:,i) > 0), i)));
    set(MeanDurBar(i), 'EdgeColor', 'k', 'FaceColor', 'w');
    hold on;
    plot(ones(size(BoutStats.IntroNotes.Intervals(find(BoutStats.IntroNotes.Intervals(:,i) > 0), i)))*i, BoutStats.IntroNotes.Intervals(find(BoutStats.IntroNotes.Intervals(:,i) > 0), i), 'k+', 'MarkerSize', 4);
end
axis tight;
temp = axis;
temp(1:2) = [0 i+1];
axis(temp);
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Position of intro note (1 is the closest to the motif)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Interval (sec)', 'FontSize', 14, 'FontWeight', 'bold');
title('Inter syllable interval', 'FontSize', 14, 'FontWeight', 'bold');

disp('Finished analysing notes');