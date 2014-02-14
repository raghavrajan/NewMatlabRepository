function [] = PlotGoodBadSyllTransitions(RawDataDir, FileList, FileType, Motif)

PresentDir = pwd;

Fid = fopen(FileList, 'r');
TempFiles = textscan(Fid, '%s', 'DeLimiter', '\n');
NoteFiles = TempFiles{1};
fclose(Fid);

for i = 1:length(NoteFiles),
    SlashIndex = find((NoteFiles{i} == '/') | (NoteFiles{i} == '\'));
    if (~isempty(SlashIndex))
        NoteFiles{i} = NoteFiles{i}(SlashIndex(end)+1:end);
    end
end

if (strfind(FileType, 'obs'))
    SongChanNo = 'obs0r';
else
    SongChanNo = 1;
end

Colours = ['rgbcmk'];
Labels = [];
AllFeats = [];
AllAmplitudes = [];
AllTimes = [];
AllSAPAmplitudes = [];
AllSAPEntropies = [];
AllSAPTimes = [];

for i = 1:length(NoteFiles),
    if (~exist([NoteFiles{i}, '.not.mat'], 'file'))
        continue;
    end
    
    Notes = load([NoteFiles{i}, '.not.mat']);
    SongFile = NoteFiles{i};
    
    [RawData, Fs] = ASSLGetRawData([RawDataDir, '/'], SongFile, FileType, SongChanNo);

    Time = (1:1:length(RawData))/Fs;

    [LogAmplitude] = ASSLCalculateLogAmplitude(RawData, Fs, Time, 8, 0.9);
    LogAmplitudeTime = linspace(Time(1), Time(end), length(LogAmplitude));
    %[Feats, EntireAmplitude] = CalculateSAPFeatsIndividualFile(RawDataDir, PresentDir, SongFile, FileType);
    [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude , m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight] = deriv(RawData, Fs);
    SAPFeature_Time = linspace(Time(1),Time(end), length(m_amplitude));
    
    Begin{1} = NaN;
    End{1} = NaN;
    clear EntireAmplitude EntireTime EntireSAPTime EntireSAPAmplitude EntireSAPEntropy;
    for j = 1:length(Notes.onsets),
        if (Notes.onsets(j)/1000 <= 0.5)
            EntireAmplitude{j} = NaN;
            EntireTime{j} = NaN;
            EntireSAPTime{j} = NaN;
            EntireSAPAmplitude{j} = NaN;
            EntireSAPEntropy{j} = NaN;
            continue;
        end
        if (Notes.onsets(j)/1000 >= (Time(end) - 0.5))
            EntireAmplitude{j} = NaN;
            EntireTime{j} = NaN;
            EntireSAPTime{j} = NaN;
            EntireSAPAmplitude{j} = NaN;
            EntireSAPEntropy{j} = NaN;
            continue;
        end
        StartIndex = find(LogAmplitudeTime > (Notes.onsets(j)/1000 - 0.5), 1, 'first');
        EndIndex = find(LogAmplitudeTime > (Notes.onsets(j)/1000 + 0.5), 1, 'first');
        EntireAmplitude{j} = LogAmplitude(StartIndex:EndIndex);
        EntireTime{j} = LogAmplitudeTime(StartIndex:EndIndex) - Notes.offsets(j)/1000;

        StartIndex = find(SAPFeature_Time > (Notes.onsets(j)/1000 - 0.5), 1, 'first');
        EndIndex = find(SAPFeature_Time > (Notes.onsets(j)/1000 + 0.5), 1, 'first');
        EntireSAPAmplitude{j} = m_amplitude(StartIndex:EndIndex);
        EntireSAPTime{j} = SAPFeature_Time(StartIndex:EndIndex) - Notes.offsets(j)/1000;
        EntireSAPEntropy{j} = m_Entropy(StartIndex:EndIndex);
    end
    AllAmplitudes = [AllAmplitudes Begin EntireAmplitude End];
    AllTimes = [AllTimes Begin EntireTime End];
    
    AllSAPAmplitudes = [AllSAPAmplitudes Begin EntireSAPAmplitude End];
    AllSAPEntropies = [AllSAPEntropies Begin EntireSAPEntropy End];
    AllSAPTimes = [AllSAPTimes Begin EntireSAPTime End];
    Labels = [Labels 'Q' Notes.labels 'q'];
end

for i = 1:length(Motif),
    figure;
    hold on;
    Matches = find(Labels == Motif(i));
    NextSylls = Labels(Matches+1);
    UniqueNextSylls = unique(NextSylls);
    NextSyllIndices = Matches + 1;
    for j = 1:length(UniqueNextSylls),
        subplot(length(UniqueNextSylls), 1, j);
        hold on;
        disp(['Syll ', Motif(i), ': Next Syll ', UniqueNextSylls(j), ': Colour ', Colours(j)]);
        TempMatches = find(Labels(NextSyllIndices) == UniqueNextSylls(j));
        TempMatches = NextSyllIndices(TempMatches) - 1;
        for k = 1:length(TempMatches),
            plot(AllSAPTimes{TempMatches(k)}, AllSAPAmplitudes{TempMatches(k)}, Colours(j));
        end
        axis tight;
    end
end

for i = 1:length(Motif),
    figure;
    hold on;
    Matches = find(Labels == Motif(i));
    NextSylls = Labels(Matches+1);
    UniqueNextSylls = unique(NextSylls);
    NextSyllIndices = Matches + 1;
    for j = 1:length(UniqueNextSylls),
        subplot(length(UniqueNextSylls), 1, j);
        hold on;
        disp(['Syll ', Motif(i), ': Next Syll ', UniqueNextSylls(j), ': Colour ', Colours(j)]);
        TempMatches = find(Labels(NextSyllIndices) == UniqueNextSylls(j));
        TempMatches = NextSyllIndices(TempMatches) - 1;
        for k = 1:length(TempMatches),
            plot(AllSAPTimes{TempMatches(k)}, AllSAPEntropies{TempMatches(k)}, Colours(j));
        end
        axis tight;
    end
end

for i = 1:length(Motif),
    figure;
    hold on;
    Matches = find(Labels == Motif(i));
    NextSylls = Labels(Matches+1);
    UniqueNextSylls = unique(NextSylls);
    NextSyllIndices = Matches + 1;
    
    TempAmplitudes = [];
    
    AmplitudeLens = cellfun(@length, AllSAPAmplitudes);
    MinLen = min(AmplitudeLens(find(AmplitudeLens > 1)));
    
    for j = 1:length(UniqueNextSylls),
        if (i ~= length(Motif))
            if ((UniqueNextSylls(j) == '?') || (UniqueNextSylls(j) == Motif(i+1)))
                hold on;
                disp(['Syll ', Motif(i), ': Next Syll ', UniqueNextSylls(j), ': Colour ', Colours(j)]);
                TempMatches = find(Labels(NextSyllIndices) == UniqueNextSylls(j));
                TempMatches = NextSyllIndices(TempMatches) - 1;
                for k = 1:length(TempMatches),
                    TempAmplitudes(k,:) = AllSAPAmplitudes{TempMatches(k)}(1:MinLen);
                    
%                    plot(AllTimes{TempMatches(k)}, AllAmplitudes{TempMatches(k)}, Colours(j));
                end
                errorbar(AllSAPTimes{TempMatches(1)}(1:MinLen), mean(TempAmplitudes), std(TempAmplitudes), Colours(j));
                axis tight;
            end
        else
            if ((UniqueNextSylls(j) == '?') || (UniqueNextSylls(j) == Motif(1)))
                hold on;
                disp(['Syll ', Motif(i), ': Next Syll ', UniqueNextSylls(j), ': Colour ', Colours(j)]);
                TempMatches = find(Labels(NextSyllIndices) == UniqueNextSylls(j));
                TempMatches = NextSyllIndices(TempMatches) - 1;
                for k = 1:length(TempMatches),
                    TempAmplitudes(k,:) = AllSAPAmplitudes{TempMatches(k)}(1:MinLen);
                end
                errorbar(AllSAPTimes{TempMatches(1)}(1:MinLen), mean(TempAmplitudes), std(TempAmplitudes), Colours(j));
                axis tight;
            end
        end
    end
end