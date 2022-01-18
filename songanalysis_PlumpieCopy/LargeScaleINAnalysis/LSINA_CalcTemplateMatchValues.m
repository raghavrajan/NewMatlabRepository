function [TemplateMatchValues] = LSINA_CalcTemplateMatchValues(BirdParameters)

StretchValues = [-3:3:3];
SyllableIndex = 0;

Templates = load(BirdParameters.TemplateFile);
for i = 1:length(Templates.SyllableTemplates),
    SyllablesWithTemplates(i) = Templates.SyllableTemplates{i}{1}.MotifTemplate(1).Label;
end

if (BirdParameters.Continuousdata == 0)
    for i = 1:length(BirdParameters.SongFileNames),
        if (isempty(BirdParameters.NoteInfo{i}.onsets))
            TemplateMatchValues{i} = [];
            continue;
        end
        [RawData, Fs] = ASSLGetRawData(BirdParameters.DataDirectory, BirdParameters.SongFileNames{i}, BirdParameters.FileType, 0);

        Time = (1:1:length(RawData))/Fs;

        WinSize = Templates.SyllableTemplates{1}{1}.MotifTemplate(1).FFTWinSize;
        WinOverlap = Templates.SyllableTemplates{1}{1}.MotifTemplate(1).FFTWinOverlap;

        PreDur = 0.5; % Pre duration before syllable onset in sec;
        PostDur = 0.5; % Post duration after syllable onset in sec;

        for j = 1:length(BirdParameters.NoteInfo{i}.onsets),
            SyllableOnsetTimeIndex = round(BirdParameters.NoteInfo{i}.onsets(j)/1000 * Fs);
            SyllableOnsetWindow = [(SyllableOnsetTimeIndex - round(PreDur * Fs)) (SyllableOnsetTimeIndex + round(PostDur * Fs))];
            if (SyllableOnsetWindow(1) <= 0)
                SyllableOnsetWindow(1) = 1;
            end
            if (SyllableOnsetWindow(2) > length(RawData))
                SyllableOnsetWindow(2) = length(RawData);
            end

            SyllableRawData = RawData(SyllableOnsetWindow(1):SyllableOnsetWindow(2));
            [P, F, S1, T] = CalculateMultiTaperSpectrogram(SyllableRawData, Fs, WinSize, WinOverlap, 1.5);

            Freq1 = find((F >= 860) & (F <= 8600));
            S = log10(abs(S1(Freq1,:)));

            SyllableMatchFs = 1/(T(2) - T(1));

            clear TempMatch SyllableMatch MatchLen;

            TemplateMatchIndex = find(SyllablesWithTemplates == BirdParameters.NoteInfo{i}.labels(j));

            if (~isempty(TemplateMatchIndex))
                TempSyllTemplate = Templates.SyllableTemplates{TemplateMatchIndex}{min(length(Templates.SyllableTemplates{TemplateMatchIndex}),3)};
                TempMatch = ASSLTemplateMatch(S, TempSyllTemplate.MotifTemplate, StretchValues);

                ActualPreDur = 0.015; % Pre dur for locating match
                if ((SyllableOnsetTimeIndex - round(PreDur * Fs)) > 0)
                    SyllableOnsetWindow = [(round((PreDur - ActualPreDur) * SyllableMatchFs)) (round((PreDur + ActualPreDur) * SyllableMatchFs))];
                else
                    ExtraTime = (SyllableOnsetTimeIndex - round(PreDur * Fs))/Fs;
                    SyllableOnsetWindow = [(round((PreDur - ExtraTime - ActualPreDur) * SyllableMatchFs)) (round((PreDur - ExtraTime + ActualPreDur) * SyllableMatchFs))];
                end

                if (SyllableOnsetWindow(1) <= 0)
                    SyllableOnsetWindow(1) = 1;
                end

                if ((SyllableOnsetWindow(1) > size(TempMatch,2)) || (SyllableOnsetWindow(2) > size(TempMatch,2)))
                    TemplateMatchValues{i}(j) = NaN;
                else
                    [MaxVal, MaxIndex] = max(TempMatch(SyllableOnsetWindow(1):SyllableOnsetWindow(2)));
                    TemplateMatchValues{i}(j) = MaxVal;
                end
            else
                TemplateMatchValues{i}(j) = NaN;
            end
        end
    end
end
