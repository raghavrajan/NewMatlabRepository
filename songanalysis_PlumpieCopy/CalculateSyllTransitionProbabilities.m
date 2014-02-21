function [SyllTransitionProb, AllLabels] = CalculateSyllTransitionProbabilities(NoteFileList, NoteFileDir, InterBoutInterval)

Fid = fopen(NoteFileList, 'r');
Files = textscan(Fid, '%s', 'DeLimiter', '\n');
Files = Files{1};

FileSep = filesep;
if (NoteFileDir(end) ~= FileSep)
    NoteFileDir(end+1) = FileSep;
end

AllLabels = [];
for i = 1:length(Files),
    Notes{i} = load([NoteFileDir, Files{i}, '.not.mat']);
    
    InterSyllIntervals = (Notes{i}.onsets(2:end) - Notes{i}.offsets(1:end-1))/1000;
    
    LongInts = find(InterSyllIntervals >= InterBoutInterval);
    
    if (isempty(LongInts))
        Bouts = [1 length(Notes{i}.onsets)];
    else
        Bouts = [1 LongInts(1)];
        for j = 1:length(LongInts)-1,
            Bouts = [Bouts; [(LongInts(j)+1) LongInts(j+1)]];
        end
        Bouts = [Bouts; [(LongInts(end)+1) length(Notes{i}.onsets)]];
    end
    
    for j = 1:size(Bouts,1),
        AllLabels = [AllLabels 'Q' Notes{i}.labels(Bouts(j,1):Bouts(j,2)) 'q'];
    end
end

UniqueLabels = unique(AllLabels);

for i = 1:length(UniqueLabels),
    SyllTransitionProb{i}.Syll = UniqueLabels(i);
    if (UniqueLabels(i) ~= 'q')
        Matches = find(AllLabels == UniqueLabels(i));
        for j = 1:length(UniqueLabels),
            SyllTransitionProb{i}.NextSyll(j) = UniqueLabels(j);
            SyllTransitionProb{i}.Prob(j) = length(find(AllLabels(Matches + 1) == UniqueLabels(j)))/length(Matches);
        end
    else
        SyllTransitionProb{i}.Prob = zeros(1, length(UniqueLabels));
    end
end

fprintf('Syllable transition probabilities \n');
for i = 1:length(SyllTransitionProb),
    NonZeroTrans = find(SyllTransitionProb{i}.Prob > 0);
    for j = 1:length(NonZeroTrans),
        fprintf('%s(%g), ', [SyllTransitionProb{i}.Syll, SyllTransitionProb{i}.NextSyll(NonZeroTrans(j))], SyllTransitionProb{i}.Prob(NonZeroTrans(j)));
    end
    fprintf('\n');
end

fprintf('Syllable frequencies\n');
for i = 1:length(SyllTransitionProb),
    fprintf('%c-%i, ', SyllTransitionProb{i}.Syll, length(find(AllLabels == SyllTransitionProb{i}.Syll)));
end
fprintf('\n');
disp('Finished');