function [Data, PropNoMotifBouts, BoutLen, BoutPropMotifSylls] = MicrolesionAnalysesNoteFiles(NoteFileList, InterBoutInterval, Motif, IntroNotes, Calls, BoutLogFileFid, SeqLogFileFid)

Fid = fopen(NoteFileList, 'r');

Temp = textscan(Fid, '%s', 'delimiter', '\n');

fclose(Fid);

NoteFiles = Temp{1};

AllLabels = [];
BoutNo = 0;
for i = 1:min([25 length(NoteFiles)]),
   Notes = load([NoteFiles{i}, '.not.mat']); 
   StartIndex = 1;
   for j = 2:length(Notes.onsets)-1,
       if ((Notes.onsets(j) - Notes.offsets(j-1))/1000 > InterBoutInterval)
           BoutNo = BoutNo + 1;
           EndIndex = j-1;
           Data{BoutNo}.Labels = ['Q', Notes.labels(StartIndex:EndIndex), 'q'];
           Data{BoutNo}.Onsets = [Notes.onsets(StartIndex)/1000; Notes.onsets(StartIndex:EndIndex)/1000; Notes.offsets(EndIndex)/1000];
           Data{BoutNo}.Offsets = [Notes.onsets(StartIndex)/1000; Notes.offsets(StartIndex:EndIndex)/1000; Notes.offsets(EndIndex)/1000];
           
           NonLabels = find(Data{BoutNo}.Labels == '0');
           Data{BoutNo}.Labels(NonLabels) = [];
           Data{BoutNo}.Onsets(NonLabels) = [];
           Data{BoutNo}.Offsets(NonLabels) = [];
           
           AllLabels = [AllLabels Data{BoutNo}.Labels];
           StartIndex = j;
       end
   end
   BoutNo = BoutNo + 1;
   EndIndex = length(Notes.labels);
   Data{BoutNo}.Labels = ['Q', Notes.labels(StartIndex:EndIndex), 'q'];
   Data{BoutNo}.Onsets = [Notes.onsets(StartIndex)/1000; Notes.onsets(StartIndex:EndIndex)/1000; Notes.offsets(EndIndex)/1000];
   Data{BoutNo}.Offsets = [Notes.onsets(StartIndex)/1000; Notes.offsets(StartIndex:EndIndex)/1000; Notes.offsets(EndIndex)/1000];
   
   NonLabels = find(Data{BoutNo}.Labels == '0');
   Data{BoutNo}.Labels(NonLabels) = [];
   Data{BoutNo}.Onsets(NonLabels) = [];
   Data{BoutNo}.Offsets(NonLabels) = [];
   
   AllLabels = [AllLabels Data{BoutNo}.Labels];
end

NoMotifBouts = 0;
ValidBout = 0;
TransProb = zeros(length(Motif));
TotalSylls = zeros(length(Motif),1);
TotalBouts = 0;

for i = 1:length(Data),
    SequenceIndex = 0;
    GoodBout = 0;
    for j = 1:length(IntroNotes),
        if (~isempty(find(Data{i}.Labels == IntroNotes(j))))
            GoodBout = 1;
            break;
        end
    end
    if (GoodBout == 1)
        YesMotif = 0;
        TotalBouts = TotalBouts + 1;
        for j = 1:length(Motif),
            if (~isempty(find(Data{i}.Labels == Motif(j))))
                YesMotif = 1;
                break;
            end
        end
        
        if (YesMotif == 0)
            if (isempty(find(Data{i}.Labels == '?')))
                NoMotifBouts = NoMotifBouts + 1;
            end
        end
        
        ValidBout = ValidBout + 1;
        BoutStarting(ValidBout) = Data{i}.Labels(2);
        BoutLength(ValidBout) = Data{i}.Offsets(end) - Data{i}.Onsets(1);
        MotifSylls = 0;
        for j = 1:length(Motif),
            Sylls = find(Data{i}.Labels == Motif(j));
            if (~isempty(Sylls))
                TotalSylls(j) = TotalSylls(j) + length(Sylls);
                for k = 1:length(Motif),
                    TransProb(j,k) = TransProb(j,k) + length(find(Data{i}.Labels(Sylls + 1) == Motif(k)));
                end
                MotifSylls = MotifSylls + length(Sylls);
            end
        end
        PropMotifSylls(ValidBout) = MotifSylls/length(Data{i}.Labels);
        
        NonMotifSylls = find((Data{i}.Labels == 'q') | (Data{i}.Labels == 'Q'));

        for k = 1:length(IntroNotes),
            NonMotifSylls = [NonMotifSylls find(Data{i}.Labels == IntroNotes(k))];
        end    

        for k = 1:length(Calls),
            NonMotifSylls = [NonMotifSylls find(Data{i}.Labels == Calls(k))];
        end
        NonMotifSylls = sort(NonMotifSylls);

        for k = 1:length(NonMotifSylls) - 1,
            TempSequence = Data{i}.Labels(NonMotifSylls(k)+1:NonMotifSylls(k+1)-1);
            TempSequence(find((TempSequence == '?'))) = [];
            if (~isempty(TempSequence))
                SequenceIndex = SequenceIndex + 1;
                Data{i}.Sequences{SequenceIndex} = TempSequence;
                Data{i}.SequenceLengths(SequenceIndex) = Data{i}.Offsets(NonMotifSylls(k+1)-1) - Data{i}.Onsets(NonMotifSylls(k)+1);
                Data{i}.SeqNoofSylls(SequenceIndex) = length(TempSequence);
            end    
        end
        Data{i}.BoutNoofSequences = SequenceIndex;
    end
end

TransProb = TransProb./repmat(TotalSylls, 1, length(Motif));


% Outputs

PropNoMotifBouts = NoMotifBouts/TotalBouts;
BoutLen = [mean(BoutLength) std(BoutLength)];
BoutPropMotifSylls = [mean(PropMotifSylls) std(PropMotifSylls)];

disp(NoteFileList);
disp(['Proportion of bouts without motif sylls is ', num2str(NoMotifBouts/TotalBouts)]);
disp(['Mean proportion of motif sylls in a bout is ', num2str(mean(PropMotifSylls))]);
disp(['Mean bout length is ', num2str(mean(BoutLength)), ' and std is ', num2str(std(BoutLength))]);
disp(unique(AllLabels));
save([NoteFileList, '.notefileanalysis.mat']);

clear TempLabels TempLabelLens
for i = 1:length(Data),
    TempLabels{i} = Data{i}.Labels;
    TempLabelLens(i) = length(Data{i}.Labels);
end

fprintf(BoutLogFileFid, '%s\n\n', NoteFileList);
[SortedLabels, SortedIndices] = sort(TempLabelLens);
for i = 1:length(SortedIndices),
    fprintf(BoutLogFileFid, '%s\n', TempLabels{SortedIndices(i)});
end
fprintf(BoutLogFileFid, '\n\n===================================================================\n\n');


clear TempSequences TempSeqLens
SeqIndex = 0;
for i = 1:length(Data),
    if (isfield(Data{i}, 'Sequences'))
        for j = 1:length(Data{i}.Sequences),
            SeqIndex = SeqIndex + 1;
            TempSequences{SeqIndex} = Data{i}.Sequences{j};
            TempSeqLens(SeqIndex) = length(Data{i}.Sequences{j});
        end
    end
end

fprintf(SeqLogFileFid, '%s\n\n', NoteFileList);
[SortedLabels, SortedIndices] = sort(TempSeqLens);
for i = 1:length(SortedIndices),
    fprintf(SeqLogFileFid, '%s\n', TempSequences{SortedIndices(i)});
end
fprintf(SeqLogFileFid, '\n\n===================================================================\n\n');

disp('Finished');