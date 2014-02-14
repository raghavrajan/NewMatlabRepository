function [] = SongSyllableDurationAnalysis(FileName,FileExt,MotifFileName,motif)


figure;
MotifDurations = axes('position',[0.1 0.675 0.8 0.3]);
hold on;
IntroNoteNos = axes('position',[0.1 0.35 0.8 0.3]);
hold on;
SyllableDurations = axes('position',[0.1 0.025 0.8 0.3]);
hold on;

fid = fopen(MotifFileName,'r');
motif_count = 0;
while ~(feof(fid))
    tline = fgetl(fid);
    comma_index = strfind(tline,',');
    hyphen_index = strfind(tline,'-');
    if (length(comma_index) > 0)
        first_no = str2num(tline(1:(comma_index - 1)));
        second_no = str2num(tline((comma_index + 1):(hyphen_index - 1)));
        if ((first_no) < 10)
            first_no = ['00',num2str(first_no)];
        else
            if ((first_no) < 100)
                first_no = ['0',num2str(first_no)];
            else
                first_no = [num2str(first_no)];
            end
        end

        if ((second_no) < 10)
            second_no = ['00',num2str(second_no)];
        else
            if ((second_no) < 100)
                second_no = ['0',num2str(second_no)];
            else
                second_no = [num2str(second_no)];
            end
        end        
        NoteFileName = strcat(FileName,'.',first_no,'.',second_no,FileExt,'.not.mat');
    else
        first_no = str2num(tline(1:(hyphen_index - 1)));

        if ((first_no) < 10)
            first_no = ['00',num2str(first_no)];
        else
            if ((first_no) < 100)
                first_no = ['0',num2str(first_no)];
            else
                first_no = [num2str(first_no)];
            end
        end

        NoteFileName = strcat(FileName,'.',first_no,FileExt,'.not.mat');
    end
    
    Notes = load(NoteFileName);
        
    IntroNoteIndices = find(Notes.labels == 'i');
    SyllableIndices = setdiff((1:1:length(Notes.onsets)),IntroNoteIndices);
    
    Syllables.labels = Notes.labels(SyllableIndices);
    Syllables.onsets = Notes.onsets(SyllableIndices);
    Syllables.offsets = Notes.offsets(SyllableIndices);
    
    ActualSongOnsetIndices = find(Notes.labels == 'a');
    
    SongOnsetsIndices = (find(Syllables.labels == 'a'));
    colour_string = ['krgbcmy'];
    
    for Songs = 1:length(SongOnsetsIndices),
        if (Songs == length(SongOnsetsIndices))
            SongLength = Syllables.offsets(end) - Syllables.onsets((SongOnsetsIndices(Songs)));
            if (Songs > 1)
                NoofIntroNotes = length(find(Notes.labels((ActualSongOnsetIndices(Songs-1)):(ActualSongOnsetIndices(Songs))) == 'i'));
            else
                NoofIntroNotes = length(find(Notes.labels == 'i'));
            end
            NoofSyllables = length(Syllables.onsets(SongOnsetsIndices(Songs):end));
            for SyllableNo = SongOnsetsIndices(Songs):length(Syllables.onsets),
                axes(SyllableDurations);
                SyllableLength = Syllables.offsets(SyllableNo) - Syllables.onsets(SyllableNo);
                plot(motif_count,SyllableLength,[colour_string(SyllableNo - SongOnsetsIndices(Songs) + 1),'s']);
            end
        else
            SongLength = Syllables.offsets((SongOnsetsIndices(Songs + 1) - 1)) - Syllables.onsets((SongOnsetsIndices(Songs)));
            NoofSyllables = length(Syllables.onsets((SongOnsetsIndices(Songs)):(SongOnsetsIndices(Songs + 1) - 1)));
            for SyllableNo = SongOnsetsIndices(Songs):(SongOnsetsIndices(Songs + 1) - 1),
                axes(SyllableDurations);
                SyllableLength = Syllables.offsets(SyllableNo) - Syllables.onsets(SyllableNo);
                plot(motif_count,SyllableLength,[colour_string(SyllableNo - (SongOnsetsIndices(Songs)) + 1),'s']);
            end
            if (Songs == 1);
                NoofIntroNotes = length(find(Notes.labels(1:ActualSongOnsetIndices(Songs)) == 'i'));
            else
                NoofIntroNotes = length(find(Notes.labels((ActualSongOnsetIndices(Songs-1)):(ActualSongOnsetIndices(Songs))) == 'i'));
            end
        end
        axes(MotifDurations);
        plot(motif_count,SongLength,'ks');
        axes(IntroNoteNos);
        plot(motif_count,NoofIntroNotes,'ks');
        motif_count = motif_count + 1;
    end

    disp(['Loaded Notes file ',NoteFileName]);
    
end
fclose(fid);
