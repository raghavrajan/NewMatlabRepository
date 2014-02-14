function [] = BatchSAPFeatureCalculation(BatchFileName,Motif)

fid = fopen(BatchFileName,'r');

while (~(feof(fid)))
    tline = fgetl(fid);
    NoteFile = [tline,'.not.mat'];
    SoundFile = tline;
    
    if ~(exist(NoteFile,'file'))
        continue;
    end
    
    Notes = load(NoteFile);
    
    for i = 1:length(Motif),
        MotifNotes = find(Notes.labels == Motif(i));
        if (length(MotifNotes) > 0)
            CalculateSAPFeatures(SoundFile);
            disp(['Calculated SAP Features for ',tline]);
            break;
        end
    end
end

fclose(fid);
