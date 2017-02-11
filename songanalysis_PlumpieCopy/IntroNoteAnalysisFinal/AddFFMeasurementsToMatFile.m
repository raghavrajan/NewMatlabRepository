function [INR] = AddFFMeasurementsToMatFile(INR, FF_File)

FFData = load(FF_File);

for i = 1:length(INR.BoutDetails),
    FF = zeros(size(INR.BoutDetails(i).Feats,1), 1);
    IN_Indices = find(INR.BoutDetails(i).labels == 'i');
    
    SongFileFlags = zeros(size(FFData.ffreq,1),1);
    for j = 1:size(FFData.ffreq,1),
        SongFileFlags(j) = length(strfind(FFData.ffreq{j,1}, INR.BoutDetails(i).SongFile));
    end
    SongFileIndices = find(SongFileFlags);
    for j = 1:length(IN_Indices),
        for k = 1:length(SongFileIndices),
            if (FFData.ffreq{SongFileIndices(k),16} == (IN_Indices(j) + INR.BoutDetails(i).BoutIndices(1) - 1))
                FF(IN_Indices(j)) = FFData.ffreq{SongFileIndices(k),2};
            end
        end
    end
    INR.BoutDetails(i).Feats = [INR.BoutDetails(i).Feats FF];
end

disp('Finished adding ff measurements to INR results file');
    