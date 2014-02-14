function [] = FindMatchingSyllables(Bout, Threshold, OutputFileName)

for i = 1:length(Bout.MaxBoutSeqMatch),
    Matches = Bout.MaxBoutSeqMatch{i} > Threshold;
    h=[1 -1];
    temp = zeros(size(Matches,1), size(Matches,2));
    temp(find(Matches > 0)) = 1;
    trans=conv(h,temp);
    onsets=find(trans > 0);
    offsets=find(trans < 0);
    for j = 1:length(onsets),
    end
end
