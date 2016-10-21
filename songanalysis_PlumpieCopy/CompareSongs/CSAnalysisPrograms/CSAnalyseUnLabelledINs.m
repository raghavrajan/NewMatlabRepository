function [] = CSAnalyseUnLabelledINs(CSData)

figure;
Colors = 'rgbcmk';
Symbols = 'o+d';
for i = 1:CSData.NoofDays,
   BoutStarts = find(CSData.AllLabels{i} == 'Q');
   BoutEnds = find(CSData.AllLabels{i} == 'q');
   for j = 1:length(BoutStarts),
       MotifSylls = regexp(CSData.AllLabels{i}(BoutStarts(j):BoutEnds(j)), ['[', CSData.MotifSyllLabels, ']']);
       if (isempty(MotifSylls))
           continue;
       end
       Durations = CSData.AllOffsets{i}(BoutStarts(j)+1:BoutEnds(j)-2) - CSData.AllOnsets{i}(BoutStarts(j)+1:BoutEnds(j)-2);
       Intervals = CSData.AllOnsets{i}(BoutStarts(j)+2:BoutEnds(j)-1) - CSData.AllOffsets{i}(BoutStarts(j)+1:BoutEnds(j)-2);
       plot(Intervals, [Colors(mod(i-1, length(Colors)) + 1), Symbols(ceil(i/length(Colors))), '-']);
       hold on;
   end
end

figure;
for i = 1:CSData.NoofDays,
   BoutStarts = find(CSData.AllLabels{i} == 'Q');
   BoutEnds = find(CSData.AllLabels{i} == 'q');
   for j = 1:length(BoutStarts),
       MotifSylls = regexp(CSData.AllLabels{i}(BoutStarts(j):BoutEnds(j)), ['[', CSData.MotifSyllLabels, ']']);
       if (isempty(MotifSylls))
           continue;
       end
       CSData.NormAllFeats{i} = (CSData.AllFeats{i} - repmat(mean(CSData.Data{i}.FeatValues), size(CSData.AllFeats{i},1), 1))./repmat(std(CSData.Data{i}.FeatValues), size(CSData.AllFeats{i},1), 1);
       [Coeff{i}, Score{i}, Latent{i}] = princomp(CSData.NormAllFeats{i});
       for k = (BoutStarts(j) + 1):(BoutEnds(j) - 2),
           Distances{i}{j}(k - BoutStarts(j)) = pdist(Score{i}(k:k+1,:));
       end
       plot(Distances{i}{j}, [Colors(mod(i-1, length(Colors)) + 1), Symbols(ceil(i/length(Colors))), '-']);
       hold on;
   end
end

disp('Finished analysing data');