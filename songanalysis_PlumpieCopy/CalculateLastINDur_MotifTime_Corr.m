function [] = CalculateLastINDur_MotifTime_Corr(INR)

for i = 1:length(INR.INs),
    LastINDur = (INR.BoutDetails(i).offsets(INR.INs{i}(end)) - INR.BoutDetails(i).onsets(INR.INs{i}(end)));
    LastINToMotifOnset = (INR.BoutDetails(i).onsets(INR.INs{i}(end) + 1) - INR.BoutDetails(i).onsets(INR.INs{i}(end)));
    INDur_TimeToStartMotif(i,:) = [LastINDur LastINToMotifOnset];
end

for i = 1:length(INR.WithinBoutINs),
    if (~isempty(INR.WithinBoutINs{i}))
        LastINDur = (INR.BoutDetails(INR.WithinBoutINBoutIndices(i)).offsets(INR.WithinBoutINs{i}(end)) - INR.BoutDetails(INR.WithinBoutINBoutIndices(i)).onsets(INR.WithinBoutINs{i}(end)));
        LastINToMotifOnset = (INR.BoutDetails(INR.WithinBoutINBoutIndices(i)).onsets(INR.WithinBoutINs{i}(end) + 1) - INR.BoutDetails(INR.WithinBoutINBoutIndices(i)).onsets(INR.WithinBoutINs{i}(end)));
        INDur_TimeToStartMotif(end+1,:) = [LastINDur LastINToMotifOnset];
    end
end

figure;
plot(INDur_TimeToStartMotif(:,1), INDur_TimeToStartMotif(:,2), 'ko');
[r, p] = corrcoef(INDur_TimeToStartMotif(:,1), INDur_TimeToStartMotif(:,2))
 