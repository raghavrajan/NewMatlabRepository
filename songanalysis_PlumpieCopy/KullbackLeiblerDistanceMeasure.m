function [KLDist] = KullbackLeiblerDistanceMeasure(Data1, Data2, NumBinsX, NumBinsY)

Limits = [0 max([max(max(Data1)) max(max(Data2))])];

EdgesX = linspace(Limits(1), Limits(2), NumBinsX);
EdgesY = linspace(Limits(1), Limits(2), NumBinsY);

for i = 2:length(EdgesX),
    Indices = find((Data1(:,1) >= EdgesX(i-1)) & (Data1(:,1) < EdgesX(i)));
    Dist1(i-1,:) = histc(Data1(Indices,2), EdgesY);
    
    Indices = find((Data2(:,1) >= EdgesX(i-1)) & (Data2(:,1) < EdgesX(i)));
    Dist2(i-1,:) = histc(Data2(Indices,2), EdgesY);
end

Dist1 = Dist1/size(Data1,1);
Dist2 = Dist2/size(Data2,1);

KLDist = 0;
for i = 1:size(Dist1,1),
    for j = 1:size(Dist1,2),
        if ((Dist1(i,j) > 0) && (Dist2(i,j) > 0))
            KLDist = KLDist + Dist1(i,j)*log2(Dist1(i,j)/Dist2(i,j));
        end
    end
end


disp('Finished KL analysis');