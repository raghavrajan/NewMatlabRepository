function [] = JitterPoints(XData, YData, DirectionofJitter, MaxJitter)

NumPartitions = 10;
NeighbourRange = (max(YData) - min(YData))/NumPartitions;

Data = [XData(:) YData(:)];
Data = sortrows(Data, 2);
for i = 1:size(Data,1),
    NumNeighbours = length(find((Data(:,2) > (Data(i,2) - NeighbourRange/2)) & (Data(:,2) < (Data(i,2) + NeighbourRange/2))));
    if (NumNeighbours == 0)
        Data(i,1) = Data(i,1);
    else
        if (NumNeighbours <= 5)
            Data(i,1) = Data(i,1) + (rand*2*MaxJitter/3 - MaxJitter/3);
        else
            if (NumNeighbours <= 10)
                Data(i,1) = Data(i,1) + (rand*2*2*MaxJitter/3 - 2*MaxJitter/3);
            else
                Data(i,1) = Data(i,1) + (rand*2*MaxJitter - MaxJitter);
            end
        end
    end
end
plot(Data(:,1), Data(:,2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');