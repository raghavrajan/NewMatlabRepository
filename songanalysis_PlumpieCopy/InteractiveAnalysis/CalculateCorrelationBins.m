function [Corr, Temp] = CalculateCorrelationBins(SpikeTrain, Range, FullRange, Width)

Edges = 0:Width:(FullRange(2) - FullRange(1));

Correlation = 0;
Correlation2 = 0;

for i = 1:length(SpikeTrain),
    Spikes = SpikeTrain{i}(find(([SpikeTrain{i}] >= Range(1)) & ([SpikeTrain{i}] <= Range(2))));
    Spikes = Spikes - Range(1);
    ST(i,:) = histc(Spikes, Edges);
    FR(i,:) = histc(Spikes, Edges);
end

Index = 0;
for i = 1:(size(FR,1)),
    for j = (i+1):size(FR,1),
        ST1 = ST(i,:);
        ST2 = ST(j,:);
        if ((norm(ST1 - mean(ST1)) * norm(ST2 - mean(ST2))) == 0)
            Temp(i,j) = 0;
        else
            Correlation2 = Correlation2 + ((ST1 - mean(ST1)) * (ST2 - mean(ST2))')/(norm(ST1 - mean(ST1)) * norm(ST2 - mean(ST2)));
            Temp(i,j) = ((ST1 - mean(ST1)) * (ST2 - mean(ST2))')/(norm(ST1 - mean(ST1)) * norm(ST2 - mean(ST2)));
            Index = Index + 1;
            Correlation = Correlation + (sum((ST1 - mean(ST1)).*(ST2 - mean(ST2))))/(sqrt(sum((ST1 - mean(ST1)).^2) * sum((ST2 - mean(ST2)).^2)));
        end
    end
end

Correlation = (Correlation * 2)/(size(FR,1) * (size(FR,1) - 1));
Correlation2 = (Correlation2 * 2)/(size(FR,1) * (size(FR,1) - 1));
Indices = find(Temp);

Corr(1,1) = mean(Temp(Indices));
Corr(1,2) = std(Temp(Indices));

% figure;
% set(gcf,'Color','w');
% 
% plot(Time,ST(1,(NoOfPoints:end)),'b');
% hold on;
% plot(Time,ST(2,(NoOfPoints:end)),'r');
% plot(Time,ST(3,(NoOfPoints:end)),'k');

disp(['Gaussian Width = ', num2str(Width),' ms and correlation = ',num2str(Corr(1)), ' +/- ', num2str(Corr(2))]);
Corr = [Width Corr];