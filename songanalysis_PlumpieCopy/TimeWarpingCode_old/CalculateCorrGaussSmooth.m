function [Correlation] = CalculateCorrGaussSmooth(SpikeTrain, MedianMotif, Latency, Width)

Time = 0:0.0001:(MedianMotif.Length + Latency);

FR = zeros(length(SpikeTrain),length(Time));

for i = 1:length(SpikeTrain),
    Indices = round([SpikeTrain{i}]/0.0001);
    Indices(find(Indices == 0)) = 1;
    FR(i,Indices) = 1;   
end

GaussianWidth = Width/0.0001;

NoOfPoints = 1 + round(2 * 4 * GaussianWidth); 
% The 4 in the previous line refers to the fact that the gaussian window is
% terminated at 4 times the standard deviation

GaussWin = gausswin(NoOfPoints, 7.5);

Correlation = 0;

for i = 1:size(FR,1),
    ST(i,:) = conv(FR(i,:),GaussWin);
end

for i = 1:size(FR,1),
    for j = (i+1):size(FR,1),
        ST1 = ST(i,:);
        ST2 = ST(j,:);
        Correlation = Correlation + ((ST1 - mean(ST1)) * (ST2 - mean(ST2))')/(norm(ST1 - mean(ST1)) * norm(ST2 - mean(ST2)));
    end
end

Correlation = (Correlation * 2)/(size(FR,1) * (size(FR,1) - 1));

% figure;
% set(gcf,'Color','w');
% 
% plot(Time,ST(1,(NoOfPoints:end)),'b');
% hold on;
% plot(Time,ST(2,(NoOfPoints:end)),'r');
% plot(Time,ST(3,(NoOfPoints:end)),'k');

disp(['Gaussian Width = ', num2str(Width * 1000),' ms and Correlation = ',num2str(Correlation)]);