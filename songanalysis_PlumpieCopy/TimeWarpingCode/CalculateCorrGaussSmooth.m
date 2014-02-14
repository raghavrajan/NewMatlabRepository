function [Correlation] = CalculateCorrGaussSmooth(SpikeTrain, MedianMotif, Latency, Width)

Time = 0:0.0001:(MedianMotif.Length + 0.1 + Latency);

FR = zeros(length(SpikeTrain),length(Time));

for i = 1:length(SpikeTrain),
    Indices = round(([SpikeTrain{i}] + 0.1)/0.0001);
    Indices(find(Indices == 0)) = 1;
    FR(i,Indices) = 1;   
end

XGauss = 1:1:(1 + round(2 * 4 * Width * (1/0.0001)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = 1/((Width/0.0001) * sqrt(2 * pi)) * exp(-(XGauss.*XGauss)/(2 * (Width/0.0001) * (Width/0.0001)));

% The 4 in the previous line refers to the fact that the gaussian window is
% terminated at 4 times the standard deviation

Correlation = 0;

for i = 1:size(FR,1),
    ST(i,:) = conv(FR(i,:),GaussWin);
end

for i = 1:size(FR,1),
    for j = (i+1):size(FR,1),
        ST1 = ST(i,:);
        ST2 = ST(j,:);
        Correlation = Correlation + ((ST1 - mean(ST1)) * (ST2 - mean(ST2))')/(norm(ST1 - mean(ST1)) * norm(ST2 - mean(ST2)));
        %Correlation = Correlation + (sum((ST1 - mean(ST1)).*(ST2 - mean(ST2))))/(sqrt(sum((ST1 - mean(ST1)).^2) * sum((ST2 - mean(ST2)).^2)));
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