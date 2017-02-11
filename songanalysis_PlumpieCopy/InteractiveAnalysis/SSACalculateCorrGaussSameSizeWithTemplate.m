function [Corr] = SSACalculateCorrGaussSameSizeWithTemplate(SpikeTrain, MedianMotif, Width, PreSongStartDuration, PreSongEndDuration, ContextString, GaussianLen, Template)

Fs = 10000;

Time = -PreSongStartDuration:1/Fs:(MedianMotif.Length -PreSongEndDuration);

FR = zeros(length(SpikeTrain),length(Time));

for i = 1:length(SpikeTrain),
    Indices = round(([SpikeTrain{i}(find((SpikeTrain{i} >= Time(1)) & (SpikeTrain{i} < Time(end))))] + PreSongStartDuration) * Fs);
    Indices(find(Indices == 0)) = 1;
    FR(i,Indices) = 1;   
end

FR1 = zeros(length(Template),length(Time));

for i = 1:length(Template),
    Indices = round(([Template{i}(find((Template{i} >= Time(1)) & (Template{i} < Time(end))))] + PreSongStartDuration) * Fs);
    Indices(find(Indices == 0)) = 1;
    FR1(i,Indices) = 1;   
end

XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));

% The 4 in the previous line refers to the fact that the gaussian window is
% terminated at 4 times the standard deviation

Correlation = 0;
Correlation2 = 0;

for i = 1:size(FR,1),
    ST(i,:) = conv(FR(i,:),GaussWin, 'same');
end

for i = 1:size(FR1,1),
    ST1(i,:) = conv(FR1(i,:),GaussWin, 'same');
end

Index = 0;
for i = 1:(size(ST1,1)),
    ST_A = ST1(i,:);
    for j = 1:size(ST,1),
        Index = Index + 1;
        ST_B = ST(j,:);
        if ((norm(ST_A - mean(ST_A)) * norm(ST_B - mean(ST_B))) == 0)
            Temp(Index) = 0;
        else
            Temp(Index) = ((ST_A - mean(ST_A)) * (ST_B - mean(ST_B))')/(norm(ST_A - mean(ST_A)) * norm(ST_B - mean(ST_B)));
        end
    end
end

if (size(FR,1) == 1)
    Temp = [0.01];
end
Indices = find(Temp);

Corr(1,1) = mean(Temp(Indices));
Corr(1,2) = std(Temp(Indices));

% figure;
% Edges = -1:0.01:1;
% PST = histc(Temp(find(Temp ~= 0)), Edges);
% plot(Edges, PST);

% figure;
% set(gcf,'Color','w');
% 
% plot(Time,ST(1,(NoOfPoints:end)),'b');
% hold on;
% plot(Time,ST(2,(NoOfPoints:end)),'r');
% plot(Time,ST(3,(NoOfPoints:end)),'k');

disp(['Correlation with same size: ', ContextString, ': Gaussian Width = ', num2str(Width * 1000),' ms and correlation = ',num2str(Corr(1)), ' +/- ', num2str(Corr(2))]);
Corr = [Width Corr];