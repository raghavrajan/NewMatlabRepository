function [Corr] = SSACalculateCorrGaussWithPreData(FileInfo, SpikeTrain, MedianMotif, Width, PreSongStartDuration, PreSongEndDuration, ContextString, GaussianLen)

Fs = 10000;

Time = -(PreSongStartDuration + 1):1/Fs:(MedianMotif.Length -PreSongEndDuration + 1);
ActualTimeIndex(1) = find(Time <= -PreSongStartDuration, 1, 'last');
ActualTimeIndex(2) = find(Time <= (MedianMotif.Length - PreSongEndDuration), 1, 'last');

FR = zeros(length(SpikeTrain),length(Time));

for i = 1:length(SpikeTrain),
    Indices = round(([SpikeTrain{i}] + PreSongStartDuration + 1) * Fs);
    Indices(find(Indices == 0)) = 1;
    FR(i,Indices) = 1;   
    Indices = find((FileInfo.SpikeData.Times{FileInfo.Syllables.Index(i)} >= (FileInfo.Syllables.Start(i,1) - PreSongStartDuration - 1)) & (FileInfo.SpikeData.Times{FileInfo.Syllables.Index(i)} < (FileInfo.Syllables.Start(i,1) - PreSongStartDuration)));
    Indices = round((FileInfo.SpikeData.Times{FileInfo.Syllables.Index(i)}(Indices) - FileInfo.Syllables.Start(i,1) + PreSongStartDuration + 1) * Fs);
    Indices(find(Indices == 0)) = 1;
    
    FR(i,Indices) = 1;
    Indices = find((FileInfo.SpikeData.Times{FileInfo.Syllables.Index(i)} >= (FileInfo.Syllables.End(i,end) - PreSongEndDuration)) & (FileInfo.SpikeData.Times{FileInfo.Syllables.Index(i)} < (FileInfo.Syllables.End(i,end) - PreSongEndDuration + 1)));
    Indices = round((FileInfo.SpikeData.Times{FileInfo.Syllables.Index(i)}(Indices) - FileInfo.Syllables.End(i,end) + PreSongStartDuration + 1 + MedianMotif.Length) * Fs);
    FR(i,Indices) = 1;
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

Index = 0;
for i = 1:(size(FR,1)),
    for j = (i+1):size(FR,1),
        ST1 = ST(i,ActualTimeIndex(1):ActualTimeIndex(2));
        ST2 = ST(j,ActualTimeIndex(1):ActualTimeIndex(2));
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

disp(['Correlation with pre data: ', ContextString, ': Gaussian Width = ', num2str(Width * 1000),' ms and correlation = ',num2str(Corr(1)), ' +/- ', num2str(Corr(2))]);
Corr = [Width Corr];