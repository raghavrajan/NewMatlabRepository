function [IFR] = CalculateIFR(SpikeTrain, Time)

SpikeTrain = sort(SpikeTrain);

Fs = 1/(Time(2) - Time(1));

IFR = zeros(size(Time));

for i = 2:length(SpikeTrain),
    IFR(ceil(SpikeTrain(i-1)*Fs)+1:ceil(SpikeTrain(i)*Fs)) = 1/(SpikeTrain(i) - SpikeTrain(i-1));
end


% disp('Finished calculating instantaneous firing rate');
