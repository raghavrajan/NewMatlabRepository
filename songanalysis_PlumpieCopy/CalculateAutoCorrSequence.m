function [R] = CalculateAutoCorrSequence(Data, WindowSize, Lags)

for i = 1:length(Data)-WindowSize,
    for j = 1:length(Lags),
        R(i,j) = 0;
        for k = i:i+WindowSize,
            R(i,j) = R(i,j) + Data(k)*Data(k+Lags(j));
        end
    end
end
