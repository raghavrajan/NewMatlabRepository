function [entropy] = CalculateEntropy(rawsong,Fs);

entropy = [];
index = 1;
entropy_index = 1;
flag = 1;
window_size = round(0.008 * Fs);
window_step = window_size - round(0.9*window_size);

while (flag)
if ((index + window_size) > length(rawsong))
flag = 0;
break;
end
temp = rawsong(index:(index + window_size));
tempfft = fft(temp);
Pyy = temp.*conj(temp)/(window_size);
entropy(entropy_index,1) = geomean(Pyy(1:(window_size + 1)))/mean(Pyy(1:(window_size + 1)));
entropy(entropy_index,2) = index/Fs;
entropy_index = entropy_index + 1;
index = index + window_step;
end