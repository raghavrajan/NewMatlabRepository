%reads DCP songfile (40kHz) into a vector and plot it.

function song = readDCP40k_plot(songname);

SAMP_RATE = 40000;

fid = fopen(songname,'r','b');
fread(fid,11,'int16');
song = fread(fid,'int32');
fclose(fid);


figure
subplot(2,1,1)
specgram(song,300,40000,280,200)
title(songname);
subplot(2,1,2)
t = [1:length(song)]/SAMP_RATE;
plot(t, song);
xlabel('Time (seconds)');