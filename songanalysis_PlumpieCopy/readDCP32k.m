%reads DCP songfile (32kHZ) into a vector

function song = readDCP32k(songname);

SAMP_RATE = 32000;

fid = fopen(songname,'r','b');
fread(fid,11,'int16');
song = fread(fid,'int16');
fclose(fid);

%t = [1:length(song)]/SAMP_RATE;
%figure
%plot(t, song);
%title(songname);
%xlabel('Time (seconds)');