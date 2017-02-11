%reads DCP songfile (40kHz) into a vector

function song = readDCP40k(songname);

SAMP_RATE = 40000;

fid = fopen(songname,'r','b');
fread(fid,11,'int16');
song = fread(fid,'int32');
fclose(fid);

%t = [1:length(song)]/SAMP_RATE;
%figure
%plot(t, song);
%title(songname);
%xlabel('Time (seconds)');