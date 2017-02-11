% This script read a DCP 40k song file and make the reverse song.
function reverseDCP40k;

songname = input('Enter name of dcp song file: ', 's');
song = readDCP40k(songname);
revsong=flipud(song);

figure
subplot(2,1,1)
plot(song)
title('original song')
subplot(2,1,2)
plot(revsong)
title('reversed song')
writeDCP40k_new(songname,revsong,'rev');
