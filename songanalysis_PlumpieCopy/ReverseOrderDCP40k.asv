% this script reads a 40kHz dcp song and segment files, then make a reverse ordered song.
function ReverseOrderDCP40k(songname, segonsetfile, segendfile);

% read a song file 
song = readDCP40k(songname);
songlen = length(song);

%read segmentation files
%(unit is # samples from the song onset, not msec!)
segonset = readtextfile(segonsetfile);
segonset = str2num(segonset);
segend = readtextfile(segendfile);
segend = str2num(segend);
nseg = length(segonset);

% make a reverse order song
roSong = song(1:segonset(1)-1);
for n=1:nseg-1
    roSong = [song(segonset(n):segend(n));roSong];
    roSong = [song(segend(n)+1:segonset(n+1)-1);roSong];
end
roSong = [song(segonset(nseg):segend(nseg));roSong];
roSong = [song(segend(nseg)+1:songlen);roSong];

% plot the roSong
figure
subplot(2,1,1)
plot(song)
title('original song')
subplot(2,1,2)
plot(roSong)
title('reverse ordered song')
writeDCP40k_new(songname,roSong,'ro-');
