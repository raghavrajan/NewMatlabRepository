% this script reads a 32kHz dcp song and segment files, then make a reverse ordered song.
function ReverseOrderDCP_new(songname, segmentfile, samplerate);

% read a song file and segmentation files
if samplerate == 40000
    song = readDCP40k(songname);
end
if samplerate == 32000
    song = readDCP32k(songname);
end
songlen = length(song);

seg = str2num(readtextfile(segmentfile));
nseg = length(seg);
seg2(1) = seg(1);
t = seg(1);
for n=2:nseg
    t = t+seg(n);
    seg2(n) = t;
end

% make a reverse order song
roSong = song(1:round(seg2(1)*samplerate/1000)-1);
if round(seg2(nseg)*samplerate/1000)-1 <= length(song)
    for n=1:nseg-1
        roSong = [song(round(seg2(n)*samplerate/1000):round(seg2(n+1)*samplerate/1000)-1);roSong];
    end
end
if round(seg2(nseg)*samplerate/1000)-1 > length(song)
    for n=1:nseg-2
        roSong = [song(round(seg2(n)*samplerate/1000):round(seg2(n+1)*samplerate/1000)-1);roSong];
    end
    roSong = [song(round(seg2(nseg-1)*samplerate/1000):length(song));roSong];
end
% plot the roSong
figure
subplot(2,1,1)
plot(song)
title('original song')
subplot(2,1,2)
plot(roSong)
title('reverse ordered song')
if samplerate == 40000
    writeDCP40k_new(songname,roSong,'ro-');
end
if samplerate == 32000
    writeDCP32k_new(songname,roSong,'ro-');
end


