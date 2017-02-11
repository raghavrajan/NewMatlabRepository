% segment a song with a threshold.  
%this script is for songfiltnorm_wav2dcp_new.m

% Do not change!
% Do not change! 
% Do not change!
% Do not change!
% Do not change! 
% Do not change!
% Do not change!
% Do not change! 
% Do not change!
% Do not change!
% Do not change! 
% Do not change!
% Do not change!
% Do not change! 
% Do not change!

function [segonsettime3,segendtime3] = segment4songfiltnorm(songname,threshold,minseg,minamp,save);

%minseg is minmum duration of intersyllable intervals (msec)
minseg2 = 8; %minmum duration of syllables (msec)
%minamp is minimum amplitude of a sound segment, with which the segment will be considered as a syllable. A segment with smaller amplitude than minamp will be discarded.

%song = readDCP40k(songname);
song = songname;
songlen = length(song);

% rectify song signal
recsong = zeros(songlen,1);
      for i=1:songlen
          if song(i) < 0
              recsong(i) = song(i)*(-1);
          elseif song(i) >= 0
              recsong(i) = song(i);
          end
      end

%segment song with threshold
segonset = zeros(songlen,1);
segend = zeros(songlen,1);
      for i=1:songlen-1
       if recsong(i)<threshold & recsong(i+1)>=threshold
           segonset(i)=1;
       end
       if recsong(i)>=threshold & recsong(i+1)<threshold
           segend(i)=1;
       end
      end
segonsettime = find(segonset>0);
segendtime = find(segend>0);

% discard intersyllable intervals shorter than minseg.
segonsettime2 = segonsettime;
segendtime2 = segendtime;
for i=1:length(segonsettime)-1
    if segonsettime(i+1)-segendtime(i)<minseg/1000*40000
        segonsettime2(i+1)=-1;
        segendtime2(i)=-1;
    end
end

segonsettime2(find(segonsettime2<0))=[];
segendtime2(find(segendtime2<0))=[];

% discard syllables shorter than minseg2.
segonsettime3 = segonsettime2;
segendtime3 = segendtime2;
for i=1:length(segonsettime2)-1
    if segendtime2(i)-segonsettime2(i)<minseg2/1000*40000
        segonsettime3(i)=-1;
        segendtime3(i)=-1;
    end
end

segonsettime3(find(segonsettime3<0))=[];
segendtime3(find(segendtime3<0))=[];

% discard syllables quieter than minamp.
segonsettime4 = segonsettime3;
segendtime4 = segendtime3;
for i=1:length(segonsettime3)-1
    if max(recsong(segonsettime3(i):segendtime3(i)))<minamp
        segonsettime4(i)=-1;
        segendtime4(i)=-1;
    end
end

segonsettime4(find(segonsettime4<0))=[]
segendtime4(find(segendtime4<0))=[]

%save segmetation times

if save == 'save';

filename1 = [songname,'_SegOnsetTime_MinSeg',num2str(minseg),'msec'];
fid = fopen(filename1, 'w');
for n=1:length(segonsettime4)  
    fprintf(fid,'%d\n',segonsettime4(n));
end
fclose(fid);

filename2 = [songname,'_SegEndTime_MinSeg',num2str(minseg),'msec'];
fid = fopen(filename2, 'w');
for n=1:length(segendtime4)
    fprintf(fid,'%d\n',segendtime4(n));
end
fclose(fid);

end

% plot song
t = [1:length(recsong)];
thr1 = ones((length(t)),1).*threshold;
plot(t,recsong,'k');hold;plot(t,thr1,'g')
axis([0,length(recsong),0,max(recsong)])
%title(songname);
%xlabel('Time (seconds)');


% plot setments
for i=1:length(segonsettime4)
    x=[segonsettime4(i),segonsettime4(i)];
    y=[0,max(recsong)];
    plot(x,y,'r');
end

for i=1:length(segendtime4)
    x=[segendtime4(i),segendtime4(i)];
    y=[0,max(recsong)];
    plot(x,y,'b');
end

    