% segment a song with a threshold

function [sylonset3,sylend3,recsong] = segmentDCP(songname,threshold,minseg,minamp);

%minseg is minmum duration of segment (msec)
%minseg2 = minseg/1000*40000;

%minamp is minmum amplitude of segment. a note in which max amplitude is smaller than minamp is discarded.

song = readDCP40k(songname);

lvnwinsize=40; % lvn window size: 1 msec for 40k song
stim = song-mean(song);
win = lvnwinsize+1:length(stim);
lvnsong = zeros(length(stim),1); 

for i = (lvnwinsize+1):length(stim)
    lvnsong(i-lvnwinsize/2) = std(stim(i-lvnwinsize:i));
end

figure
 plot(stim);
 hold on
 plot(lvnsong,'r');
 xlabel('time');
 ylabel('amplitude');
 legend('original song','LVN normalizer');
 title(songname)

%rectify the song signal
%recsong = zeros((length(song)),1);
%for i=1:length(song);
%    if song(i)<0
%        recsong(i)=song(i)*(-1);
%    elseif song(i)>=0
%        recsong(i)=song(i);
%    end
%end


%segmentation
sylonset = zeros((length(lvnsong)),1);
sylend = zeros((length(lvnsong)),1);
for i=1:length(lvnsong)-1;
    if lvnsong(i)<threshold & lvnsong(i+1)>=threshold
        sylonset(i) = 1;
    end
    if lvnsong(i)>=threshold & lvnsong(i+1)<threshold
        sylend(i) = 1;
    end
end

sylonsettime = find(sylonset>0);
sylendtime = find(sylend>0);


%duiscard segments shorter than 2 msec
    sylonsettime2 = sylonsettime;
    sylendtime2 = sylendtime;
    for i=1:length(sylonsettime)-1
        if sylonsettime(i+1)-sylendtime(i) < 2/1000*40000;
            sylendtime2(i)=-1;
            sylonsettime2(i+1)=-1;
        end
    end

    sylonsettime2(find(sylonsettime2<0))=[];
    sylendtime2(find(sylendtime2<0))=[];

    
% discard segment in which max amplitude is smaller than minamp
    sylonsettime3 = sylonsettime2;
    sylendtime3 = sylendtime2;

    for i=1:length(sylonsettime2);
        if max(lvnsong(sylonsettime2(i):sylendtime2(i)))<minamp
            sylonsettime3(i)=-1;
            sylendtime3(i)=-1;
        end 
    end

    sylonsettime3(find(sylonsettime3<0))=[];
    sylendtime3(find(sylendtime3<0))=[];

    
%duiscard segments shorter than minseg    
for j=2:minseg    
    
    sylonsettime4 = sylonsettime3;
    sylendtime4 = sylendtime3;
    for i=1:length(sylonsettime3)-1
        if sylonsettime3(i+1)-sylendtime3(i) < j/1000*40000;
            sylendtime4(i)=-1;
            sylonsettime4(i+1)=-1;
        end
    end

    sylonsettime4(find(sylonsettime4<0))=[];
    sylendtime4(find(sylendtime4<0))=[];


%    sylonsettime3 = sylonsettime2;
%    sylendtime3 = sylendtime2;
%    for i=1:length(sylonsettime2)-1
%        if sylendtime2(i)-sylonsettime2(i) < j/1000*40000;
%            sylendtime3(i)=-1;
%            sylonsettime3(i)=-1;
%        end
%    end
%
%    sylonsettime3(find(sylonsettime3<0))=[];
%    sylendtime3(find(sylendtime3<0))=[];

    sylonsettime3 = sylonsettime4;
    sylendtime3 = sylendtime4;

end
  
sylonsettime4 = sylonsettime3;
sylendtime4 = sylendtime3;

% plot song
figure
t = [1:length(song)]/40000;
thr = ones((length(t)),1).*threshold;
plot(t,lvnsong,'k');hold
plot(t,thr,'g');
axis([0,length(song)/40000,0,max(song)])
threshold2=num2str(threshold);
title(songname);
xlabel('Time (seconds)');


% plot setments
for i=1:length(sylonsettime4)
    x=[sylonsettime4(i)/40000,sylonsettime4(i)/40000];
    y=[0,max(song)];
    plot(x,y,'r');
end

for i=1:length(sylendtime4)
    x=[sylendtime4(i)/40000,sylendtime4(i)/40000];
    y=[0,max(song)];
    plot(x,y,'b');
end

% Output segmentation time
syllable_onset = sylonsettime4./40000
syllable_end = sylendtime4./40000

%Output for segment2
%sylonset4 = zeros((length(recsong)),1);
%sylend4 = zeros((length(recsong)),1);
%for i=1:length(sylonset3)
%    sylonset4(sylonset3(i))=1;
%end
%for i=1:length(sylend3)
%    sylend4(sylend3(i))=1;
%end