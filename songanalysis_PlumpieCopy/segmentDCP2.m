% segment using threshold increasing from threshold 1 to tureshold2 (this script needs to be modified)

function segmentDCP2(songname,threshold1,threshold2,minseg)

[sylonset3,sylend3,recsong] = segment (songname,threshold1,minseg);

for i=threshold1+1:threshold2
    [sylonset4,sylend4] = segment (songname,i,minseg);
    
    sylonset3_4 = sort([sylonset3;sylonset4]);
    sylonset3_4b = sylonset3_4;
    for j=1:length(sylonset3_4)-1
        if sylonset3_4(j+1)-sylonset3_4(j) < minseg/1000*40000;
            sylonset3_4b(j)=-1;
            sylonset3_4b(j+1)=-1;
        end
    end
    sylonset3_4b(find(sylonset3_4b<0))=[];
    sylonset5 = sort([sylonset3;sylonset3_4b]);
    sylonset3 = sylonset5;
    
    sylend3_4 = sort([sylend3;sylend4]);
    sylend3_4b = sylend3_4;
    for j=1:length(sylend3_4)-1
        if sylend3_4(j+1)-sylend3_4(j) < minseg/1000*40000;
            sylend3_4b(j)=-1;
            sylend3_4b(j+1)=-1;
        end
    end
    sylend3_4b(find(sylend3_4b<0))=[];
    sylend5 = sort([sylend3;sylend3_4b]);
    sylend3 = sylend5;
    
end


% plot song
t = [1:length(recsong)]/40000;
thr1 = ones((length(t)),1).*threshold1;
thr2 = ones((length(t)),1).*threshold2;
plot(t,recsong,'k');hold;
plot(t,thr1,'g');plot(t,thr2,'g');
axis([0,length(recsong)/40000,0,max(recsong)])
%threshold2=num2str(threshold);
title(songname);
xlabel('Time (seconds)');


% plot setments
for i=1:length(sylonset5)
    x=[sylonset5(i)/40000,sylonset5(i)/40000];
    y=[0,max(recsong)];
    plot(x,y,'r');
end

for i=1:length(sylend5)
    x=[sylend5(i)/40000,sylend5(i)/40000];
    y=[0,max(recsong)];
    plot(x,y,'b');
end

% Output segmentation time
syllable_onset = sylonset5./40000
syllable_end = sylend5./40000
