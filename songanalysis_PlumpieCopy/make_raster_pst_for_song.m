function [] = make_raster_pst_for_song(SpikeFile,NoteFile,bin_size)

SpikeTimes = load(SpikeFile);
Notes = load(NoteFile);

Notes.onsets = Notes.onsets/1000;
Notes.offsets = Notes.offsets/1000;

MotifStarts = find(Notes.labels == 'a');
MotifEnds = find(Notes.labels == 'c');

figure(3);
hold on;
% figure(4);
% hold on;
for i = 1:length(MotifStarts),
    figure(3);
    if (i == 1)
        TrialSpikeTimes = SpikeTimes(find((SpikeTimes > (Notes.onsets(MotifStarts(i)) - 0.5)) & (SpikeTimes < Notes.offsets(MotifEnds(i))))) - Notes.onsets(MotifStarts(i));
    else
        if ((Notes.onsets(MotifStarts(i)) - Notes.offsets(MotifEnds(i-1))) > 0.5)
            TrialSpikeTimes = SpikeTimes(find((SpikeTimes > (Notes.onsets(MotifStarts(i)) - 0.5)) & (SpikeTimes < Notes.offsets(MotifEnds(i))))) - Notes.onsets(MotifStarts(i));
        else
            TrialSpikeTimes = SpikeTimes(find((SpikeTimes > (Notes.offsets(MotifEnds(i-1)))) & (SpikeTimes < Notes.offsets(MotifEnds(i))))) - Notes.onsets(MotifStarts(i));
        end
    end
    y_value = ones(length(TrialSpikeTimes),1) * i/5;
    plot(TrialSpikeTimes,y_value,'w+');
    marker_string = repmat('|',size(TrialSpikeTimes,1),size(TrialSpikeTimes,2));
    text(TrialSpikeTimes,y_value,marker_string,'HorizontalAlignment','center','VerticalAlignment','middle');
    
%     figure(4);
%     for j = 1:3,
%         SyllableLength = [Notes.onsets(MotifStarts(i) - 1 + j) Notes.offsets(MotifStarts(i) - 1 + j)] - Notes.onsets(MotifStarts(i));
%         y_value = ones(length(SyllableLength),1)*i/5;
%         plot(SyllableLength,y_value,'r');
%     end
    
end

figure(3)
xlabel('Time (sec)','FontWeight','bold','FontSize',12);
ylabel('Rendition no','FontWeight','bold','FontSize',12);
axis tight;
temp = axis;
plot([0 temp(2)],[temp(4) + 0.2 temp(4) + 0.2],'r','LineWidth',4);
xlabel('Time (sec)','FontSize',14,'FontWeight','bold');
axis([temp(1) temp(2) 0 (temp(4) + 0.3)]);