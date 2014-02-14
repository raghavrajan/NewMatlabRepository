function [] = plot_raw_trial_data(DataFile,ChannelNo,NoteFile)

channel_string = strcat('obs',num2str(ChannelNo),'r');
pathname = pwd;
pathname = strcat(pathname,'/');
[SpikeTimes,Fs] = soundin_copy(pathname,DataFile,channel_string);
Times = 0:1/Fs:length(SpikeTimes)/Fs;
Times(end) = [];

Notes = load(NoteFile);

Notes.onsets = Notes.onsets/1000;
Notes.offsets = Notes.offsets/1000;

MotifStarts = find(Notes.labels == 'a');
MotifEnds = find(Notes.labels == 'c');

figure(2);
hold on;

max_y_value = 0;
if (length(MotifStarts) > 5),
    TotalTrials = 5;
else
    TotalTrials = length(MotifStarts);
end

for i = 1:TotalTrials,
    figure(2);
    if (i == 1)
        TrialSpikeTimes = (find((Times > (Notes.onsets(MotifStarts(i)) - 0.5)) & (Times < Notes.offsets(MotifEnds(i)))));
    else
        if ((Notes.onsets(MotifStarts(i)) - Notes.offsets(MotifEnds(i-1))) > 0.5)
            TrialSpikeTimes = (find((Times > (Notes.onsets(MotifStarts(i)) - 0.5)) & (Times < Notes.offsets(MotifEnds(i)))));
        else
            TrialSpikeTimes = (find((Times > (Notes.offsets(MotifEnds(i-1)))) & (Times < Notes.offsets(MotifEnds(i)))));
        end
    end
    plot((Times(TrialSpikeTimes) - Notes.onsets(MotifStarts(i))),(SpikeTimes(TrialSpikeTimes) + max_y_value));
    max_y_value = max_y_value + (max(SpikeTimes(TrialSpikeTimes)) - (min(SpikeTimes(TrialSpikeTimes)))) + 5000;
end

axis tight;

temp = axis;
plot([0 temp(2)],[temp(4) temp(4)],'r','LineWidth',5);
xlabel('Time (sec)','FontSize',14,'FontWeight','bold');
