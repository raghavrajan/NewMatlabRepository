function [] = plot_raw_trial_song_playback_data(DataFile,ChannelNo)

channel_string = strcat('obs',num2str(ChannelNo),'r');
pathname = pwd;
pathname = strcat(pathname,'/');
[SpikeTimes,Fs] = soundin_copy(pathname,DataFile,channel_string);
Times = 0:1/Fs:length(SpikeTimes)/Fs;
Times(end) = [];

NoOfTrials = length(SpikeTimes)/(5.1 * Fs);

for i = 1:NoOfTrials,
    TrialOnsets(i) = ((i-1) * (5.1)) + 2;
end

figure(2);
hold on;

max_y_value = 0;
if (NoOfTrials > 10),
    TotalTrials = 10;
else
    TotalTrials = length(NoOfTrials);
end

for i = 1:TotalTrials,
    figure(2);
    TrialSpikeTimes = (find((Times > (TrialOnsets(i) - 0.5)) & (Times < (TrialOnsets(i) + 1.5))));
    plot((Times(TrialSpikeTimes) - TrialOnsets(i)),(SpikeTimes(TrialSpikeTimes) + max_y_value));
    max_y_value = max_y_value + (max(SpikeTimes(TrialSpikeTimes)) - (min(SpikeTimes(TrialSpikeTimes)))) + 5000;
end

plot([0 0.51],[max_y_value + 100 max_y_value + 100],'r');
plot([0.57 1.1],[max_y_value + 100 max_y_value + 100],'r');

axis tight;