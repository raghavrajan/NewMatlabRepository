function [] = MakeWaveFile(pathname,filename,varargin)

ChannelNo = 0;

% Get Song Data from observer file
channel_string = strcat('obs',num2str(0),'r');
[rawsong,Fs] = soundin_copy(pathname,filename,channel_string);

% Convert to uV - 5V on the data acquisition is 32768
rawsong = rawsong * 1/32768;

% Time axis
time = 0:1/32000:(length(rawsong)/32000);
time(end) = [];

if (length(varargin) > 0)
    start_time = varargin{1};
    end_time = varargin{2};
    start_time_index = find(time <= start_time,1,'last');
    end_time_index = find(time <= end_time,1,'last');
else
    start_time_index = 1;
    end_time_index = length(time);
end

% Plot song spectrogram

time = time(start_time_index:end_time_index);
rawsong = rawsong(start_time_index:end_time_index);

temp = resample(rawsong,44100,Fs);
OutputFileName = [filename,'.',num2str(start_time),'_',num2str(end_time),'.wav'];
wavwrite(temp,44100,16,OutputFileName);

disp('Finished');