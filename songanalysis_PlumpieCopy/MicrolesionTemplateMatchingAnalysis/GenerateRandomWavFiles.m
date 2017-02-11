function [] = GenerateRandomWavFiles(Duration, Fs, OutputDir, NoofFiles)

% First set the random number generator to the same state by using the same
% seed

s = RandStream('mt19937ar', 'Seed', 0);
RandStream.setDefaultStream(s);

% Check if OutputDir has the a slash at the end - if not add it
FileSep = filesep;
if (OutputDir(end) ~= FileSep)
    OutputDir(end+1) = FileSep;
end

% Check if OutputDir exists - if not make it
if (~exist(OutputDir, 'dir'))
    mkdir(OutputDir);
end

PresentDir = pwd;

cd(OutputDir);

% Now generate random wav files
for i = 1:NoofFiles,
    OutputFileName = ['Random.', num2str(Fs), '.', num2str(i), '.wav'];
    RandomSound = rand(round(Fs*Duration), 1); % these random numbers lie between 0 and 1
    
    % Shift the random sounds to lie between -1 and 1
    RandomSound = -1 + RandomSound*2;
    wavwrite(RandomSound, Fs, OutputFileName);
end

cd(PresentDir);
disp('Finished writing random sound files');
    