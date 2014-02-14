function [] = Resample_and_Convert_to_wav(DirectoryName,BirdName, FileType)

if (ispc)
    if ~(DirectoryName(end) == '\')
        DirectoryName(end + 1) = '\';
    end
else
    if ~(DirectoryName(end) == '/')
        DirectoryName(end + 1) = '/';
    end
end

cd(DirectoryName);

FilesToConvert = dir([BirdName, '*']);

for i = 1:length(FilesToConvert),
    try
        if (strfind(FileType, 'obs'))
            [rawsong,Fs] = soundin_copy(DirectoryName,FilesToConvert(i).name,'obs0r');
            rawsong = rawsong/32768;
        else
            if (strfind(FileType, 'okrank'))
                [rawsong,Fs] = ReadOKrankData(DirectoryName, FilesToConvert(i).name, 1);
                rawsong = rawsong/10;
            end
        end
        
        temp = resample(rawsong,44100,Fs);
        OutputFileName = [FilesToConvert(i).name,'.wav'];
        wavwrite(temp,44100,16,OutputFileName);
        disp(['Finished converting ',FilesToConvert(i).name]);
    catch
        disp(['Could not convert ',FilesToConvert(i).name]);
    end
end

disp(['Finished converting ',num2str(i),' files']);
