function [Data, Fs] = ReadOKrankData(DirectoryName, FileName, ChanNo)

Data = [];
Fs = [];
PresentDirectory = pwd;

cd(DirectoryName);

[datafid, message] = fopen(FileName, 'r');

[recfid, message] = fopen([FileName, '.rec'], 'r');

if ((recfid) > 0)
    while (~feof(recfid))
        tline = fgetl(recfid);
        if (strfind(tline, 'ai_freq'))
            ColonIndex = find(tline == ':');
            Fs = str2double(tline((ColonIndex + 1):end));
        end

        if (strfind(tline, 'n_ai_chan'))
            ColonIndex = find(tline == ':');
            NoOfChannels = str2double(tline((ColonIndex + 1):end));
        end

        if (strfind(tline, 'n_samples'))
            ColonIndex = find(tline == ':');
            NoOfSamples = str2double(tline((ColonIndex + 1):end));
            break;
        end
    end

    fclose(recfid);

    fseek(datafid, (ChanNo - 1) * 2, 'bof');
    [Data, num_read] = fread(datafid, inf, 'uint16', (NoOfChannels - 1) * 2);
    Data = (Data - 32768) * 10/32768;
    if (num_read ~= NoOfSamples)
        disp(['No of samples does not match that of recfile: ',FileName]);
    end
    
end
    
if (recfid > 0)
    disp(['Read ', num2str(NoOfSamples), ' from ', FileName]);
else
    disp(FileName);
end

if ((datafid) > 0)
    fclose(datafid);
end
   
cd(PresentDirectory);

