function [rawsong, Fs] = read_okrank_data(sound_file, chan_spec)

rawsong = [];
Fs = [];

[datafid, message] = fopen(sound_file, 'r');

[recfid, message] = fopen([sound_file, '.rec'], 'r');

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

    ChanNo = str2double(chan_spec) + 1;
    fseek(datafid, (ChanNo - 1) * 2, 'bof');
    [rawsong, num_read] = fread(datafid, inf, 'uint16', (NoOfChannels - 1) * 2);
    rawsong = (rawsong - 32768) * 10/32768;
    if (num_read ~= NoOfSamples)
        disp(['No of samples does not match that of recfile: ',FileName]);
    end
end

if ((datafid) > 0)
    fclose(datafid);
end
   



