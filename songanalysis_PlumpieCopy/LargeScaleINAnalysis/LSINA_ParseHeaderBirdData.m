function [BirdParameters] = LSINA_ParseHeaderBirdData(HeaderLine, BirdDetails)

% Go through the header line and check what it corresponds to and then
% assign bird detail variables appropriately

for i = 1:length(HeaderLine),
    ParameterName = HeaderLine{i}(find(HeaderLine{i} ~= ' '));
    switch (ParameterName)
        case {'Continuousdata', 'Interboutinterval', 'InterINinterval', 'SerialNo', 'Fulldaydata', 'Templateindex', 'Onsettimeinms', 'Offsettimeinms', 'SongChanNo', 'NeuronNo', 'RevisedNeuronNo', 'FinalNeuronNo', 'SpikeChanNo', 'Bout', 'BoutOnset', 'BoutOffset', 'Tutor', 'NestNo'}
            for j = 1:length(BirdDetails),
                eval(['BirdParameters(', num2str(j), ').', ParameterName, '= str2double(BirdDetails{', num2str(j), '}{', num2str(i), '});']);
            end
            
        case 'SpikeClusterNos'
            for j = 1:length(BirdDetails),
                Temp = textscan(BirdDetails{j}{i}, '%f', 'DeLimiter', ',');
                eval(['BirdParameters(', num2str(j), ').', ParameterName, '= Temp{1};']);
            end
            
        case 'Optionalclusterstotest'
            for j = 1:length(BirdDetails),
                Temp = textscan(BirdDetails{j}{i}, '%f', 'DeLimiter', ',');
                eval(['BirdParameters(', num2str(j), ').', ParameterName, '= Temp{1};']);
            end
            
        case 'CommonMotifs'
            for j = 1:length(BirdDetails),
                eval(['Temp = BirdDetails{', num2str(j), '}{', num2str(i), '};']);
                if (~isempty(find(Temp == ',')))
                    Temp = textscan(Temp, '%s', 'DeLimiter', ',');
                    Temp = Temp{1};
                    for k = 1:length(Temp),
                        Temp{k} = Temp{k}(find(Temp{k} ~= ' '));
                    end
                    BirdParameters(j).CommonMotifs = Temp;
                else
                    BirdParameters(j).CommonMotifs{1} = Temp;
                end
            end
            
        otherwise
            for j = 1:length(BirdDetails),
                eval(['BirdParameters(', num2str(j), ').', ParameterName, '= BirdDetails{', num2str(j), '}{', num2str(i), '};']);
            end
    end
end