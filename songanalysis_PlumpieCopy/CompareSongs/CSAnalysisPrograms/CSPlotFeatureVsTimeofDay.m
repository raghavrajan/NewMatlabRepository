function [] = CSPlotFeatureVsTimeofDay(CSData)

PresentDir = pwd;

Syllable = inputdlg('Choose the syllable label that you want to plot (Please choose only one syllable)');

Feature = inputdlg('Choose the feature that you want to plot vs. time of day (Please choose only one feature)');

FileSep = filesep;

for i = 1:CSData.NoofDays,
    if (CSData.Data{i}.DirName(end) ~= FileSep)
        TimesFile = dir([CSData.Data{i}.DirName, FileSep, '*.times']);
    else
        TimesFile = dir([CSData.Data{i}.DirName, '*.times']);
    end
    
    if (~isempty(TimesFile))
        Fid = fopen(TimesFile(1).name, 'r');
        Temp = textscan(Fid, '%s', 'DeLimiter', '\n');
        Temp = Temp{1};
        fclose(Fid);
    else    
        cd(CSData.Data{i}.DirName);
        Files = dir('*.wav');
        for j = 1:length(Files),
            Temp{j} = Files(j).name;
        end
        cd(PresentDir);
    end
    for j = 1:length(CSData.Data{i}.FileName),
        FileIndex = find(cellfun(@length, strfind(Temp, CSData.Data{i}.FileName{j})));
        TempLocationIndex = strfind(Temp{FileIndex}, 'Recorded:');
        if (~isempty(TempLocationIndex))
            TempFileTime = Temp{FileIndex}((TempLocationIndex + 21):(TempLocationIndex + 28));
            CSData.Data{i}.FileTime{j} = str2double(TempFileTime(1:2)) + str2double(TempFileTime(4:5))/60 + str2double(TempFileTime(7:8))/3600; 
        else
            if (FileIndex ~= length(Temp))
                TempFileTime = Temp{FileIndex + 1}((end - 9):(end - 4));
            else
                TempFileTime = Temp{FileIndex}((end - 9):(end - 4));
            end
            CSData.Data{i}.FileTime{j} = str2double(TempFileTime(1:2)) + str2double(TempFileTime(3:4))/60 + str2double(TempFileTime(5:6))/3600; 
        end
    end
end
Colors = 'rgbcmk';
Symbols = 'o+d';

FeatureValueMatrix = [];

for i = 1:CSData.NoofDays,
    FeatureIndex = find(cellfun(@length, strfind(CSData.Data{i}.ToBeUsedFeatures, Feature{1})));
    
    Indices = find(CSData.Data{i}.SyllIndexLabels == Syllable{1});
    for j = 1:length(Indices),
        FeatureValueMatrix(end+1, :) = [CSData.Data{i}.FileTime{CSData.Data{i}.SyllIndices(Indices(j),1)} CSData.Data{i}.FeatValues(Indices(j), FeatureIndex) i];
    end
end
figure;
for i = 1:FeatureValueMatrix(end, 3),
    Indices = find(FeatureValueMatrix(:,3) == i);
    plot(FeatureValueMatrix(Indices, 1), FeatureValueMatrix(Indices, 2), [Colors(mod(i-1, length(Colors)) + 1), Symbols(ceil(i/length(Colors)))]);
    hold on;
    Indices = Indices(find(~isnan(FeatureValueMatrix(Indices, 2))));
    
    errorbar(FeatureValueMatrix(Indices(end), 1) + 1, mean(FeatureValueMatrix(Indices, 2)), std(FeatureValueMatrix(Indices, 2)), [Colors(mod(i-1, length(Colors)) + 1), Symbols(ceil(i/length(Colors)))], 'MarkerSize', 9, 'LineWidth', 2);
end

xlabel('Time of day (hour)', 'FontSize', 16);

ylabel(Feature{1}, 'FontSize', 16);

title(['Syllable ', Syllable{1}], 'FontSize', 16);
disp('Finished plotting data');