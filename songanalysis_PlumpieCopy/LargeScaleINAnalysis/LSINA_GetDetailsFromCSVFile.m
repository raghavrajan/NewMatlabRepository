function [HeaderLine, BirdDetails] = LSINA_GetDetailsFromCSVFile(BirdDetailsTextFile)

% The text file is a csv file with semi-colon as the delimiter
% So, first load up all the lines of the file

Fid = fopen(BirdDetailsTextFile, 'r');
TextLines = textscan(Fid, '%s', 'DeLimiter', '\n');
TextLines = TextLines{1};
fclose(Fid);

% Now for each line, get the data using semi-colon as the delimiter
for i = 1:length(TextLines),
    IndividualLines = textscan(TextLines{i}, '%s', 'DeLimiter', ';');
    if (i== 1)
        HeaderLine = IndividualLines{1};
    else
        BirdDetails{i-1} = IndividualLines{1};
    end
end
