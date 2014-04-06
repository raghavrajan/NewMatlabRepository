function [Parameters] = MA_ParseParametersFile(ParametersFile)

%==========================================================================
% Function to read and extract all the different parameters from a text
% for a bird that has gone through some treatment
% Raghav Rajan - 29th November 2013
%==========================================================================

Fid = fopen(ParametersFile, 'r');
Temp = textscan(Fid, '%s', 'DeLimiter', '\n');
Temp = Temp{1};
fclose(Fid);

% Extract bird name
Index = find(cellfun(@length, strfind(Temp, 'Bird Name:')));
ColonIndex = find(Temp{Index} == ':');
Parameters.BirdName = strtrim(Temp{Index}(ColonIndex+1:end));

% Extract file type
Index = find(cellfun(@length, strfind(Temp, 'File Type:')));
ColonIndex = find(Temp{Index} == ':');
Parameters.FileType = strtrim(Temp{Index}(ColonIndex+1:end));

% Extract surgery date
Index = find(cellfun(@length, strfind(Temp, 'Surgery Date:')));
ColonIndex = find(Temp{Index} == ':');
Parameters.SurgeryDate = strtrim(Temp{Index}(ColonIndex+1:end));

% First, pre-treatment
% Read # of pre-treatment days
Index = find(cellfun(@length, strfind(Temp, 'No of pre-treatment days:')));
ColonIndex = find(Temp{Index} == ':');
Parameters.NoPreDays = str2double(strtrim(Temp{Index}(ColonIndex+1:end)));

% Read data directories and file lists of all pre-treatment days
for i = 1:Parameters.NoPreDays,
    % First date
    Index = find(cellfun(@length, strfind(Temp, ['Pre-treatment Day ', num2str(i), ' X-label:'])));
    ColonIndex = find(Temp{Index} == ':');
    Parameters.PreDate{i} = strtrim(Temp{Index}(ColonIndex+1:end));
    
    % First data directory
    Index = find(cellfun(@length, strfind(Temp, ['Pre-treatment Day ', num2str(i), ' raw data directory:'])));
    ColonIndex = find(Temp{Index} == ':');
    Parameters.PreDataDir{i} = strtrim(Temp{Index}(ColonIndex+1:end));
    
    % Next pre directed song filelists
    Index = find(cellfun(@length, strfind(Temp, ['Pre-treatment Day ', num2str(i), ' directed song filelist:'])));
    ColonIndex = find(Temp{Index} == ':');
    Parameters.PreDirSongFileList{i} = strtrim(Temp{Index}(ColonIndex+1:end));
    
    % Next pre undirected song filelists
    Index = find(cellfun(@length, strfind(Temp, ['Pre-treatment Day ', num2str(i), ' undirected song filelist:'])));
    ColonIndex = find(Temp{Index} == ':');
    Parameters.PreUnDirSongFileList{i} = strtrim(Temp{Index}(ColonIndex+1:end));
end
   
% Next post treatment
% Read # of post-treatment days
Index = find(cellfun(@length, strfind(Temp, 'No of post-treatment days:')));
ColonIndex = find(Temp{Index} == ':');
Parameters.NoPostDays = str2double(strtrim(Temp{Index}(ColonIndex+1:end)));

% Read data directories and file lists of all post-treatment days
for i = 1:Parameters.NoPostDays,
    % First date
    Index = find(cellfun(@length, strfind(Temp, ['Post-treatment Day ', num2str(i), ' X-label:'])));
    ColonIndex = find(Temp{Index} == ':');
    Parameters.PostDate{i} = strtrim(Temp{Index}(ColonIndex+1:end));
    
    % First data directory
    Index = find(cellfun(@length, strfind(Temp, ['Post-treatment Day ', num2str(i), ' raw data directory:'])));
    ColonIndex = find(Temp{Index} == ':');
    Parameters.PostDataDir{i} = strtrim(Temp{Index}(ColonIndex+1:end));
    
    % Next directed song filelists
    Index = find(cellfun(@length, strfind(Temp, ['Post-treatment Day ', num2str(i), ' directed song filelist:'])));
    ColonIndex = find(Temp{Index} == ':');
    Parameters.PostDirSongFileList{i} = strtrim(Temp{Index}(ColonIndex+1:end));
    
    % Next undirected song filelists
    Index = find(cellfun(@length, strfind(Temp, ['Post-treatment Day ', num2str(i), ' undirected song filelist:'])));
    ColonIndex = find(Temp{Index} == ':');
    Parameters.PostUnDirSongFileList{i} = strtrim(Temp{Index}(ColonIndex+1:end));
end

%==========================================================================
% Extract % remaining HVC
% total
Index = find(cellfun(@length, strfind(Temp, '% total HVC remaining:')));
ColonIndex = find(Temp{Index} == ':');
Parameters.PercentTotalHVCremaining = str2double(strtrim(Temp{Index}(ColonIndex+1:end)));
    
% left
Index = find(cellfun(@length, strfind(Temp, '% left HVC remaining:')));
ColonIndex = find(Temp{Index} == ':');
Parameters.PercentLeftHVCremaining = str2double(strtrim(Temp{Index}(ColonIndex+1:end)));

% right
Index = find(cellfun(@length, strfind(Temp, '% right HVC remaining:')));
ColonIndex = find(Temp{Index} == ':');
Parameters.PercentRightHVCremaining = str2double(strtrim(Temp{Index}(ColonIndex+1:end)));
%==========================================================================

%==========================================================================
% Now get the template files
% motif template
Index = find(cellfun(@length, strfind(Temp, 'Motif template file:')));
ColonIndex = find(Temp{Index} == ':');
Parameters.MotifTemplateFileName = strtrim(Temp{Index}(ColonIndex+1:end));

% syllable templates
Index = find(cellfun(@length, strfind(Temp, 'Syllable templates file:')));
ColonIndex = find(Temp{Index} == ':');
Parameters.SyllableTemplatesFileName = strtrim(Temp{Index}(ColonIndex+1:end));
%==========================================================================

%==========================================================================
% Now get the motif
Index = find(cellfun(@length, strfind(Temp, 'Motif:')));
ColonIndex = find(Temp{Index} == ':');
Parameters.Motif = strtrim(Temp{Index}(ColonIndex+1:end));
%==========================================================================