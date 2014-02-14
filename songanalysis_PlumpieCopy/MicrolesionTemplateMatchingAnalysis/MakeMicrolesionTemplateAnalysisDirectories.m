function [handles] = MakeMicrolesionTemplateAnalysisDirectories(handles, DirString, VariableString, ComparisonType)

FileSep = filesep;

if (strfind(ComparisonType, 'Normal'))
    OutputDir = handles.OutputDir;
else
    if (handles.OutputDir(end) == FileSep)
        OutputDir = [handles.OutputDir, ComparisonType, FileSep];
    else
        OutputDir = [handles.OutputDir, FileSep, ComparisonType, FileSep];
    end
end

eval(['handles.', ComparisonType, '_', VariableString, 'TemplateMatchOutputFilesDir = ''', OutputDir, get(eval(['handles.', VariableString, 'TemplateFileNameText']), 'String'), '.', DirString, FileSep, ''';']);
mkdir(eval(['handles.', ComparisonType, '_', VariableString, 'TemplateMatchOutputFilesDir']));

eval(['handles.Dir', ComparisonType, '_', VariableString, 'TemplateMatchOutputFilesDir = ''', OutputDir, get(eval(['handles.', VariableString, 'TemplateFileNameText']), 'String'), '.', DirString, FileSep, 'Dir'';']);
mkdir(eval(['handles.Dir', ComparisonType, '_', VariableString, 'TemplateMatchOutputFilesDir']));

eval(['handles.UnDir', ComparisonType, '_', VariableString, 'TemplateMatchOutputFilesDir = ''', OutputDir, get(eval(['handles.', VariableString, 'TemplateFileNameText']), 'String'), '.', DirString, FileSep, 'UnDir'';']);
mkdir(eval(['handles.UnDir', ComparisonType, '_', VariableString, 'TemplateMatchOutputFilesDir']));

