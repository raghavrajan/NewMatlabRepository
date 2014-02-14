function [] = CalculateFFsAroundBursts(DirectoryName)

Files = dir('*SSAOutput*');

for i = 1:length(Files),
    load(Files(i).name);
    DataDir = AnalysisOutput.RawDataDirectory;
    RecDir = AnalysisOutput.RecFileDirectory;
    