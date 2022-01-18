function [ROCArea] = PlotGaussianSmoothedMotifActivity(DirectoryName, FileList, FileType, GaussWidth, MotifStartSyll, INLabel, LagToSound)

PresentDir = pwd;

Fid = fopen(FileList, 'r');
Files = textscan(Fid, '%s', 'DeLimiter', '\n');
Files = Files{1};
fclose(Fid);

AmplitudeFs = 1000; % resolution of 1ms for syllable durations

Width = GaussWidth; % in seconds - Gaussian width or std
GaussianLen = 3; % length of gaussian

XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (AmplitudeFs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * AmplitudeFs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * AmplitudeFs) * (Width * AmplitudeFs)));

% For an exponential function
XExponential = (0:1:(3*Width)*AmplitudeFs)/AmplitudeFs;
ExponentialWin = exp(-XExponential/Width);
ExponentialWin = ExponentialWin/sum(ExponentialWin);

StartVal = [];
FirstStartVal = [];
INVal = [];

for i = 1:length(Files)/2,
    [RawData, Fs] = GetData(DirectoryName, Files{i}, FileType, 0);
    Notes = load(fullfile(DirectoryName, 'ASSLNoteFiles', [Files{i}, '.not.mat']));
    SyllDurationPlot = zeros(length(0:1/AmplitudeFs:length(RawData)/Fs), 1);
    Onsets = Notes.onsets/1000;
    Offsets = Notes.offsets/1000;
    
    for j = 1:length(Onsets),
        SyllDurationPlot(ceil(Onsets(j)*AmplitudeFs):ceil(Onsets(j)*AmplitudeFs)) = 1;
    end
    
    % SmoothedSyllDurationPlot = conv(SyllDurationPlot, GaussWin, 'same');
    SmoothedSyllDurationPlot = conv(SyllDurationPlot, ExponentialWin);
    SmoothedSyllDurationPlot = SmoothedSyllDurationPlot(1:length(SyllDurationPlot));
    
    Matches = find(Notes.labels == MotifStartSyll);
    for j = Matches(:)',
        if (j == Matches(1))
            FirstStartVal(end+1) = SmoothedSyllDurationPlot(round((Onsets(j) - LagToSound) * AmplitudeFs));
        else
            StartVal(end+1) = SmoothedSyllDurationPlot(round((Onsets(j) - LagToSound) * AmplitudeFs));
        end
    end
    
    INs = find(Notes.labels == INLabel);
    for j = INs(:)',
        INVal(end+1) = SmoothedSyllDurationPlot(round((Onsets(j) - LagToSound) * AmplitudeFs));
    end
end

CriterionVal = prctile(FirstStartVal(:), 5);
%disp(['Fraction of IN values >= CriterionVal = ', num2str(length(find(INVal(:) >= CriterionVal))/length(INVal))]);


% Now to predict in the rest of the data

PredictionSetStartVal = [];
PredictionSetINVal = [];

for i = (1+length(Files)/2):length(Files),
    [RawData, Fs] = GetData(DirectoryName, Files{i}, FileType, 0);
    Notes = load(fullfile(DirectoryName, 'ASSLNoteFiles', [Files{i}, '.not.mat']));
    SyllDurationPlot = zeros(length(0:1/AmplitudeFs:length(RawData)/Fs), 1);
    Onsets = Notes.onsets/1000;
    Offsets = Notes.offsets/1000;
    
    for j = 1:length(Onsets),
        SyllDurationPlot(ceil(Onsets(j)*AmplitudeFs):ceil(Offsets(j)*AmplitudeFs)) = 1;
    end
    
    % SmoothedSyllDurationPlot = conv(SyllDurationPlot, GaussWin, 'same');
    SmoothedSyllDurationPlot = conv(SyllDurationPlot, ExponentialWin);
    SmoothedSyllDurationPlot = SmoothedSyllDurationPlot(1:length(SyllDurationPlot));
    
    Matches = find(Notes.labels == MotifStartSyll);
    for j = Matches(:)',
        PredictionSetStartVal(end+1) = SmoothedSyllDurationPlot(round((Onsets(j) - LagToSound) * AmplitudeFs));
    end
    
    INs = find(Notes.labels == INLabel);
    for j = INs(:)',
        PredictionSetINVal(end+1) = SmoothedSyllDurationPlot(round((Onsets(j) - LagToSound) * AmplitudeFs));
    end

end

%disp(['Fraction of Prediction set start syllables that are above criterion = ', num2str(length(find(PredictionSetStartVal >= CriterionVal))/length(PredictionSetStartVal))]);
%disp(['Fraction of Prediction set INs that are above criterion = ', num2str(length(find(PredictionSetINVal >= CriterionVal))/length(PredictionSetINVal))]);

AllINVals = [PredictionSetINVal(:); INVal(:)];
FractionINsAboveCriterion = length(find(AllINVals >= CriterionVal))/length(AllINVals);
AllStartVals = [FirstStartVal(:); StartVal(:); PredictionSetINVal(:)];

Index = 1;
for Threshold = 0:0.05:1,
    FalsePositive(Index) = length(find(AllINVals >= Threshold))/(length(AllINVals));
    TruePositive(Index) = length(find(AllStartVals >= Threshold))/(length(AllStartVals));
    Index = Index + 1;
end

% figure
% plot(FalsePositive, TruePositive, 'bs-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
% hold on;
% plot([0 1], [0 1], 'k--', 'LineWidth', 1);
ROCArea = polyarea([FalsePositive(:); 1], [TruePositive(:); 0]);
%disp(['ROC area = ', num2str(ROCArea)]);