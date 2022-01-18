function [] = ExperimentingWithLME()

% This is my script for experimenting and getting some intuition about Linear Mixed Effects Models.


% First, I'm going to try with an n=10 birds and 4 distances. I will
% generate random numbers for the response variable for each of these
% distances and take different slopes and different intercepts for each of
% the birds from a gaussian

%% Different slopes and different intercepts for each bird
n = 10; % Number of birds
N_Mean = 20; % Mean number of trials
N_Sigma = 4; % Sigma for number of trials

Distances = {'L0' 'L1' 'L2' 'L3'}; % 4 different distances
InterceptMean_Sigma = [2 1];
SlopeMean_Sigma = [0.5 2];

ResponseVar = [];
Groups = [];
for i = 1:n,
    Intercept(i) = normrnd(InterceptMean_Sigma(1), InterceptMean_Sigma(2));
    Slope(i) = normrnd(SlopeMean_Sigma(1), SlopeMean_Sigma(2));
    
    for j = 1:length(Distances),
        NumTrials = ceil(normrnd(N_Mean, N_Sigma));
        if (NumTrials < 3)
            NumTrials = 3;
        end
        ResponseVar = [ResponseVar; rand(NumTrials,1) + Intercept(i) + Slope(i)*j];
        Groups = [Groups; [ones(NumTrials,1)*i ones(NumTrials,1)*j]];
    end
end

DataTable = table(ResponseVar, Distances(Groups(:,2))', Groups(:,1), 'VariableNames', {'Response', 'Distance', 'Bird'});

% Model this with a lme with distance as fixed effect and bird slope and
% bird intercept as random effects
LME1 = fitlme(DataTable, 'Response ~ 1 + Distance'); % Only distance as fixed effect
LME2 = fitlme(DataTable, 'Response ~ 1 + Distance + (1|Bird)'); % Distance as fixed effect and random intercept for each bird
LME3 = fitlme(DataTable, 'Response ~ 1 + Distance + (1|Bird) + (-1 + Distance|Bird)'); % Distance as fixed effect and random intercept and slope for each bird

disp('Finished');
    