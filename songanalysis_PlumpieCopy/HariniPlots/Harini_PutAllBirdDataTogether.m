function [AllBirdFinalStats] = Harini_PutAllBirdDataTogether(FinalBirdStats, InterBoutInterval)

% =============================================================================================
% This is a script to combine all of Harini's data together separately for
% # of INs, # of motifs/bout and first motif duration - these are all the
% bout level features. I will do the same for syllable level features
% afterwards
% There will be an array with the data for each of the features. In
% addition there will be grouping arrays to indicate the following:
% 1) Bird identity
% 2) Time of day
% 3) Recording day #
% 4) Session # in the day
% 5) Overall session #
% 6) Type of song based on video scoring
% 7) Female response again based on video scoring
% 8) Distance 
% 9) Microphone type
% 10) Bout # within the session
%                                                   Raghav Rajan 22.05.2019
% =============================================================================================

%% Indices for various different variables and conditions

% Indices for distances
DistanceLabels = {'L0' 'L1' 'L2' 'L3' 'L4' 'UN'};
Distances = [0 20 60 110 165 185];

% Indices for Microphones
MicrophoneTypes = {'HF' 'BP' 'MM'};
MicrophoneTypeIndices = [1 2 3];

% Indices for song type
SongTypes = {'D' 'DUN' 'UN' 'NA'};  % based on video scoring
SongTypeIndices = [1 2 3 NaN];

% Indices for female response
FemaleResponseTypes = {'R' 'NR' 'NA'};  % based on observation and sometimes video scoring
SongTypeIndices = [1 2 NaN];

%% Actual loop for putting together all the data
for i = 1:length(FinalBirdStats),
    