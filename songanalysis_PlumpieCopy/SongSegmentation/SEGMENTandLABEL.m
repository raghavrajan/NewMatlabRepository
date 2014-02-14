function [] = SEGMENTandLABEL()

% SEGMENTandLABEL is a program with GUI to segment and label songs
% Raghav - February 2009

% Global Variables
RawSong = [];
Fs = 32000;

SaL = figure('Visible','on','Name','SEGEMENT and LABEL','Position',[360 500 450 285]);
