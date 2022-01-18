function [] = PlotSpecificSyllsLogAmplitudeKao(DirName, FileList, FileType, SyllLabels)

Fid = fopen(fullfile(DirName, FileList), 'r');
Files = textscan(Fid, '%s', 'DeLimiter', '\n');
Files = Files{1};
fclose(Fid);

Colours = distinguishable_colors(length(SyllLabels));

figure;

for i = 1:length(SyllLabels),
    MeanAmp{i} = [];
    MeanAmp_Middle{i} = [];
end

for i = 1:length(Files),
    [RawData, Fs] = GetData(DirName, Files{i}, FileType, 1);
    Notes{i} = load(fullfile(DirName, 'ASSLNoteFiles', [Files{i}, '.not.mat']));
    [LogAmplitude] = ASSLCalculateLogAmplitudeKao(RawData, Fs, [], [], []);
    for j = 1:length(SyllLabels),
        Matches = find(Notes{i}.labels == SyllLabels(j));
        if (~isempty(Matches))
            for k = 3:length(Matches),
                subplot(4,4,j);
                hold on;
                Onset = round(Notes{i}.onsets(Matches(k))*Fs/1000);
                Offset = round(Notes{i}.offsets(Matches(k))*Fs/1000);
                
                if (Onset <= 0)
                    Onset = 1;
                end
                
                if (Offset > length(LogAmplitude))
                    Offset = length(LogAmplitude);
                end
                
                plot(LogAmplitude(Onset:Offset), 'k', 'Color', Colours(j,:));
                MeanAmp{j}(end+1) = mean(LogAmplitude(Onset:Offset));
                MeanAmp_Middle{j}(end+1) = mean(LogAmplitude(Onset+3500:Onset+5500));
            end
        end
    end
end

for i = 1:length(SyllLabels),
   subplot(4,4,i);
   axis tight;
   AllAxis(i,:) = axis;
end
% axis tight;
% AllAxis = axis;
% axis([0 1.05*AllAxis(2) 0 1.05*AllAxis(4)]);
% 
for i = 1:length(SyllLabels),
    subplot(4,4,i);
    title(['Syll ', SyllLabels(i)]);
    axis([0 1.05*max(AllAxis(:,2)) 0 1.05*max(AllAxis(:,4))]);
    Temp = axis;
    plot([3500 3500], Temp(3:4), 'k--');
    plot([5500 5500], Temp(3:4), 'k--');
end
% 

disp('Finished');