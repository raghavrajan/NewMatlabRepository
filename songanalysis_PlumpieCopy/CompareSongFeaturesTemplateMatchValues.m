function [] = CompareSongFeaturesTemplateMatchValues(MatchFile, FFFile, TemplateFeaturesFile, TemplateSyllNo, SAPFeaturesFile, TemplateSAPFile)

Syll_FF = load(FFFile);
Matches = load(MatchFile);
TemplateSyll = load(TemplateFeaturesFile);
TemplateFF = TemplateSyll.ffreq{TemplateSyllNo, 2};

SAPFeats = load(SAPFeaturesFile);
TemplateSAPFeats = load(TemplateSAPFile);
TemplateFreq = TemplateSAPFeats.Feats.Freq(TemplateSyllNo);
TemplateAmp = TemplateSAPFeats.Feats.Amp(TemplateSyllNo);
TemplateEntropy = TemplateSAPFeats.Feats.Entropy(TemplateSyllNo);

MatchValue = [];
FFs = [];

for i = 1:size(Syll_FF.ffreq),
    FileIndex = find(cellfun(@length, strfind(Matches.UnDirBout.FileName, Syll_FF.ffreq{i,1}(1:end-5))));
    if (~isempty(FileIndex))
        FFs = [FFs Syll_FF.ffreq{i,2}];
        MatchValue = [MatchValue max(Matches.UnDirBout.MaxBoutSeqMatch{FileIndex}(find((Matches.UnDirBout.T{FileIndex} > (Syll_FF.ffreq{i,14} - 0.01)) & (Matches.UnDirBout.T{FileIndex} < (Syll_FF.ffreq{i,14} + 0.01)))))];
    end
end
Entropy = SAPFeats.Feats.Entropy;
MeanFreq = SAPFeats.Feats.Freq;
Amplitude = SAPFeats.Feats.Amp;

figure;
subplot(2,2,1);
plot((FFs - TemplateFF)/TemplateFF * 100, MatchValue, 'k+');
xlabel('% change in ff', 'FontSize', 12, 'FontName', 'Arial');
ylabel('Match to template', 'FontSize', 12, 'FontName', 'Arial');
axis tight;

subplot(2,2,2);
plot((MeanFreq - TemplateFreq)/TemplateFreq * 100, MatchValue, 'k+');
xlabel('% change in mean frequency', 'FontSize', 12, 'FontName', 'Arial');
ylabel('Match to template', 'FontSize', 12, 'FontName', 'Arial');
axis tight;

subplot(2,2,3);
plot((Entropy - TemplateEntropy)/TemplateEntropy * 100, MatchValue, 'k+');
xlabel('% change in entropy', 'FontSize', 12, 'FontName', 'Arial');
ylabel('Match to template', 'FontSize', 12, 'FontName', 'Arial');
axis tight;

subplot(2,2,4);
plot((Amplitude - TemplateAmp)/TemplateAmp * 100, MatchValue, 'k+');
xlabel('% change in logamplitude', 'FontSize', 12, 'FontName', 'Arial');
ylabel('Match to template', 'FontSize', 12, 'FontName', 'Arial');
axis tight;

[r, p] = corrcoef(abs((FFs - TemplateFF)/TemplateFF * 100)', MatchValue');
disp(['Correlation co-efficient between ff and match value is ', num2str(r(1,2)), ' and the p-value is ', num2str(p(1,2))]);

[r, p] = corrcoef(abs((Entropy - TemplateEntropy)/TemplateEntropy * 100)', MatchValue');
disp(['Correlation co-efficient between entropy and match value is ', num2str(r(1,2)), ' and the p-value is ', num2str(p(1,2))]);

[r, p] = corrcoef(abs((Amplitude - TemplateAmp)/TemplateAmp * 100)', MatchValue');
disp(['Correlation co-efficient between amplitude and match value is ', num2str(r(1,2)), ' and the p-value is ', num2str(p(1,2))]);

[r, p] = corrcoef(abs((MeanFreq - TemplateFreq)/TemplateFreq * 100)', MatchValue');
disp(['Correlation co-efficient between mean frequency and match value is ', num2str(r(1,2)), ' and the p-value is ', num2str(p(1,2))]);

disp('Finished');