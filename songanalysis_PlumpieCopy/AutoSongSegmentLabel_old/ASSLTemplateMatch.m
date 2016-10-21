function [Match] = ASSLTemplateMatch(TempS, MotifTemplate,  StretchValues)

WarpIndex = 1;
for i = 1:length(MotifTemplate),
    WMotif = MotifTemplate(i).MotifTemplate;
    if (size(WMotif,2) > size(TempS,2))
        continue;
    end
    
    TempMeanSTD = CalculateMeanSTDforSpectralMatch(TempS(1:size(TempS,1)*size(TempS,2)), size(WMotif,1)*size(WMotif,2), (size(TempS,2) - size(WMotif,2) + 1), size(WMotif,1));

    WinMean = TempMeanSTD(1:length(TempMeanSTD)/2);
    WinSTD = TempMeanSTD((length(TempMeanSTD)/2 + 1):end);
    
    [TempMatch{WarpIndex}] = CalTemplateMatch(WMotif, TempS, WinMean, WinSTD);
    TempMatch{WarpIndex} = TempMatch{WarpIndex} * size(WMotif,1)*size(WMotif,2);
	
    WarpIndex = WarpIndex + 1;
end

for MatchNo = 1:length(TempMatch),
    Match(MatchNo,:) = TempMatch{MatchNo}(1:length(TempMatch{end}));
end
Match = max(Match);
