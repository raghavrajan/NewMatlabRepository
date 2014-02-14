function [Match] = TempRTemplateMatch(TempS, WMotif)

Match = 1./(sum(bsxfun(@minus, TempS, WMotif).^2)/size(TempS,1));