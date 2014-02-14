function [WarpedMotifTemplate] = StretchCompressTemplatesSAPFeatures(RawData, Fs, TimeStretch, FreqStretch)

S = RawData;
T = (1:1:length(RawData))/Fs;

x = T;
xx = linspace(T(1), T(end), round(TimeStretch*length(T)));
for i = 1:size(S, 1),
    WarpedS(i,:) = spline(x, S(i,:), xx);
    WarpedS(i,:) = (WarpedS(i,:) - mean(WarpedS(i,:)))/std(WarpedS(i,:));
%    WarpedS(i,:) = interp1(x, S(i,:), xx);
end

%WarpedS = (WarpedS - mean(WarpedS(:)))/std(WarpedS(:));
WarpedMotifTemplate = WarpedS;

