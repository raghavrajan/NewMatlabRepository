function [WarpedMotifTemplate] = StretchCompressTemplates(S, T, TimeStretch, Normalization)


x = T;
xx = linspace(T(1), T(end), round(TimeStretch*length(T)));
for i = 1:size(S, 1),
    WarpedS(i,:) = spline(x, S(i,:), xx);
%    WarpedS(i,:) = interp1(x, S(i,:), xx);
end

if (Normalization == 1)
    WarpedS = (WarpedS - mean(WarpedS(:)))/std(WarpedS(:));
else
    if (Normalization == 2)
        WarpedS = (WarpedS - median(WarpedS(:)))/mad(WarpedS(:));
    end
end

WarpedMotifTemplate = WarpedS;


