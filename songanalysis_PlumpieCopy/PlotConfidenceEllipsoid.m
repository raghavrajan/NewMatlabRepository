function [] = PlotConfidenceEllipsoid(Data, Colour, STDLen)

conf = 2*normcdf(STDLen)-1;     %# covers the appropriate std as defined by STDLen
scale = chi2inv(conf,2);     %# inverse chi-squared with dof=#dimensions

Cov = cov(Data) * scale;

Mu = mean(Data);
[V D] = eig(Cov);
[D order] = sort(diag(D), 'descend');
D = diag(D);
V = V(:, order);

t = linspace(0,2*pi,100);
phi = linspace(0, pi, 100);
e = [sin(phi).*cos(t) ; sin(phi).*sin(t); cos(phi)];        %# unit sphere
VV = V*sqrt(D);               %# scale eigenvectors
e = bsxfun(@plus, VV*e, Mu'); %#' project sphere back to orig space

%# plot cov and major/minor axes
fill3(e(1,:), e(2,:), e(3,:), Colour, 'FaceColor', Colour, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
