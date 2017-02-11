function [] = PlotEllipsoid(Data, Colour, STDLen)

conf = 2*normcdf(STDLen)-1;     %# covers the appropriate std as defined by STDLen
scale = chi2inv(conf,size(Data,2));     %# inverse chi-squared with dof=#dimensions

C = cov(Data);
M = mean(Data);

[U, L] = eig(C);
radii = sqrt(scale*diag(L));
[xc, yc, zc] = ellipsoid(0, 0, 0, radii(1), radii(2), radii(3));
a = kron(U(:,1),xc); b= kron(U(:,2),yc); c = kron(U(:,3),zc);
data = a+b+c;
n = size(data,2);
x = data(1:n,:)+M(1);y = data(n+1:2*n,:)+M(2); z = data(2*n+1:end,:)+M(3);
 
sc = surface(x, y, z, 'FaceColor', Colour, 'EdgeColor', 'none', 'FaceAlpha', 0.1);

disp(['Volume is ', num2str(sqrt(det(C)))]);
