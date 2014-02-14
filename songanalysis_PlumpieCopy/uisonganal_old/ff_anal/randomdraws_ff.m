function [ratio_alonevdir,percent_greater,number_greater,ndraws]=randomdraws_ff(ff_alone,ff_dir,ndraws);

%this function will randomly assign ff values to a dir/undir category and
%calculate the ratio of the variance for the two randomly assigned categories
%it will compare this with the ratio of the variances for undir/dir song 
%(actual data) and find the percent of draws where the ratio >= that of the
%real data

%Inputs:
%ff_alone        array of FF for undir renditions
%ff_dir          array of FF for dir renditions
%ndraws          number of times the categories are randomly assigned
%
%Outputs:
%
%ratio_alonevdir    %ratio of variance alone/variance dir (real data)
%percent_greater    % of random draws where the ratio>=that of the real data

number_greater=0;

if nargin < 3
    ndraws=input('How many random draws?');
end

ratio_alonevdir=var(ff_alone)/var(ff_dir);

[nrows, ncols]=size(ff_alone);
if nrows < ncols
    ff_alone=ff_alone';
end

[nrows, ncols]=size(ff_dir);
if nrows < ncols
    ff_dir=ff_dir';
end

total_array=[ff_alone;ff_dir];

n_alone=length(ff_alone);       %number of renditions in the undir cond
n_dir=length(ff_dir);           %number of renditions in the dir cond
n_total=n_alone+n_dir;          


for j=1:ndraws
    
    x=randperm(n_total);
    
    for i=1:n_alone
        alone(i)=total_array(x(i));
    end
    
    for m=n_alone+1:n_alone+n_dir;
        dir(m-n_alone)=total_array(x(m));
    end
    
    tempvar=var(alone)/var(dir);

    if tempvar>=ratio_alonevdir;
        number_greater=number_greater+1;
    end
        
end

percent_greater=(number_greater/ndraws)*100;
disp(strcat('% greater/equal to observer ratio of variance = ',num2str(percent_greater)))
disp(strcat('ratio of variance(alone)/variance(dir) = ', num2str(ratio_alonevdir)))

