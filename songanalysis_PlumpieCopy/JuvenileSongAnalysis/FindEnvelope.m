function [Envelope] = FindEnvelope(x,sf,y) 

threshold = 0;
p_maxima = find(sf > threshold);

maxima = find((diff(p_maxima)) > 1);

%Not finding maxima for the first resp cycle
%[maximum,resp_maxima_index(1)] = max(sf(p_maxima(1):p_maxima(maxima(1))));
%resp_maxima_index(end) = resp_maxima_index(end) + p_maxima(1) - 1;

for i = 1:(length(maxima)-2),
    [maximum(i), resp_maxima_index(i)] = max(sf(p_maxima(maxima(i) + 1):(p_maxima(maxima(i+1))))); 
    resp_maxima_index(end) = resp_maxima_index(end) + p_maxima(maxima(i) + 1) - 1;
end

%Not finding maxima for the last two resp cycles
%[maximum, resp_maxima_index(end + 1)] = max(sf(p_maxima(maxima(end) + 1)):sf(p_maxima(end)));
%resp_maxima_index(end) = resp_maxima_index(end) + p_maxima(maxima(end) + 1) - 1;

pre_zero = find(sf > (0.08 * median(maximum)));
pre_zero_actual = 1;
pre_zero_actual = [pre_zero_actual; find((diff(pre_zero)) > 1)];
pre_zero_actual(end + 1) = length(pre_zero);

post_zero = find(sf < (0.08 * median(maximum)));
post_zero_actual(1) = 1;
post_zero_actual = [post_zero_actual; find((diff(post_zero)) > 1)];
post_zero_actual(end + 1) = length(post_zero);

for i = 1:length(resp_maxima_index),
%    temp = find((sf(1:resp_maxima_index(i)) < 0),1,'last');
    temp = find((post_zero(post_zero_actual) < resp_maxima_index(i)),1,'last');
    resp_zero_crossing(i,1) = post_zero(post_zero_actual(temp));
%    temp = find((sf(resp_maxima_index(i):end) < 0),1,'first');
    temp = find((pre_zero(pre_zero_actual) > resp_maxima_index(i)),1,'first');
    resp_zero_crossing(i,2) = pre_zero(pre_zero_actual(temp));
end

resp_zero_crossing = unique(resp_zero_crossing,'rows');

resp_maxima_index = [];
for zero_crossing = 1:size(resp_zero_crossing(:,1)),
    [maximum, resp_maxima_index(zero_crossing)] = max(sf(resp_zero_crossing(zero_crossing,1):resp_zero_crossing(zero_crossing,2)));
    resp_maxima_index(end) = resp_maxima_index(end) + resp_zero_crossing(zero_crossing,1) - 1;
end

MaximaX = x(resp_maxima_index);
MaximaY = y(resp_maxima_index);

if (size(x,1) == 1)
    x = x';
end
if (size(MaximaX,1) == 1)
    MaximaX = MaximaX';
end
if (size(MaximaY,1) == 1)
    MaximaY = MaximaY';
end
Envelope = interp1q(MaximaX, MaximaY, x);

disp(['Calculated envelope']);