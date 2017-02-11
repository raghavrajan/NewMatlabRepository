function zs = interpzeros(x)

zs = [];

for i = 1:length(x)-1
    %if there is a zero crossing on this interval
    if ((x(i) > 0) - (x(i+1) > 0)) ~= 0 
        dy = x(i+1) - x(i);
        b = x(i) - dy*i;
        z = -b/dy;
        zs = [zs;z];
    end
end