function ptrdata = textptr
% Pointer data for the letter "a"

ptrdata = NaN*ones(16,16);

ptrdata(4,7:11) = 1;
ptrdata(5,6:12) = 1;
ptrdata(6,6) = 1;
ptrdata(6,11:13) = 1;
ptrdata(7,12:13) = 1;
ptrdata(8,7:13) = 1;
ptrdata(9,5:13) = 1;
ptrdata(10,5:6) = 1;
ptrdata(10,12:13) = 1;
ptrdata(11,5:6) = 1;
ptrdata(11,11:13) = 1;
ptrdata(12,5:13) = 1;
ptrdata(13,6:10) = 1;
ptrdata(13,12:13) = 1;

