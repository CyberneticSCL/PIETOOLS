function logval = poly2double(in)
% checks if a given input in can be converted to double and returns 1 if
% successful
try
    tmp = double(in);
    logval = 1;
catch
    logval = 0;
end
end