function p = sosgetsolmat(sosp,VAR,digits)

if nargin == 2
    digits = 5;   % Default
end;

p = VAR*0;
for i=1:size(VAR,1)
    for j=1:size(VAR,2)
        p(i,j) = sosgetsol(sosp,VAR(i,j),digits);
    end
end
