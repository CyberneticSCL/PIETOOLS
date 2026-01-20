function C = plus(A,B)
C = B;
% need error handling: add checks to ensure A and B are compatible
% dimensions of params in A and B, vars_in and vars_out

for i=1:size(A.params,1)
    for j=1:size(A.params,2)
        C.params{i,j} = A.params{i,j}+C.params{i,j};
    end
end
end
