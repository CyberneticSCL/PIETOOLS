function [prog,P] = sosmatrvar(prog,Z,sizePQ,wscoeff)
%[prog,P]=SOSMATRVAR(prog,Z,sizePQ,wscoeff) declares a n x m matrix-valued
%polynomial in monomials Z. 
%
% Inputs: 
% prog - The SOSprogram to which to attach the variable
% Z - The list of monomials which appear in the elements of the matrix.
% sizePQ - The dimension of the resulting matrix will be sizePQ so that
%               sizePQ(1)=#rows, sizePQ(2)=#cols;
% wscoeff - if wscoeff='wscoeff', elements of the matrix will be declared
% as pvars in the workspace
%
% NOTES:
% Distributed with DelayTOOLS
% Compatable with MULTIPOLY and SOSTOOLS as of June 2013
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% Version control:
% 7_26_2019 - Added additional support for matrices of zero dimension.
%
% NOTES:
% Needs to be improved, but still not worse that sosmatrixpoly, which is a
% blatent copy



sizeP = sizePQ(1);
sizeQ = sizePQ(2);

if sizeP==0 || sizeQ==0
    P=polynomial(zeros(sizePQ));
else
    
    if nargin > 3 & wscoeff == 'wscoeff'
        for i = 1:sizeP*sizeQ
            [prog,matrPi(i)] = sospolyvar(prog,Z,'wscoeff');
        end
    else;
        for i = 1:sizeP*sizeQ
            [prog,matrPi(i)] = sospolyvar(prog,Z);
        end
    end
    for i = 1:sizeP;
        for j = 1:sizeQ;
            P(i,j) = matrPi((i-1)*sizeQ+j);
        end
    end
end