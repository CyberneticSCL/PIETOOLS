function D = sdvar2dpvar(P,dims,vars,ZL,ZR,Zd)
% D = SDVAR2DPVAR(P,DIMS,VARS,ZL,ZR,ZD) takes an 'sdvar' object P,
% representing a decision variable matrix in terms of the 'sdopvar' format,
% and returns a 'dvpar' object D representing the same decision variable,
% so that
%   D = (Im o ZL{1}(s1) o ... o ZL{M}(sM))^T unvec(P.A + P.B'*d)
%           (In o ZR{1}(t1) o ... o ZR{N}(tN))
%
% INPUTS
% - P:  struct with fields 'A' and 'B', specifying the coefficients
%       defining the decision variable matrix in the 'sdopvar' form
% - dims:   1 x 2 array [m,n], specifying the row and column dimensions of
%           the decision variable matrix D;
% - vars:   struct with fields 'in' and 'out', each a row cell array
%           specifying the names of the input and output independent 
%           variables in terms of which the object is defined;
% - ZL:     1 x M cell (for M output variables), specifying the degrees of
%           the monomials in each output variable;
% - ZR:     1 x N cell (for N input variables), specifying the degrees of
%           the monomials in each input variable;
% - Zd:     nZd x 1 cellstr specifying the names of the decision variables,
%           'd', in terms of which the object is defined.
%
% OUTPUTS
% - D:  m x n 'dpvar' object representing the decision variable matrix
%       defined by coefficients P.A and P.B.
%


% Check that the 'sdvar' object is properly specified
if ~isstruct(P) || ~isfield(P,'A') || ~isfield(P,'B')
    error("'sdvar' object must be specified as 'struct' with fields 'A' and 'B'.")
end
A = P.A;
B = P.B;
if size(A,2)~=1
    error("A should be specified as column vector.")
elseif size(A,1)~=size(B,2)
    error("Number of columns of B should match number of rows of A.")
elseif size(B,1)~=numel(Zd)
    error("Row dimension of B must match the number of decision variables.")
end
ncoeffs = size(A,1);

% Extract row and column dimensions
m = dims(1);
n = dims(2);

% Extract variable names
vars1 = vars.out;       M = numel(vars1);
vars2 = vars.in;        N = numel(vars2);

% Determine the number of monomials in each input and output variable
if numel(ZL)~=M
    error('Number of left-monomials should mathc number of output variables.')
elseif numel(ZR)~=N
    N = numel(ZR);
end
nZL_arr = cellfun(@(a)numel(a),ZL);
nZR_arr = cellfun(@(a)numel(a),ZR);
nnZL_arr = [1,cumprod(nZL_arr)];
nnZR_arr = [1,cumprod(nZR_arr)];
nZL = nnZL_arr(end);
nZR = nnZR_arr(end);
if ncoeffs~=m*n*nZL*nZR
    error("Number of coefficients should equal to the total number of elements times the total number of monomials.")
end

% Determine the full list of independent and decision variables
Zd = Zd(:);
nZd = numel(Zd);
varname = [vars1(:); vars2(:)];

% Construct the full degmat associated with the Kronecker product of the
% monomials,
%   ZL{1} o ZL{2} o ... o ZL{M} o ZR{1} o ... o ZR{N}
degmat = zeros(1,0);
for i=1:M
    degmat = [repelem(degmat,nZL_arr(i),1),repmat(ZL{i},nnZL_arr(i),1)];
end
for i=1:N
    degmat = [repelem(degmat,nZR_arr(i),1),repmat(ZR{i},nZL*nnZR_arr(i),1)];
end

% For each row of A and B', determine which row and column in the 
% matrix-valued object it corresponds to
ridcs1 = (1:m*nZL)';        cidcs1 = (1:n*nZR)';
ridcs = ceil(ridcs1/nZL);   cidcs = ceil(cidcs1/nZR);
% Also determine which monomial in ZL and ZR it corresponds to
ZLidcs = ridcs1 - (ridcs-1)*nZL;
ZRidcs = cidcs1 - (cidcs-1)*nZR;
% Account matrix structure
ridcs = repmat(ridcs,n*nZR,1);
cidcs = repelem(cidcs,m*nZL,1);
Zidcs = (repmat(ZLidcs,n*nZR,1)-1)*nZR + repelem(ZRidcs,m*nZL,1);

% Extract nonzero coefficients from A and B
[ridcsA,~,valsA] = find(A);
[didcsB,ridcsB,valsB] = find(B);        % note that we store B, not B^T
% Determine row and column numbers of coefficients associated with constant
% term
ridcsC1 = (ridcs(ridcsA)-1)*nZd + 1;    % add 1 for constant term
cidcsC1 = (cidcs(ridcsA)-1)*nZL*nZR + Zidcs(ridcsA);
% Determine row and column numbers of coefficients associated with decision
% variables
ridcsC2 = (ridcs(ridcsB)-1)*nZd + 1 + didcsB;    % add 1 for constant term
cidcsC2 = (cidcs(ridcsB)-1)*nZL*nZR + Zidcs(ridcsB);
% Declare the sparse coefficient matrix defining the dpvar object
ridcsC = [ridcsC1(:);ridcsC2(:)];
cidcsC = [cidcsC1(:);cidcsC2(:)];
valsC = [valsA(:); valsB(:)];
C = sparse(ridcsC,cidcsC,valsC,m*(nZd+1),n*nZL*nZR);

% Construct the dpvar object
D = dpvar(C,degmat,varname,Zd,[m,n]);

end