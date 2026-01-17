function varargout = spIkron(q,varargin)
% IC = SPIKRON(Q,C) takes the Kronecker product of the identity matrix of
% dimension QxQ and the sparse matrix defined by C
%
% INPUTS
% - q:  scalar integer specifying the dimension of the identity matrix
% - C:  Either:
%       1. A sparse matrix for which to compute
%       1x5 cell representing the sparse matrix as
%           Cmat = sparse(C{1}, C{2}, C{3}, C{4}, C{5}),
%       with elements:
%       C{1} - m x 1 array of nonnegative integers, specifying the row
%               indices in Cmat
%       C{2} - m x 1 array of nonnegative integers, specifying the column
%               indices in Cmat
%       C{3} - m x 1 array specifying the values of the nonzero elements in
%               rows C{1} and columns C{2} of Cmat
%       C{4} - scalar integer specifying the row dimension of Cmat
%       C{5} - scalar integer specifying the column dimension of Cmat
%       C{6} - (optional) number of nonzero elements in Cmat
%
% OUTPUTS
% - IC: 1x5 cell representing the sparse matrix I o Cmat as
%           ICmat = sparse(IC{1}, IC{2}, IC{3}, IC{4}, IC{5}).
%

% If C is specified as an actual array, we compute the kronecker product
% using the built-in Matlab routine
if nargin==2
    C = varargin{1};
    if ~isa(C,'double')
        error("For single input case, matrix must be specified as sparse array.")
    end
    if nargout>1
        error("Too many outputs for single input case.")
    end
    IC = sparse(kron(speye(q),C));
    varargout{1} = IC;    
    return
end

% Otherwise, assume the input C:=varargin represents a sparse matrix as
%   Cmat = sparse(C{1}, C{2}, C{3}, C{4}, C{5}).
if nargin<6
    error("Insufficient input arguments.")
elseif nargin>7
    error("The sparse matrix can be represented by a cell of at most 6 elements.")
end
% Establish nonzero elements of Cmat, and the row and column indices where
% they appear.
r_idcs = varargin{1}(:);   
c_idcs = varargin{2}(:);
Cvals = varargin{3}(:);
m = varargin{4};           
n = varargin{5};
if numel(r_idcs)~=numel(c_idcs)
    error("Number of row and column indices should match.")
end
if isscalar(Cvals)
    Cvals = Cvals*ones(size(r_idcs));
elseif numel(r_idcs)~=numel(Cvals)
    error("Number of values should match the number of row and column indices.")
end

% Determine the row and column indices of nonzero elements of I_{q} o Cmat, 
% and construct the cell representing the associated sparse matrix.
r_idcs_I = r_idcs + (0:q-1)*m;
c_idcs_I = c_idcs + (0:q-1)*n;
Cvals_I = repmat(Cvals,[q,1]);
if nargout==1
    if nargin==7
        varargout{1} = sparse(r_idcs_I(:), c_idcs_I(:), Cvals_I(:), m*q, n*q,varargin{6});
    else
        varargout{1} = sparse(r_idcs_I(:), c_idcs_I(:), Cvals_I(:), m*q, n*q);
    end
elseif nargout<=6
    varargout{1} = r_idcs_I;
    varargout{2} = c_idcs_I;
    varargout{3} = Cvals_I;
    varargout{4} = m*q;
    varargout{5} = n*q;
    if varargout==6
        if nargin==7
            varargout{6} = varargin{6}*q;
        else
            varargout{6} = nnz(Cvals_I);
        end
    end
end

end