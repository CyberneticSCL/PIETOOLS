function varargout = spIkron(q,varargin)
% IC = SPIKRON(Q,C) takes the Kronecker product of the identity matrix of
% dimension QxQ and the sparse matrix defined by C
%
% INPUTS
% - q:  scalar integer specifying the dimension of the identity matrix
% - C:  Either:
%       1. a sparse matrix, or
%       2. the elements C{1},...,C{5} (and optionally C{6}) of a sparse
%          matrix, passed as separate arguments, with
%           Cmat = sparse(C{1}, C{2}, C{3}, C{4}, C{5}),
%       and elements:
%       C{1} - m x 1 array of nonnegative integers, specifying the row
%               indices in Cmat
%       C{2} - m x 1 array of nonnegative integers, specifying the column
%               indices in Cmat
%       C{3} - m x 1 array specifying the values of the nonzero elements in
%               rows C{1} and columns C{2} of Cmat
%       C{4} - scalar integer specifying the row dimension of Cmat
%       C{5} - scalar integer specifying the column dimension of Cmat
%       C{6} - (optional) number of nonzero elements in Cmat; read for
%               IC{6} only (IC{6} = q*C{6}; without it, nnz(IC{3}))
%
% OUTPUTS
% - IC: with one output, the sparse matrix I o Cmat; with 5 or 6 outputs,
%       the elements IC{1},...,IC{6} with
%           ICmat = sparse(IC{1}, IC{2}, IC{3}, IC{4}, IC{5}).
%
% Initial coding MP, SS, DJ, AT, DB - 10/24/2026
% MMP, 11/06/2026: Comment pass (maintainer, 11/06/2026): INPUTS/OUTPUTS
%                  corrected to the argument form the code reads (separate
%                  arguments, not one cell); code unchanged. The branch
%                  'if varargout==6' compares a cell with 6 and never runs
%                  (nargout==6 is meant); left as is.
% MMP, 11/08/2026: 'if varargout==6' compared a cell with 6, so every call
%                  with 2 to 6 outputs errored; replaced by 'if nargout==6'
%                  (the 11/06 note that it was left as is no longer applies).
% MMP, 11/08/2026: The one-output form no longer passes C{6} as the nzmax of
%                  I o Cmat: sparse requires nzmax at least the number of
%                  triplets, q times those of Cmat (duplicates included),
%                  so the call errored for q > 1. C{6} now scales into IC{6}
%                  only. Help of C{6} to match.

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
%   if nargin==7                                                            % MMP, 11/08/2026 (was)
%       varargout{1} = sparse(r_idcs_I(:), c_idcs_I(:), Cvals_I(:), m*q, n*q,varargin{6}); % MMP, 11/08/2026 (was)
%   else                                                                    % MMP, 11/08/2026 (was)
%       varargout{1} = sparse(r_idcs_I(:), c_idcs_I(:), Cvals_I(:), m*q, n*q); % MMP, 11/08/2026 (was)
%   end                                                                     % MMP, 11/08/2026 (was)
    % No nzmax from C{6}: sparse needs nzmax at least the number of         % MMP, 11/08/2026
    % triplets, q times those of Cmat (duplicates included), so the count   % MMP, 11/08/2026
    % of Cmat errored for q > 1; sparse allocates the exact size itself.    % MMP, 11/08/2026
    varargout{1} = sparse(r_idcs_I(:), c_idcs_I(:), Cvals_I(:), m*q, n*q);  % MMP, 11/08/2026
elseif nargout<=6
    varargout{1} = r_idcs_I;
    varargout{2} = c_idcs_I;
    varargout{3} = Cvals_I;
    varargout{4} = m*q;
    varargout{5} = n*q;
%   if varargout==6                                                         % MMP, 11/08/2026 (was)
    if nargout==6                                                           % MMP, 11/08/2026
        if nargin==7
            varargout{6} = varargin{6}*q;
        else
            varargout{6} = nnz(Cvals_I);
        end
    end
end

end