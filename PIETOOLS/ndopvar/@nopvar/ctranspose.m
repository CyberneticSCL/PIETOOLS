function Pt = ctranspose(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pt = ctranspose(P) returns the conjuagate transpose of 'nopvar' object P.
% Note that current implementation returns NaNs for dim due to transposing
% all elements in cell.
%
% Version: 1.0
% 
% INPUT
% P:  nopvar object;
% 
% OUTPUT
% Pt: nopvar object representing adjoint of P;
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding CR - 1/16/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isa(P,'nopvar')
    error('Input must be an nopvar variable.');
end

% test if \alpha \in {1,2}^N.
N = ndims(P.C); % number of leading 3-sized dimensions.
idx = cell(1, N);
[idx{:}] = ndgrid(1:3);
combos = cellfun(@(x) x(:), idx, 'UniformOutput', false);
combos = [combos{:}]; % size = (3^N) × N array containing all index combinations.
mask = all(combos~=1,2);
combos = combos(mask, :); % contains all indices of N leading dimensions which do not have a 1.

temp = P.C;
linIdx = sub2ind(size(temp), combos(:,1), combos(:,2));
temp(linIdx) = {[]}; % set all permitted non-zero elements to empty.
allZeroOrEmpty = all(cellfun(@(x) isempty(x) || (isnumeric(x) && nnz(x) == 0), temp(:)));

if ~allZeroOrEmpty
    error('Input must contain no delta terms!');
end

% Instantiate adjoint operator.
Pt = nopvar();
Pt.deg = P.deg;
Pt.dom = P.dom;
Pt.vars = P.vars;
Pt.C = P.C;

% transpose all matrices - this leads to incorrect sizes for cells corresponding to k=0 (see nopvar class).
for ii=1:3^N
    Pt.C{ii} = Pt.C{ii}.';
end

% swap 2 and 3 in first N dimensions.
idx = repmat({[1 3 2]}, 1, N);
Pt.C = Pt.C(idx{:});

end