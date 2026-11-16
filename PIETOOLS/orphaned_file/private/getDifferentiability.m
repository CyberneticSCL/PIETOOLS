function n = getDifferentiability(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% not used currently

if nargin==1
    nBC = varargin{1};
elseif nargin==2
    nBC = varargin{1};
    nTotal = varargin{2};
end

nmax = cell(1,floor(sqrt(nBC)));

for i = 1:floor(sqrt(nBC)) % there need to be atleast i BCs for i times differentiable state
    nmax{i} = 0:floor(nBC/i);
end

% generate all permutations for n using nmax
nGrid = cell(1,floor(sqrt(nBC)));
[nGrid{:}] = ndgrid(nmax{:});
nGrid = cellfun(@(x) x(:), nGrid,'uniformoutput',false);
nGrid = [nGrid{:}];
nGrid = nGrid.*(1:floor(sqrt(nBC)));

if exist('nTotal')
    n_possible = nGrid((sum(nGrid,2)==nTotal),:);
end
n = n_possible;
end