function G = mapAlphaBetaToGamma(alpha, beta)
% mapAlphaBetaToGamma  Map multi-indices alpha,beta (values in {1,2,3})
% to all possible gamma multi-indices under the specified per-component rule.
%
% Inputs:
%   alpha, beta : 1xd or d x 1 integer vectors with entries in {1,2,3}
%
% Output:
%   G : N x d integer matrix. Each row is one gamma multi-index.
%       N = product over i of number of possible gamma_i values.

    % Ensure row vectors
    if iscell(alpha)
        alpha = cell2mat(alpha);
    end
    if iscell(beta)
        beta = cell2mat(beta);
    end

    alpha = alpha(:).';
    beta  = beta(:).';
    d = numel(alpha);

    if numel(beta) ~= d
        error('alpha and beta must have the same length.');
    end

    % Build per-dimension option lists for gamma_i
    opts = cell(1, d);
    for i = 1:d
        opts{i} = gammaOptions(alpha(i), beta(i));
    end

    % Cartesian product across dimensions to form all gamma multi-indices
    grids = cell(1, d);
    [grids{:}] = ndgrid(opts{:});   % each grids{i} is an array of size |opts1|x...x|optsd|

    N = numel(grids{1});
    G = zeros(N, d);
    for i = 1:d
        G(:, i) = grids{i}(:);
    end
end

function gi = gammaOptions(a, b)
% gammaOptions  Return the vector of possible gamma values for one component.

    % Mapping table as cell array: row=a, col=b
    T = { ...
        1,   2,    3; ...
        2,   2,   [2 3]; ...
        3,  [2 3], 2 ...
    };

    gi = T{a, b};
end
