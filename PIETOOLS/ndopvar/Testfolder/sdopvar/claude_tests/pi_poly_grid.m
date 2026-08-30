function V = pi_poly_grid(P,vars,pts)
% V = PI_POLY_GRID(P,VARS,PTS) evaluates a matrix-valued 'polynomial' at a
% set of sample points, returning the values as a matrix.
%
% INPUTS
% - P:      matrix-valued 'polynomial' object in a subset of VARS;
% - vars:   1 x n 'cellstr' object naming the variables, indexing the
%           columns of PTS;
% - pts:    G x n array of sample points;
%
% OUTPUTS
% - V:      G x numel(P) array, with V(g,:) the entries of P evaluated at
%           PTS(g,:), in column-major order;
%
% Evaluation is vectorized over the sample points, so this is much cheaper
% than repeated calls to 'subs' when a whole grid is needed.
%
% MP, 08/22/2026: Initial coding

P = polynomial(P);
G = size(pts,1);
T = size(P.degmat,1);

pvars = P.varname;
if isempty(pvars)
    V = repmat(reshape(full(sum(P.coefficient,1)),1,[]),G,1);
    return
end

[tf,loc] = ismember(pvars,vars);
if ~all(tf)
    error("Variable '"+string(pvars{find(~tf,1)})+"' is not among the sample variables.")
end

M = ones(G,T);
for v = 1:numel(pvars)
    % Note the parentheses: transpose binds tighter than power, so without
    % them this would be read as (pts.^degmat).' and fail on the dimensions.
    M = M .* (pts(:,loc(v)).^(full(P.degmat(:,v)).'));
end

V = M*full(P.coefficient);

end
