function [xin, pts, Y, nb] = p0_apply(X, q, P, seed)
% [XIN,PTS,Y,NB] = P0_APPLY(X,Q,P,SEED) applies the fixed container X to Q
% random polynomial test functions per input component, by the definition
% of each sopvar block (apply_sopvar), and evaluates each output component
% at P random points of its domain.
%
% INPUT
% - X:      'copvar' without decision variables;
% - q:      double, number of test functions (columns);
% - P:      double, number of evaluation points per output component;
% - seed:   double, seed of the random stream (test functions and points).
%
% OUTPUT
% - xin:    numel(X.dim_in) x 1 cell; xin{j} is a dim_in(j) x q 'polynomial'
%           in the variables of input component j of total degree <= 2, or
%           a 'double' for a component without variables;
% - pts:    numel(X.dim_out) x 1 cell; pts{i} is P x (number of variables of
%           output component i), the evaluation points, in the order of
%           X.vars;
% - Y:      numel(X.dim_out) x 1 cell; Y{i} is dim_out(i) x q x P, the output
%           component i at the points;
% - nb:     double, the number of nonempty blocks applied.
rng(seed,'twister');
vars = X.vars;  dom = X.dom;
N = numel(X.dim_in);  M = numel(X.dim_out);
xin = cell(N,1);
for j = 1:N
    xin{j} = p0_randpoly(X.dim_in(j), q, vars(X.space_in(j,:)));
end
Y = cell(M,1);  pts = cell(M,1);  nb = 0;
for i = 1:M
    acc = polynomial(zeros(X.dim_out(i), q));
    for j = 1:N
        B = X.C{i,j};
        if isempty(B),  continue,   end
        acc = acc + apply_sopvar(B, xin{j});
        nb = nb + 1;
    end
    vi = vars(X.space_out(i,:));
    di = dom(X.space_out(i,:),:);
    pts{i} = di(:,1)' + rand(P, numel(vi)).*(di(:,2)-di(:,1))';
    Y{i} = p0_evalpoly(acc, vi, pts{i});
end
end
