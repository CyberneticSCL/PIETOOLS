function X = cx_on_registry(X,R)
% X = CX_ON_REGISTRY(X,R) restates container X over the variable registry
% of container R (a superset of X's): same blocks, same operator, metadata
% remapped. Bridges the one real gap for executives: @copvar/@cdopvar plus
% and mtimes demand IDENTICAL registries, so an R^n -> R^n operator (empty
% registry, e.g. -gam*Iw or Dzu) cannot meet an L_2 one in a sum or product
% until restated. Same as on_registry in test_copvar_blockops.
%
% Initial coding MMP, 09/25/2026
% MMP, 09/26/2026: Redundant since 09/26/2026: plus and mtimes merge
%                  registries themselves (merge_copvar_registry); calls are kept,
%                  and return at once when the registries already agree.
if isequal(X.vars,R.vars) && isequal(X.dom,R.dom),  return,     end
meta = metadata(X);
[tf,loc] = ismember(X.vars,R.vars);
if ~all(tf)
    error('cx_on_registry:notSubset','X''s registry is not contained in R''s.')
end
so = false(size(X.C,1),numel(R.vars));  so(:,loc) = X.space_out;
si = false(size(X.C,2),numel(R.vars));  si(:,loc) = X.space_in;
meta.vars = R.vars;     meta.dom = R.dom;
meta.space_out = so;    meta.space_in = si;
if isa(X,'cdopvar'),    X = cdopvar(X.C,meta);
else,                   X = copvar(X.C,meta);
end
end
