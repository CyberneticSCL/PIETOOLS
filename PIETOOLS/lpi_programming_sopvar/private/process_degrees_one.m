function d = process_degrees_one(deg,n3)
% Expand a single degree specification.

if isnumeric(deg) && isscalar(deg)
    deg = struct('int',deg);
elseif ~isa(deg,'struct')
    error("Monomial degrees should be specified as a scalar or 'struct' object.")
end

if ~isfield(deg,'int') || isempty(deg.int)
    deg.int = 1;
end
if ~isfield(deg,'mult') || isempty(deg.mult)
    deg.mult = deg.int;
end

d = struct();
d.int = expand_caps(deg.int,n3,'int');
d.mult = expand_caps(deg.mult,n3,'mult');

if ~isfield(deg,'joint') || isempty(deg.joint)
    d.joint = sum(d.int)+sum(d.mult);
elseif ~isscalar(deg.joint)
    error("The joint degree should be specified as a scalar.")
else
    d.joint = deg.joint;
end

% Optional cap per variable SUBSET, over the order the basis is built in,    % MMP, 09/21/2026
% [theta_1..theta_n3, s_1..s_n3], indexed as 'build_exponent_grid'           % MMP, 09/21/2026
% documents. 'int', 'mult' and 'joint' are the singleton and full-set        % MMP, 09/21/2026
% special cases, so this subsumes them; empty leaves the basis determined    % MMP, 09/21/2026
% by those three alone, unchanged from before. Needed because a             % MMP, 09/21/2026
% 'poslpivar_2d' basis is cut out by 2^nvars subset caps, which per-variable % MMP, 09/21/2026
% plus total cannot express -- see 'settings2possopvar'.                     % MMP, 09/21/2026
if ~isfield(deg,'subset') || isempty(deg.subset)                            % MMP, 09/21/2026
    d.subset = [];                                                          % MMP, 09/21/2026
elseif numel(deg.subset)~=2^(2*n3)                                          % MMP, 09/21/2026
    error("A 'subset' degree array should have 2^(2*n3) entries for n3 "...
          +"spatial variables.")                                            % MMP, 09/21/2026
else                                                                        % MMP, 09/21/2026
    d.subset = reshape(deg.subset,1,[]);                                    % MMP, 09/21/2026
end                                                                         % MMP, 09/21/2026

end


