function PIE = pielr_norm_pie(PIE)                                          % CC, 09/23/2026
% PIELR_NORM_PIE  Normalise the PIE input for pielr_solve.
%
% Generalises pielr_certify's local norm_pie, which requires opvar2d and
% errors "this package is 2-D only".  Here 1-D and 2-D are both allowed, and
% a pde_struct is converted, so the entry point takes the same range of input
% lpiscript does.
%
% A pie_struct is returned untouched: the l2gain builder needs Bw, Cz, Dzw and
% Tw as well as T and A, and initialize() on a pie_struct is what fills them.
% Any other input is an error.
%
% MMP, 11/26/2026: the repackaging of a bare struct of opvar or opvar2d
%   operators deleted, with its help (see the end of the code).

if isa(PIE,'pde_struct')
    PIE = convert(PIE,'pie');   return
end
if isa(PIE,'pie_struct')
    PIE = initialize(PIE);      return
end
if isstruct(PIE) && isfield(PIE,'type') && strcmpi(PIE.type,'pie') && isfield(PIE,'params')
    PIE = PIE.params;   return                    % a 'sys' wrapper, as lpiscript accepts
end
% BEGIN MMP, 11/26/2026: deleted, CC 09/23/2026: the repackaging of a bare
% struct whose T and A are opvar or opvar2d objects into T, A, vars, dom and
% dim. Those classes left the path (f83486c8); once opvar is renamed opvar,
% a bare struct of 1-D containers would pass its class test, on a path never
% exercised.
error('pielr_solve:input',['PIE must be a ''pde_struct'', a ''pie_struct'', or a ' ...
      'struct with type ''pie'' and the field params; got %s.'],class(PIE));
% END MMP, 11/26/2026
end
