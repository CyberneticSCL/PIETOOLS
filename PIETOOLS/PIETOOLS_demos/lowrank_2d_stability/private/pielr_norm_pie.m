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
% Only a bare struct of operators is repackaged, and then only T, A and the
% vars/dom defaults -- which is all the stability builders read.

if isa(PIE,'pde_struct')
    PIE = convert(PIE,'pie');   return
end
if isa(PIE,'pie_struct')
    PIE = initialize(PIE);      return
end
if isstruct(PIE) && isfield(PIE,'type') && strcmpi(PIE.type,'pie') && isfield(PIE,'params')
    PIE = PIE.params;   return                    % a 'sys' wrapper, as lpiscript accepts
end
try
    Top = PIE.T;   Aop = PIE.A;
catch
    error('pielr_solve:input','PIE must carry T and A (opvar or opvar2d objects).');
end
is2 = isa(Top,'opvar2d') || isa(Top,'dopvar2d');
if ~(isa(Top,'opvar') || is2) || ~(isa(Aop,'opvar') || isa(Aop,'opvar2d'))
    error('pielr_solve:input','PIE.T and PIE.A must be opvar or opvar2d objects.');
end
% field assignment, not struct('T',Top,...): with an opvar2d argument that
% call dispatches to opvar2d's own struct method and errors
S = struct();
S.T = Top;   S.A = Aop;
try, S.vars = PIE.vars; catch, S.vars = []; end   %#ok<*CTCH>
try, S.dom  = PIE.dom;  catch, S.dom  = []; end
if isempty(S.vars), S.vars = [Top.var1, Top.var2]; end
if isempty(S.dom),  S.dom  = Top.I; end
S.dim = 1 + is2;
PIE = S;
end
