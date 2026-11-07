function [sol,M] = bl_b_stab2(frac,setname,exec,pz)
% 2-D reaction-diffusion, Dirichlet on all four edges. lam* = 2*pi^2.
%
% pz (optional): the eq_use_psatz vector for the negativity constraint. Omitted
% means the shipped default [0;0], i.e. no Positivstellensatz term. [3;4;5;6]
% selects the four LINEAR box generators added to poslpivar_2d (0cf1add1).
%
% Enabling them takes two steps that are easy to get wrong, both measured:
%   - shipped settings carry only TWO eq_opts_psatz entries, so entries 3 and 4
%     must be CLONED from an existing one; assigning just .psatz creates a
%     struct with no `exclude` field and poslpivar_2d dies on it.
%   - each psatz block needs the FULL eq_deg. The eq_deg-1 convention was
%     measured to destroy the certificate outright (reach 0.9999 -> 0.0000).
%
% CC, 11/07/2026: exec names a stock executive; the container executive of the
%   same LPI runs (bl_sop), since the stock ones reject a copvar PIE. pz sets
%   the settings_2d psatz fields, which the container path reads only for the
%   'stock' slack; its default slack (lpi_ineq_sop, psatz 'auto') already holds
%   the four linear face generators in 2-D, codes 3-6. pz is applied only when
%   the settings carry eq_opts_psatz: the 2-D settings set use_sosineq = 1 and
%   no eq_opts_psatz since the master merge of 10/27/2026 (998a68cc), where
%   the old code stopped on the missing field.
clear stateNameGenerator
pvar s1 s2 t
x = pde_var(1,[s1;s2],[0,1;0,1]);
PDE = [diff(x,t,1)==diff(x,s1,2)+diff(x,s2,2)+frac*2*pi^2*x;
       subs(x,s1,0)==0; subs(x,s1,1)==0;
       subs(x,s2,0)==0; subs(x,s2,1)==0];
PIE = initialize(convert(PDE));
st  = bl_settings(setname);
%if nargin >= 4 && ~isempty(pz)                                             % CC, 11/07/2026 (was)
if nargin >= 4 && ~isempty(pz) && isfield(st.settings_2d,'eq_opts_psatz')   % CC, 11/07/2026
    st.settings_2d.eq_use_psatz = pz(:);
    otmpl = st.settings_2d.eq_opts_psatz{1};
    for j = 1:numel(pz)
        o = otmpl;  o.psatz = pz(j);
        st.settings_2d.eq_opts_psatz{j} = o;
        st.settings_2d.eq_deg_psatz{j}  = st.settings_2d.eq_deg;
    end
elseif nargin >= 4 && ~isempty(pz)                                          % CC, 11/07/2026
    warning('bl_b_stab2:pz','pz not applied: no eq_opts_psatz');            % CC, 11/07/2026
end
%[sol,Pop,Qop] = feval(['PIETOOLS_' exec],PIE,st);                          % CC, 11/07/2026 (was)
[fn,fa] = bl_sop(exec);                                                     % CC, 11/07/2026
[sol,Pop] = feval(fn,PIE,st,fa{:});                                         % CC, 11/07/2026
Qop = [];               % the container executive has no third operator     % CC, 11/07/2026
M = struct('Pop',Pop,'Qop',Qop,'Top',PIE.T,'Aop',PIE.A);
end
