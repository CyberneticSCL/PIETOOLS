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
clear stateNameGenerator
pvar s1 s2 t
x = pde_var(1,[s1;s2],[0,1;0,1]);
PDE = [diff(x,t,1)==diff(x,s1,2)+diff(x,s2,2)+frac*2*pi^2*x;
       subs(x,s1,0)==0; subs(x,s1,1)==0;
       subs(x,s2,0)==0; subs(x,s2,1)==0];
PIE = initialize(convert(PDE));
st  = bl_settings(setname);
if nargin >= 4 && ~isempty(pz)
    st.settings_2d.eq_use_psatz = pz(:);
    otmpl = st.settings_2d.eq_opts_psatz{1};
    for j = 1:numel(pz)
        o = otmpl;  o.psatz = pz(j);
        st.settings_2d.eq_opts_psatz{j} = o;
        st.settings_2d.eq_deg_psatz{j}  = st.settings_2d.eq_deg;
    end
end
[sol,Pop,Qop] = feval(['PIETOOLS_' exec],PIE,st);
M = struct('Pop',Pop,'Qop',Qop,'Top',PIE.T,'Aop',PIE.A);
end
