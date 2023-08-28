clc; clear;

pvar s t theta;
x = state('pde'); z = state('out'); w = state('in');

% heat equation
eqns = [diff(x,t)==8*x+diff(x,s,2); 
             z==int(x,s,[0,1]);
             subs(x,s,0)==0; subs(x,s,1)==0];
% % transport equation
% eqns = [diff(x,t)==diff(x,s); 
%              z==int(x,s,[0,1]);
%              subs(x,s,1)==w;];
% % reaction-diffusion equation
% eqns = [diff(x,t)==2*x+diff(x,s,2)+s*w; 
%              z==int(x,s,[0,1]);
%              subs(x,s,0)==0; subs(diff(x,s),s,1)==0];
%%
pde = sys();
pde = addequation(pde,eqns);
pie = convert(pde);
PIE = pie.params;
%%
st = lpisettings('heavy');
st.sos_opts.solver = 'mosek';
st.epneg = 0.02;
% [~,P] = lpisolve(PIE,st,'stability');
[~,~,decayrate] = PIETOOLS_exponential_stability(PIE,st);
[~,~,decayrate_dual] = PIETOOLS_exponential_stability_dual(PIE,st);
