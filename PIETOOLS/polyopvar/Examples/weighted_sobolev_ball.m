function g = weighted_sobolev_ball(r, alpha, Top, v)

% g = weighted_sobolev_ball(...) returns the weighted sobolev ball (Sec. 7.1 in Automatica paper).
%
% INPUTS
% - r:         Radius of the local ball.
% - alpha:     (n+1)-dim array containing parameters of weighted Sobolev ball.
% - Top:       Map from fundamental to PDE state represented as an ndopvar object.
% - v:         Fundamental state.

%
% OUTPUTS
% - g:         Weighted sobolev ball of radius r.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n         = size(alpha) - 1; % degree of PDE.
Top_opvar = ndopvar2dopvar(Top);


g = r^2 - alpha(1)*innerprod(Top*v,Top*v);

for i=2:n+1
    Rop = dopvar2ndopvar(diff(Top_opvar,Top_opvar.var1,i-1,'pure'));
    g = g - alpha(i)*innerprod(Rop*v,Rop*v);
end