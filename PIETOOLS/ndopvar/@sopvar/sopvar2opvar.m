function obj = sopvar2opvar(objSopvar)
% this converts 1D sopvar to 1D opvar objects
pvar s1 t1 s1_dum;

opvar obj;
obj.dim = [0,0; fliplr(objSopvar.dims)];
obj.I = objSopvar.dom_3;
R0 = quadPoly.quadPoly2polynomial(objSopvar.params{1});
R1 = quadPoly.quadPoly2polynomial(objSopvar.params{2});
R2 = quadPoly.quadPoly2polynomial(objSopvar.params{3});

obj.R = struct('R0',subs(R0,t1,s1_dum), ...
                 'R1',subs(R1,t1,s1_dum), ... 
                  'R2',subs(R2,t1,s1_dum));
end