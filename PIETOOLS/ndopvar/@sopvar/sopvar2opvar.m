function obj = sopvar2opvar(objSopvar)
% this converts 1D sopvar to 1D opvar objects
pvar s1 t1 s1_dum;

opvar obj;
obj.I = objSopvar.dom.in;
obj.dim = [0, 0; objSopvar.dims];
dims = objSopvar.dims;
R0 = quadPoly(objSopvar.params{1}, objSopvar.ZL, objSopvar.ZR,  dims, {'s1'}, {'s1'}, 0)   ;
R1 = quadPoly(objSopvar.params{2}, objSopvar.ZL, objSopvar.ZR,  dims, {'s1'}, {'t1'}, 0)   ;
R2 = quadPoly(objSopvar.params{3}, objSopvar.ZL, objSopvar.ZR,  dims, {'s1'}, {'t1'}, 0)   ;

R0 = combine(R0);
R1 = combine(R1);
R2 = combine(R2);

obj.R.R0 = quadPoly.quadPoly2polynomial(R0);
obj.R.R1 = quadPoly.quadPoly2polynomial(R1);
obj.R.R2 = quadPoly.quadPoly2polynomial(R2); 

% obj.R = struct('R0',subs(R0,t1,s1_dum), ...
%                  'R1',subs(R1,t1,s1_dum), ... 
%                   'R2',subs(R2,t1,s1_dum));
end