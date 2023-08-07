function varargout = tf(PIE,w)
nx = PIE.T.dim(:,2);
nz = PIE.C1.dim(1,1);
nw = PIE.B1.dim(1,2);

C = PIE.C1; B = PIE.B1; T = PIE.T; A = PIE.A;

opvar Z Zb;
Z.dim = PIE.C1.dim;
Zb.dim = PIE.B1.dim;

bigTF = @(w) [C,Z;Z,C]*inv([-A,w*T; -w*T -A],1e-12)*[B-w*PIE.Tw,Zb;Zb,B-w*PIE.Tw];

tmp = bigTF(w);

varargout{1} = tmp(1:nz,1:nw);
varargout{2} = tmp(1:nz,nw+1:2*nw);
end