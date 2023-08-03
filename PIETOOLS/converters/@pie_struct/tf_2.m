function varargout = tf_2(PIE,w)
% assuming inputs and outputs are finite dimensional
% pde doesnt have ode coupling

nx = PIE.T.dim(:,2);
nz = PIE.C1.dim(1,1);
nw = PIE.B1.dim(1,2);

C = PIE.C1; B = PIE.B1; Tw = PIE.Tw;
T = PIE.T; A = PIE.A;

opvar Z Zb;
Z.dim = PIE.C1.dim;
Zb.dim = PIE.B1.dim;

bigC = [C,Z;Z,C]; 
bigT = [-A,w*T; -w*T,-A];
bigB = [B-w*Tw,Zb;Zb,B-w*Tw];




bigTF = @(w) bigC*inv(bigT)*bigB;

tmp = squeeze(trapz(dS, TF)); % TF is size n_grid*2*nz*2*nw;

varargout{1} = tmp(1:nz,1:nw);
varargout{2} = tmp(1:nz,nw+1:2*nw);
end