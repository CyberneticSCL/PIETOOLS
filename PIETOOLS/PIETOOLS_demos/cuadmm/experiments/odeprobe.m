% odeprobe.m -- does an ODE-coupled plant produce a nonzero MULTIPLIER cell?
%
% Variable coefficients provably do not (measured: A.R22{1,1} = 0 for six
% variants), because the PIE fundamental state is the highest derivative and
% every coefficient multiplies an already-integrated quantity.  A multiplier
% survives only where a state component carries ZERO spatial derivative order
% -- i.e. a finite-dimensional ODE state.  That also moves Top.dim(:,1) off
% [0;0;0;1], which is the regime where the degree-balance rule is separately
% known to break (multiplicity: single-factor blocks need offset e, not e/2).
cuadmm_path;
pvar s1 s2 t

fprintf('OP variant|Top.dim|T.R00|A.R00|Q.R00|Q nonzero cells (R00/R22)\n');
V = {'ODE coupled via BC', 1; 'ODE coupled via source', 2};
for k = 1:size(V,1)
    try
        [Top,Aop] = mkode(V{k,2});
        np  = Top.dim(:,1);
        Iop = opvar2d(eye(sum(np)),[np np],Top.I,[Top.var1 Top.var2]);
        Qop = clean_opvar((Aop'*(Iop*Top))' + Aop'*(Iop*Top),1e-12);
        nzc = '';
        if nz(Qop.R00), nzc = [nzc 'R00 ']; end
        for i=1:3, for j=1:3
            if nz(Qop.R22{i,j}), nzc = [nzc sprintf('R22{%d%d} ',i,j)]; end
        end, end
        fprintf('OP %s|[%s]|%d|%d|%d|%s\n', V{k,1}, num2str(np'), ...
                nz(Top.R00), nz(Aop.R00), nz(Qop.R00), strtrim(nzc));
    catch ME
        fprintf('OP %s|FAILED|%s\n', V{k,1}, strrep(ME.message,newline,' '));
    end
end
fprintf('OPDONE\n');

function [Top,Aop] = mkode(kind)
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0 1;0 1]);
X = pde_var('state',1,[],[]);
switch kind
    case 1   % ODE enters through a boundary condition (as in DEMO1)
        sys = [diff(x,t,1) == diff(x,s1,2)+diff(x,s2,2)+2*x;
               diff(X,t,1) == -2*X;
               subs(x,s1,0)==0;  subs(x,s1,1)==X;
               subs(x,s2,0)==0;  subs(x,s2,1)==0];
    case 2   % ODE driven by the PDE, PDE driven by the ODE
        sys = [diff(x,t,1) == diff(x,s1,2)+diff(x,s2,2)+2*x + X;
               diff(X,t,1) == -2*X + int(int(x,s1,[0,1]),s2,[0,1]);
               subs(x,s1,0)==0;  subs(x,s1,1)==0;
               subs(x,s2,0)==0;  subs(x,s2,1)==0];
end
PIE = initialize(convert(sys));  Top = PIE.T; Aop = PIE.A;
end

function b = nz(p)
if isempty(p), b = 0; return; end
if isa(p,'double'), b = double(any(p(:)~=0)); return; end
b = double(~isempty(p.coefficient) && any(any(p.coefficient~=0)));
end
