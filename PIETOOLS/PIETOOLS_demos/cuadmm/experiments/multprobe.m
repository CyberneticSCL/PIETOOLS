% multprobe.m -- FIND a plant whose POINTWISE MULTIPLIER cell is nonzero.
%
% Every plant measured so far has Qop.R22{1,1} = 0 exactly, because Top is a
% Green operator with no multiplier part and Qop = A'T + T'A composes each
% operator WITH T -- and multiplier-then-integral is still an integral.  So a
% variable coefficient may not be enough.  Probe candidates and report which
% cells are nonzero in Top, Aop and Qop.  Build only, no solves.
cuadmm_path;
pvar s1 s2 t

C = {};
C{end+1} = {'V0 constant (control)', @() mk(@(x,s1,s2) diff(x,s1,2)+diff(x,s2,2)+2*x, 2)};
C{end+1} = {'V1 var diffusion s1'  , @() mk(@(x,s1,s2) (1+s1)*diff(x,s1,2)+diff(x,s2,2)+2*x, 2)};
C{end+1} = {'V2 var diffusion both', @() mk(@(x,s1,s2) (1+s1)*diff(x,s1,2)+(1+s2)*diff(x,s2,2)+2*x, 2)};
C{end+1} = {'V3 var reaction'      , @() mk(@(x,s1,s2) diff(x,s1,2)+diff(x,s2,2)+2*(1+s1*s2)*x, 2)};
C{end+1} = {'V4 var advection'     , @() mk(@(x,s1,s2) diff(x,s1,2)+diff(x,s2,2)+(1+s1)*diff(x,s1,1)+2*x, 2)};
C{end+1} = {'V5 mixed order (1st s2)', @() mk1(@(x,s1,s2) diff(x,s1,2)+diff(x,s2,1)+2*x)};

fprintf('MP plant|Top.dim|T.R22{1,1}|A.R22{1,1}|Q.R22{1,1}|Q.R22 nonzero cells\n');
for k = 1:numel(C)
    try
        [Top,Aop] = C{k}{2}();
        np  = Top.dim(:,1);
        Iop = opvar2d(eye(sum(np)),[np np],Top.I,[Top.var1 Top.var2]);
        Qop = clean_opvar((Aop'*(Iop*Top))' + Aop'*(Iop*Top),1e-12);
        nzc = '';
        for i=1:3, for j=1:3
            if nz(Qop.R22{i,j}), nzc = [nzc sprintf('%d%d ',i,j)]; end
        end, end
        fprintf('MP %s|[%s]|%d|%d|%d|%s\n', C{k}{1}, num2str(np'), ...
                nz(Top.R22{1,1}), nz(Aop.R22{1,1}), nz(Qop.R22{1,1}), strtrim(nzc));
    catch ME
        fprintf('MP %s|FAILED|%s\n', C{k}{1}, strrep(ME.message,newline,' '));
    end
end
fprintf('MPDONE\n');

function [Top,Aop] = mk(rhsfun,~)
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0 1;0 1]);
sys = [diff(x,t,1) == rhsfun(x,s1,s2);
       subs(x,s1,0)==0; subs(x,s1,1)==0;
       subs(x,s2,0)==0; subs(x,s2,1)==0];
PIE = initialize(convert(sys));  Top = PIE.T; Aop = PIE.A;
end

function [Top,Aop] = mk1(rhsfun)
% first order in s2 -> only ONE boundary condition in that direction
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0 1;0 1]);
sys = [diff(x,t,1) == rhsfun(x,s1,s2);
       subs(x,s1,0)==0; subs(x,s1,1)==0;
       subs(x,s2,0)==0];
PIE = initialize(convert(sys));  Top = PIE.T; Aop = PIE.A;
end

function b = nz(p)
if isempty(p), b = 0; return; end
if isa(p,'double'), b = double(any(p(:)~=0)); return; end
b = double(~isempty(p.coefficient) && any(any(p.coefficient~=0)));
end
