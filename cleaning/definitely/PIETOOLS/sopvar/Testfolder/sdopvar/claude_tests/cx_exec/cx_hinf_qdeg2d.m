function deg = cx_hinf_qdeg2d(Rm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEG = CX_HINF_QDEG2D(RM) sizes the free operator Q of the 2-D
% non-coercive executive from the container RM, as far as lpivar_cdopvar's
% vocabulary allows. The stock line is
%
%   PTdeg = get_lpivar_degs(Rop,Top);  lpivar(prog,Top.dim,PTdeg)
%
% whose 2-D branch hands lpivar_2d R's own max-degree table per component
% and per cell (d2{i,j} = Rmaxdegs.R22{i,j}, per-variable AND joint caps,
% no -1 adjustment). lpivar_cdopvar takes one cap per variable ROLE for a
% whole block - 'mult' (left degree, multiplier direction), 'int' [left
% right] (integral direction), 'out', 'in' - with no joint and no per-cell
% caps. The best available match is the per-role MAXIMUM of R's degrees,
% read here over every structurally nonzero coefficient of every block:
% a family that CONTAINS the stock one (nested, not equal). Per-cell zero
% structure (e.g. a cell of R that is identically zero, which lpivar_2d
% still gives a degree-0 parameter) and joint caps are not expressible.
% Kernel layout (I_q kron ZL)' C (I_p kron ZR), ZL over vars.out = [S2,S3],
% ZR over vars.in = [S3,S1], each a Kronecker product over its variables,
% first variable outer (sdopvar.m header).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deg = struct('mult',0,'int',[0 0],'out',0,'in',0);
for i = 1:size(Rm.C,1)
    for j = 1:size(Rm.C,2)
        B = Rm.C{i,j};
        if isempty(B),  continue,   end
        vo = B.vars.out;    vi = B.vars.in;
        S3 = intersect(vo,vi);      n3 = numel(S3);
        EL = kron_exps(B.ZL);       ER = kron_exps(B.ZR);
        nL = size(EL,1);    nR = size(ER,1);
        m = B.dims(1)*nL;   n = B.dims(2)*nR;
        [~,pL] = ismember(S3,vo);   [~,pR] = ismember(S3,vi);
        s2 = ~ismember(vo,S3);      s1 = ~ismember(vi,S3);
        for q = 1:numel(B.params.A)
            A = B.params.A{q};  nz = false(m*n,1);
            if numel(A)==m*n,   nz = nz | full(A(:)~=0);    end
            Bq = B.params.B{q};
            if size(Bq,2)==m*n, nz = nz | full(any(Bq~=0,1)).';    end
            [r,c] = find(reshape(nz,m,n));
            if isempty(r),  continue,   end
            eL = EL(mod(r-1,nL)+1,:);   eR = ER(mod(c-1,nR)+1,:);
            g = cell(1,max(n3,1));      [g{:}] = ind2sub(3*ones(1,max(n3,2)),q);
            for k = 1:n3
                if g{k}==1
                    deg.mult = max(deg.mult,max(eL(:,pL(k))));
                else
                    deg.int = max(deg.int,[max(eL(:,pL(k))), max(eR(:,pR(k)))]);
                end
            end
            if any(s2),     deg.out = max(deg.out,max(max(eL(:,s2))));  end
            if any(s1),     deg.in  = max(deg.in, max(max(eR(:,s1))));  end
        end
    end
end
end

function E = kron_exps(Z)
% Exponent rows of the Kronecker-product monomial basis Z{1} x Z{2} x ...
E = zeros(1,0);
for k = 1:numel(Z)
    z = Z{k}(:);
    E = [kron(E,ones(numel(z),1)), repmat(z,size(E,1),1)];
end
end
