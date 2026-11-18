function [rows,offs,pd,label] = pf_terms(setname,nv,o)
% [ROWS,OFFS,PD,LABEL] = PF_TERMS(SET,NV,O) the factor rows (POSCOPVAR_DIRECT
% psatz_rows: 0 none, 1 lower face, 2 upper face, 3 box quadratic) of a
% named Positivstellensatz term set in NV registry variables, one offset
% per row (O, default 1) and the per-direction flag PD (true: the offset
% lowers the weight in the acting directions only; false: in every
% direction, the scalar rule of the repo). The plain term is always
% present (psatz = 0) and is not a row here.
%   'none'      no rows
%   'box'       N rows with one 3: the Putinar set {1, g_k}
%   'product'   the 2^N - 1 nonzero rows over {0,3}^N: the box products
%   'faces'     2N rows with one entry in {1,2}
%   'facesS'    the faces with the scalar offset rule (pd false)
%   'corners'   the faces and the 4 C(N,2) two-face rows
%   'box0', 'product0', 'faces0', 'corners0': the same at offset 0
% SCRATCH, MMP, 11/15/2026 (psatz_faces test plan, Sec. 3).
if nargin<3 || isempty(o),  o = 1;  end
name = setname;
if endsWith(name,'0'),  o = 0;  name = name(1:end-1);   end
pd = true;
switch name
    case 'none',    rows = zeros(0,nv);
    case 'prodonly',rows = 3*ones(1,nv);    % the stock 2-D psatz = 1: g_1 ... g_N alone
    case 'box',     rows = 3*eye(nv);
    case 'product'
        rows = zeros(0,nv);
        for i = 1:2^nv-1
            rows(end+1,:) = 3*(bitget(i,1:nv)>0);                           %#ok<AGROW>
        end
    case {'faces','facesS'}
        rows = zeros(0,nv);
        for k = 1:nv
            for e = 1:2,    r = zeros(1,nv);    r(k) = e;    rows(end+1,:) = r;     end %#ok<AGROW>
        end
        if strcmp(name,'facesS'),   pd = false;     end
    case 'corners'
        rows = pf_terms('faces',nv,o);
        for i = 1:nv
            for j = i+1:nv
                for e = 1:2
                    for f = 1:2,    r = zeros(1,nv);    r(i) = e;   r(j) = f;   rows(end+1,:) = r;  end %#ok<AGROW>
                end
            end
        end
    otherwise,  error('pf_terms:set','unknown term set %s',setname);
end
offs = o*ones(1,size(rows,1));
label = sprintf('%s_o%d_pd%d',name,o,pd);
end
