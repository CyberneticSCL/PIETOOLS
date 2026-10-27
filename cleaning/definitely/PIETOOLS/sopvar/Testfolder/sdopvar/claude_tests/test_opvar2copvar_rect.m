%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_OPVAR2COPVAR_RECT checks 'opvar2copvar' on RECTANGULAR opvars,
% dim = [m1 m2; n1 n2], and on opvars with zero components - the PIE input
% and output operators B1, C1, D11 are of that kind. 'Test_copvar_converters'
% covers the square case only.
%
% Per case, against the semantics (CLAUDE.md S4), never a round trip:
%   ACTION   y = P x for a polynomial test function x, from the opvar
%            DEFINITION (P, Q1, Q2, R0, R1, R2 with 'int' and 'subs'),
%            against the blockwise action of the container ('apply_sopvar');
%   SHAPE    one row per nonzero-dimensional output space, one column per
%            input space, with the opvar's dimensions;
%   BLOCKS   a nonzero component is never [], and a block holding a zero
%            component exists only where its row or column has no other;
%   VALID    'verify' passes and copvar(C) rebuilds the same metadata from
%            the blocks alone (both failed for D11 = 0 before 09/25/2026).
%
% MMP, 09/25/2026: Initial coding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','sopvar:noncanonicalMultiplier');
pvar s1 s1_dum
dom = [0,1];    tol = 1e-9;     nchk = 0;
nm = {'P','Q1';'Q2','R'};

% dim, components set to zero, label
CASES = {
  [2 1; 3 2], {},            'full rectangular'
  [2 0; 0 3], {},            'Q1 only, L2^3 -> R^2 (C1-like)'
  [0 2; 3 0], {},            'Q2 only, R^2 -> L2^3 (B1-like)'
  [0 0; 2 3], {},            'R only, L2^3 -> L2^2'
  [2 3; 0 2], {},            'P and Q1, R^3 x L2^2 -> R^2'
  [1 2; 0 0], {'P'},         'P = 0 alone (D11 = 0)'
  [2 1; 3 2], {'P','Q1'},    'zero R row'
  [2 1; 3 2], {'P','Q2'},    'zero R column'
  [2 1; 3 2], {'Q1','Q2'},   'zero off-diagonal'
  [2 1; 3 2], {'P','Q1','Q2','R'}, 'all zero'
};
for ic = 1:size(CASES,1)
    [d,zc,lbl] = deal(CASES{ic,:});
    rng(100+ic);
    Pop = rand_opvar(d,2,s1,s1_dum,dom);
    for k = 1:numel(zc)
        c = Pop.(zc{k});
        if isstruct(c)
            c.R0 = 0*c.R0;  c.R1 = 0*c.R1;  c.R2 = 0*c.R2;
        else
            c = 0*c;
        end
        Pop.(zc{k}) = c;
    end
    Pm = opvar2copvar(Pop);
    rows = [d(1,1)>0, d(2,1)>0];    cols = [d(1,2)>0, d(2,2)>0];

    % SHAPE
    assert(isequal(size(Pm),[nnz(rows),nnz(cols)]),'%s: grid %s',lbl,mat2str(size(Pm)));
    assert(isequal(Pm.dim_out(:),d(rows,1)) && isequal(Pm.dim_in(:),d(cols,2)),'%s: dims',lbl);
    isL2 = [false,true];
    assert(isequal(any(Pm.space_out,2),isL2(rows)') && isequal(any(Pm.space_in,2),isL2(cols)'),...
        '%s: which rows/columns are L2',lbl);

    % BLOCKS
    ir = cumsum(rows);  jc = cumsum(cols);
    for i = find(rows)
        for j = find(cols)
            B = Pm.C{ir(i),jc(j)};
            z = ismember(nm{i,j},zc);
            if ~z
                assert(~isempty(B),'%s: nonzero component %s is []',lbl,nm{i,j});
            elseif ~isempty(B)
                % allowed only where its row or column has no nonzero one
                rz = all(ismember(nm(i,cols),zc));  cz = all(ismember(nm(rows,j),zc));
                assert(rz || cz,'%s: zero component %s kept as a block',lbl,nm{i,j});
            end
        end
    end

    % VALID
    v = verify(Pm);
    assert(v.true,'%s: verify fails: %s',lbl,strjoin(v.flags(:)',' | '));
    assert(isequal(metadata(copvar(Pm.C)),metadata(Pm)),'%s: blocks do not rebuild the container',lbl);

    % ACTION
    e = action_err(Pop,Pm,d,rows,cols,dom,s1,s1_dum);
    assert(e<tol,'%s: action differs by %.2e',lbl,e);
    nchk = nchk+6;
    fprintf('  passed: %-34s dim %-12s grid %dx%d\n',lbl,mat2str(d),size(Pm,1),size(Pm,2));
end
fprintf('test_opvar2copvar_rect passed (%d checks).\n',nchk);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = action_err(Pop,Pm,d,rows,cols,dom,v1,v2)
% y = P x from the opvar components, against the container blockwise. A
% zero-dimensional part contributes nothing and has no test function.
m1 = d(1,1); m2 = d(1,2); n1 = d(2,1); n2 = d(2,2);
xR = polynomial((1:m2)');
xL = polynomial(zeros(n2,1));
for k = 1:n2,   xL(k) = k + (k+1)*v1 + 0.5*v1^2;   end
yR = polynomial(zeros(m1,1));   yL = polynomial(zeros(n1,1));
if m1>0 && m2>0,    yR = yR + Pop.P*xR;     end
if m1>0 && n2>0
    Q1 = polynomial(Pop.Q1);
    if ismember(v2.varname{1},Q1.varname),  Q1 = subs(Q1,v2,v1);   end
    yR = yR + int(Q1*xL,v1,dom(1),dom(2));
end
if n1>0 && m2>0,    yL = yL + Pop.Q2*xR;    end
if n1>0 && n2>0
    xLt = subs(xL,v1,v2);
    yL = yL + Pop.R.R0*xL + int(Pop.R.R1*xLt,v2,dom(1),v1) + int(Pop.R.R2*xLt,v2,v1,dom(2));
end
% Container side: row i of the grid is output space find(rows)(i).
ro = find(rows);    co = find(cols);
X = {xR, xL};       Y = {polynomial(zeros(m1,1)), polynomial(zeros(n1,1))};
for a = 1:numel(ro)
    for b = 1:numel(co)
        B = Pm.C{a,b};
        if isempty(B),  continue,   end
        Y{ro(a)} = Y{ro(a)} + apply_sopvar(B,X{co(b)});
    end
end
e = 0;
if m1>0,    e = max(e,perr(yR,Y{1}));   end
if n1>0,    e = max(e,perr(yL,Y{2}));   end
end

function e = perr(A,B)
D = polynomial(A) - polynomial(B);
e = full(max(abs(D.coefficient(:))));
if isempty(e),  e = 0;  end
end
