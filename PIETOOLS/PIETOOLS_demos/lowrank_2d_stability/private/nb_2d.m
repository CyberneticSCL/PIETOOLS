function R = nb_2d(deg,dims,sep,excl)
% nb_2d(deg,dims,sep,excl) -- per-group OPERATOR BASIS COUNT for a 2-D
% poslpivar_2d Gram block.
%
% poslpivar_2d has SIXTEEN monomial groups where poslpivar has four
% (poslpivar_2d.m lines 562-606 build them, lines 762-890 stack them):
%    g1  identity              n0   1 var-free
%    g2  Zxo(s1)               nx   1 variable
%    g3  Zxa(s1,t1)            nx   2 variables
%    g4  Zxb(s1,t1)            nx   2
%    g5  Zyo(s2)               ny   1
%    g6  Zya(s2,t2)            ny   2
%    g7  Zyb(s2,t2)            ny   2
%    g8  Z2oo(s1,s2)           n2   2
%    g9  Z2ao(s1,t1,s2)        n2   3
%    g10 Z2bo(s1,t1,s2)        n2   3
%    g11 Z2oa(s1,s2,t2)        n2   3
%    g12 Z2ob(s1,s2,t2)        n2   3
%    g13 Z2aa(s1,t1,s2,t2)     n2   4 variables
%    g14 Z2ba(s1,t1,s2,t2)     n2   4
%    g15 Z2ab(s1,t1,s2,t2)     n2   4
%    g16 Z2bb(s1,t1,s2,t2)     n2   4
% Ten of the sixteen are over three or four variables, and a k-variable
% group holds O(d^k) monomials, which is the mechanism behind "2D has lots
% of operator bases even for very small degree".
%
% Gram side length N = sum_g n_g*|Z_g| = STATES x BASES; n_bases = sum_g |Z_g|
% over the groups that survive exclusion/separability.
%
% deg is the LF_deg/eq_deg struct (.dx,.dy,.d2); dims = [n0 nx ny n2].
pvar s1 t1 s2 t2
dx = deg.dx;  dy = deg.dy;  d2 = deg.d2;
n0=dims(1); nx=dims(2); ny=dims(3); n2=dims(4);

% variable tables, matching poslpivar_2d.m lines 563-606 EXACTLY -- the
% joint-degree entries of d2 are indexed column-major against this order,
% so a permuted table would silently miscount
V = { [],   s1,        [s1;t1],      [s1;t1], ...
      s2,   [s2;t2],   [s2;t2], ...
      [s1;s2], [s1;t1;s2], [s1;t1;s2], [s1;s2;t2], [s1;s2;t2], ...
      [s1;t1;s2;t2], [s1;t1;s2;t2], [s1;t1;s2;t2], [s1;t1;s2;t2] };
D = { [],   dx{1}, dx{2}, dx{3}, ...
      dy{1}, dy{2}, dy{3}, ...
      d2{1,1}, d2{2,1}, d2{3,1}, d2{1,2}, d2{1,3}, ...
      d2{2,2}, d2{3,2}, d2{2,3}, d2{3,3} };

nZ = zeros(1,16);  nZ(1) = 1;
for g = 2:16
    nZ(g) = size(nb_monoms(V{g},D{g}),1);
end

% group inclusion, replicating poslpivar_2d.m lines 426-481
exL = excl(:)';
if numel(exL)==4
    exL = [exL(1),exL(2:4),exL(2:4),reshape((exL(2:4)').*exL(2:4),1,[])];
end
sp = sep(:)';
if numel(sp)==1, sp = sp*ones(1,6); end %#ok<ISCL>
if sp(1)==1, exL(4)  = 1; end
if sp(2)==1, exL(7)  = 1; end
if sp(3)==1, exL(10) = 1; end
if sp(4)==1, exL(12) = 1; end
if sp(5)==1 && sp(6)
    exL(14:16) = 1;
elseif sp(5)==1
    exL([14,16]) = 1;
elseif sp(6)==1
    exL([15,16]) = 1;
end
if n0==0, exL(1)    = 1; end
if nx==0, exL(2:4)  = 1; end
if ny==0, exL(5:7)  = 1; end
if n2==0, exL(8:16) = 1; end
inc = ~exL;

ng = [n0 nx nx nx ny ny ny n2*ones(1,9)];
R.nZ      = nZ;
R.inc     = inc;
R.Zdim    = inc.*ng.*nZ;
R.N       = sum(R.Zdim);
R.n_bases = sum(inc.*nZ);
R.nvars   = [0 1 2 2 1 2 2 2 3 3 3 3 4 4 4 4];
end
