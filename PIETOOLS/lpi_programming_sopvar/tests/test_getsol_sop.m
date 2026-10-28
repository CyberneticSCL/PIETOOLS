function R = test_getsol_sop(parts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_GETSOL_SOP(PARTS) tests GETSOL_LPIVAR_SOP, LPIGETSOL_SOP and
% SUBS_DVAR_SOP, which replace the decision variables of an LPI object by
% their values, and errors at the first failed check. Operators are
% compared by their action on random polynomial test functions, the
% kernels integrated by Gauss-Legendre quadrature of the 'sopvar' class
% definition (HEATND_APPLY):
%   (i)   9 random 'cdopvar' objects: 1, 2 and 3 spatial variables on
%         [0,1], [-1,2], [0.5,1.5], R^n spaces, structurally empty blocks,
%         a decision-variable list per block merged by the constructor, a
%         'sopvar' block among 'sdopvar' blocks; the values given by name,
%         shuffled, with 4 names the operator does not have and RRx 3
%         entries longer than decvartable. Each extracted block, applied to
%         a test function (4 terms of degree at most 2 per variable) at 3
%         points, equals to 1e-10 relative the action of the 'sdopvar'
%         kernels converted by SDVAR2DPVAR, fixed by SOSGETSOL and
%         integrated separately, 6 Gauss nodes per direction in both; with
%         the values deranged among the names, every block differs by more
%         than 1e-3 relative. The three
%         routines return the same 'copvar', bit for bit, which passes
%         verify, keeps the container metadata, the empty blocks and the
%         'sopvar' block of the input, and raises no canonical-form
%         warning. A random half of the names gives a 'cdopvar' on the
%         others, in order, with their B rows unchanged; the others then
%         give the full result to 1e-13.
%   (ii)  The stability LPI of CX_PIE2PDESTABILITY (T'*Q - P = 0, D + N = 0
%         with D = A'*Q + Q'*A + epneg*P, P and N positive, Q free) for
%         x_t = x_ss + pi^2/2 x, x(0) = x(1) = 0 (CX_PLANT('rd',0.5)), at
%         lpisettings('light') with MOSEK. The solve has pinf = numerr = 0,
%         the equality residual rel_b on solinfo.RRx at most 1e-6, and Gram
%         eigenvalues at least -1e-8 times the largest. The extracted P, Q,
%         D, N pass verify and agree bit for bit across the three routines.
%         On 6 test functions (5 monomials of degree at most 3 per
%         variable; kernels by 10 Gauss nodes, L2 norms by 12) the relative
%         residuals of T'*Q = P and D = -N are at most 1e-6, the extracted
%         D equals A'*Q + Q'*A + epneg*P of the extracted operators to
%         1e-10 relative, and <x,Px>/(||x|| ||Px||) and <x,Nx>/(||x||
%         ||Nx||) are at least -1e-8. RRx plus 1e-3*norm(RRx)/sqrt(n) times
%         a standard normal vector, n = numel(decvartable), raises the
%         larger residual by a factor of more than 1e3.
%   (iii) A 1-D program (s in [0,1]; RRx random, 2 entries longer than
%         decvartable) and a 2-D program on [0,1]^2. LPIGETSOL_SOP returns
%         the result of LPIGETSOL bit for bit for [], 'double',
%         'polynomial', 'opvar', 'dpvar', a 'dpvar' expression, char,
%         cellstr, a cell of 'dpvar', two 'dopvar', a 'dopvar2d' and an
%         'opvar2d' object, and so does SUBS_DVAR_SOP for the 1-D objects;
%         the comparison detects RRx + 1; GETSOL_LPIVAR_SOP returns the
%         result of GETSOL_LPIVAR for 'dopvar' and 'dopvar2d'; a permuted
%         name list gives the same 'dpvar' result; a name listed twice in
%         decvartable takes the value of its lower row index, in
%         GETSOL_LPIVAR_SOP and SOSGETSOL alike. A 'cdopvar' gives a
%         'copvar' and an 'sdopvar' a 'sopvar'; a 'copvar' and a 'sopvar'
%         are returned unchanged; a cell mixing both families is solved
%         entry by entry; a 'cdopvar' whose names the program does not list
%         is returned with the same names. A struct, logical, int8, string
%         or (when the class exists) 'nopvar' input raises
%         'lpigetsol_sop:badClass', and a cell holding a struct raises an
%         error; so do an unsolved program (empty solinfo.info) and, in
%         GETSOL_LPIVAR_SOP, a program without solinfo.RRx
%         ('getsol_lpivar_sop:notSolved').
%   (iv)  The H-infinity gain LPI for x_t = x_ss + pi^2/2 x + w,
%         z = int_0^1 x ds, x(0) = x(1) = 0, w in L2[0,1] (CX_PLANT('io1')),
%         at lpisettings('light') with MOSEK: gamma a 'dpvar' decision
%         variable, gamma >= 0, minimized; T'*Q - R = 0 with R positive and
%         Q free; N + K = 0 with N positive, for the KYP operator
%           K = [-gamma*I, Dzw', Bw'*Q; Dzw, -gamma*I, Cz;
%                Q'*Bw, Cz', A'*Q + Q'*A],
%         a program in which solinfo.x and solinfo.RRx differ in length.
%         The solve has pinf = numerr = 0, rel_b at most 1e-6 and the Gram
%         eigenvalue ratio at least -1e-8; the extracted K equals K rebuilt
%         from the gamma of LPIGETSOL and the extracted Q to 1e-10
%         relative; on 6 test functions T'*Q = R holds to 1e-5 relative and
%         <x,Kx>/(||x|| ||Kx||) is at most 1e-8.
% In (ii) and (iv) the residuals with solinfo.x in place of RRx are
% printed, not checked.
%
% INPUT
% - parts:  (optional) 'cellstr' or char, the parts to run, from
%           {'i','ii','iii','iv'}; default (or []) all four.
%
% OUTPUT
% - R:      'struct' with one field per part run:
%             i    'struct' array, one entry per case, with fields label
%                  (char), q (double, the number of decision variables),
%                  nblk (double, the blocks compared), err (double, the
%                  largest relative error), err_ctrl (double, the smallest
%                  relative error of the control) and err_part (double, the
%                  largest change of a kept B row);
%             ii   'struct' with fields ndec (double, numel(decvartable)),
%                  build, solve, extract (double, seconds), rel_b,
%                  psd_relmin, rel1 (T'*Q = P), rel2 (D = -N), relD, psdP,
%                  psdN (double, the quantities above) and controls (1 x 2
%                  'struct' array with fields name, rel1, rel2);
%             iii  'struct' with fields nlegacy, nnew, nerr (double, the
%                  legacy objects compared, the new-class checks, the
%                  rejected inputs);
%             iv   'struct' with fields gam (double, gamma), rel_b, relK,
%                  rel1, maxK, rel_x (double, the quantities above and the
%                  residual with solinfo.x), ndec, nx (double,
%                  numel(decvartable) and numel(solinfo.x));
%           printed lines; errors at the first failure.
%
% NOTES
% Parts (ii) and (iv) require MOSEK, and CX_PLANT, CX_SETTINGS,
% CX_PIE2PDESTABILITY and CX_SOLVE (sopvar/Testfolder/sdopvar/claude_tests/
% cx_exec) on the path.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_getsol_sop
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, 09/29/2026
% MMP, 10/06/2026: part (iv): cx_hinf_lf, cx_hinf_qdeg, cx_space_list and
%                  cx_hinf_slack -> library poslpivar_settings_sop ('lf';
%                  'slack' + the equality cx_hinf_slack imposed),
%                  get_lpivar_degs_sop, copvar_space_list, so the test runs
%                  on the library translators (lpi_programming_sopvar);
%                  programs unchanged (126-program bit-identity check).
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                  inputs, outputs), the licence block added, step comments
%                  in the body and headers of the local functions; code
%                  unchanged. The deleted header text was a notes
%                  reference, an author attribution, the rationale for the
%                  independent evaluation, a reference to the Tier 1 plan,
%                  and the claim that the solinfo.x control is checked in
%                  (ii) (it is printed only).

% Default: all four parts; one part given as char
if nargin<1 || isempty(parts),  parts = {'i','ii','iii','iv'};  end
if ischar(parts),   parts = {parts};    end
% Each requested part from its own fixed seed
R = struct();
if any(strcmp(parts,'i')),      rng(20260929);  R.i   = part_i();     end
if any(strcmp(parts,'ii')),     rng(20260930);  R.ii  = part_ii();    end
if any(strcmp(parts,'iii')),    rng(20261001);  R.iii = part_iii();   end
if any(strcmp(parts,'iv')),     rng(20261002);  R.iv  = part_iv();    end
fprintf('test_getsol_sop: parts {%s} passed.\n',strjoin(parts,','));
end


% ========================================================================
% (i) random operators against an independent evaluation
% ========================================================================
function out = part_i()
% OUT = PART_I() runs part (i) on the nine random 'cdopvar' cases below,
% and errors at the first failed check.
% INPUT  none.
% OUTPUT out: 'struct' array, one entry per case, with fields label, q,
%             nblk, err, err_ctrl and err_part (R.i of TEST_GETSOL_SOP).
% Tolerance, quadrature, the variables and their domains; sp(o,i) the
% spaces struct of output spaces o and input spaces i over V
tol = 1e-10;    nq = 6;     % Gauss nodes: exact to degree 11, integrands here <= 4
V = {'s1','s2','s3'};   D = [0 1; -1 2; 0.5 1.5];   % distinct domains per variable
sp = @(o,i) struct('out',{o},'in',{i},'vars',{V});
% {label, spaces, dims, deg, ndec, occ, sopvar block (i,j) or []}
cases = {
  {'1-D s1->s1',            sp({{'s1'}},{{'s1'}}),                     struct('out',2,'in',3),         2, 9,       [],  []}
  {'1-D s2->s1 (S1,S2)',    sp({{'s1'}},{{'s2'}}),                     struct('out',2,'in',1),         2, 5,       [],  []}
  {'2-D s1s2->s1s2',        sp({{'s1','s2'}},{{'s1','s2'}}),           struct('out',1,'in',2),         1, 6,       [],  []}
  {'s2s3->s1s2',            sp({{'s1','s2'}},{{'s2','s3'}}),           struct('out',2,'in',1),         1, 6,       [],  []}
  {'3 vars s1s2s3',         sp({{'s1','s2','s3'}},{{'s1','s2','s3'}}), struct('out',1,'in',1),         1, 5,       [],  []}
  {'R^n x L2 2x2',          sp({{},{'s1'}},{{},{'s1'}}),               struct('out',[2;1],'in',[3;1]), 2, 7,       [],  []}
  {'3x3 empty blocks',      sp({{},{'s1'},{'s1','s2'}},{{'s1'},{},{'s2','s3'}}), ...
                            struct('out',[2;1;1],'in',[1;2;1]),        1, 8, logical([1 0 1;1 1 0;0 1 1]), []}
  {'2x2 per-block lists',   sp({{'s1'},{}},{{'s1'},{'s2'}}),           struct('out',[1;2],'in',[2;1]), 2, [3 4;5 2],[],  []}
  {'2x2 sopvar block',      sp({{'s1'},{}},{{'s1'},{'s2'}}),           struct('out',[1;2],'in',[2;1]), 2, 6,       [],  [1 2]}
  };
out = struct('label',{},'q',{},'nblk',{},'err',{},'err_ctrl',{},'err_part',{});
for c = 1:numel(cases)
    % The random cdopvar of the case (parameter density 0.5), a fixed block
    % put in place of block sblk when the case has one
    cs = cases{c};
    [lbl,spaces,dims,deg,ndec,occ,sblk] = cs{:};
    P = rand_cdopvar(spaces,dims,D,deg,ndec,0.5,occ);
    if ~isempty(sblk)                   % a fixed block among decision blocks
        F = rand_copvar(spaces,dims,D,deg,0.5,occ);
        C = P.C;    C{sblk(1),sblk(2)} = F.C{sblk(1),sblk(2)};
        P = cdopvar(C);
        assert(isa(P.C{sblk(1),sblk(2)},'sopvar'),'case %s: sopvar block not kept',lbl)
    end
    assert(verify(P).true,'case %s: generator produced an invalid cdopvar',lbl)
    % A random value for every decision variable of P
    Zd = cellstr(P.Zd(:));      q = numel(Zd);
    d = randn(q,1);
    % Shuffled list, padded with names the operator does not have.
    extra = strcat('zz_absent_',cellstr(string((1:4)')));
    names = [Zd; extra];        vals = [d; randn(4,1)];
    pm = randperm(numel(names));    names = names(pm);  vals = vals(pm);
    % A solved program with that list; RRx 3 entries longer than it
    prog = struct('decvartable',{names},'solinfo',struct('RRx',[vals; randn(3,1)],'info',struct('pinf',0)));

    % The three entry points: copvar results, no canonical-form warning,
    % verify, identical results and the container metadata of P
    lastwarn('');
    S1 = subs_dvar_sop(P,names,vals);
    S2 = getsol_lpivar_sop(prog,P);
    S3 = lpigetsol_sop(prog,P);
    [wm,wid] = lastwarn;
    assert(~contains(wid,'noncanonical'),'case %s: result was rewritten (%s)',lbl,wm)
    assert(isa(S1,'copvar') && isa(S2,'copvar') && isa(S3,'copvar'),'case %s: not copvar',lbl)
    assert(verify(S1).true,'case %s: result fails verify',lbl)
    assert(same_ops(S1,S2,0) && same_ops(S1,S3,0),'case %s: entry points disagree',lbl)
    assert(isequal(S1.vars,P.vars) && isequal(S1.dom,P.dom) && isequal(S1.space_out,P.space_out) ...
        && isequal(S1.space_in,P.space_in) && isequal(S1.dim_out,P.dim_out) && isequal(S1.dim_in,P.dim_in), ...
        'case %s: container metadata changed',lbl)

    % Control: the same values shuffled among the operator names.
    dperm = d(derangement(q));
    Sbad = subs_dvar_sop(P,Zd,dperm);

    % Every block: empty blocks and the sopvar block kept; each sdopvar
    % block against the independent evaluation and the control
    [M,N] = size(P.C);  err = 0;    errc = inf;     nblk = 0;
    for i = 1:M
        for j = 1:N
            b = P.C{i,j};   bs = S1.C{i,j};
            if isempty(b)
                assert(isempty(bs),'case %s: zero block (%d,%d) not kept',lbl,i,j);    continue
            end
            if isa(b,'sopvar')
                assert(isequal(bs,b),'case %s: sopvar block (%d,%d) changed',lbl,i,j); continue
            end
            assert(isa(bs,'sopvar'),'case %s: block (%d,%d) is %s',lbl,i,j,class(bs))
            % the standalone sdopvar path gives the same block
            Sb = subs_dvar_sop(b,names,vals);
            assert(isa(Sb,'sopvar') && same_ops(Sb,bs,0),'case %s: sdopvar path differs at (%d,%d)',lbl,i,j)
            % A test function and points; the extracted block, the
            % independent evaluation and the control at them; relative
            % errors, on a nonzero independent value
            [f,S] = test_data(b,D,V);
            y_get = heatNd_apply(bs,f,S,nq);
            y_ind = indep_apply(b,names,vals,f,S,nq);
            y_bad = heatNd_apply(Sbad.C{i,j},f,S,nq);
            nr = norm(y_ind,'fro');
            assert(nr>1e-6,'case %s: vacuous comparison at (%d,%d)',lbl,i,j)
            err = max(err,norm(y_get-y_ind,'fro')/nr);
            errc = min(errc,norm(y_bad-y_ind,'fro')/nr);
            nblk = nblk+1;
        end
    end
    % Largest error below tol; the control detected in every block
    assert(err<tol,'case %s: extracted operator differs from the independent evaluation by %.3g',lbl,err)
    assert(errc>1e-3,'case %s: control not detected (%.3g) - the comparison cannot see a wrong d',lbl,errc)

    % Partial substitution: a random half of the names given.
    giv = randperm(q,floor(q/2));   rest = setdiff(1:q,giv);
    Sp = subs_dvar_sop(P,Zd(giv),d(giv));
    assert(isa(Sp,'cdopvar') && verify(Sp).true,'case %s: partial result invalid',lbl)
    assert(isequal(cellstr(Sp.Zd(:)),Zd(rest)),'case %s: remaining names wrong or reordered',lbl)
    % Largest change of a kept B row, over the decision blocks and cells
    ep = 0;
    for k = 1:numel(P.C)
        if ~isa(P.C{k},'sdopvar'),  continue,   end
        Bo = P.C{k}.params.B;   Bp = Sp.C{k}.params.B;
        for g = 1:numel(Bo)
            ep = max(ep,norm(Bp{g}-Bo{g}(rest,:),1));   % kept rows are the original rows
        end
    end
    Sf = subs_dvar_sop(Sp,Zd,d);                        % then the rest
    assert(isa(Sf,'copvar') && same_ops(Sf,S1,1e-13),'case %s: two-step substitution differs',lbl)
    assert(ep==0,'case %s: kept B rows changed',lbl)

    % Record and print the case
    out(end+1) = struct('label',lbl,'q',q,'nblk',nblk,'err',err,'err_ctrl',errc,'err_part',ep); %#ok<AGROW>
    fprintf('(i) %-22s q=%3d blocks=%d  rel err %.2e  control %.2e\n',lbl,q,nblk,err,errc);
end
end


% ------------------------------------------------------------------------
function [f,S] = test_data(b,D,V)
% [F,S] = TEST_DATA(B,D,V) returns a random polynomial test function on the
% block's input space (handle of points whose columns are in b.vars.in
% order) and three evaluation points of its output space (columns in
% b.vars.out order), the domains looked up by variable NAME.
% INPUT  b: 'sdopvar' object; D: nv x 2 double, the domain of each
%           variable of V; V: 1 x nv cellstr, the variable names.
% OUTPUT f: function_handle, f(X) for X an n x numel(b.vars.in) double
%           returns n x b.dims(2) double: 4 terms of degree at most 2 per
%           variable, the first constant, N(0,1) coefficients;
%        S: 3 x numel(b.vars.out) double, uniform in the domains, or 1 x 0
%           when b has no output variable.
% Exponents and coefficients of the terms
vin = b.vars.in(:)';    vout = b.vars.out(:)';
p = b.dims(2);
nt = 4;                                         % terms, degree <= 2 per variable
F.deg = randi([0 2],nt,numel(vin));     F.deg(1,:) = 0;
F.coef = randn(nt,p);
f = @(X) prod_terms(X,F);
% Three uniform points in the domain of each output variable
S = zeros(3,numel(vout));
for k = 1:numel(vout)
    dk = D(strcmp(V,vout{k}),:);
    S(:,k) = dk(1) + (dk(2)-dk(1))*rand(3,1);
end
if isempty(vout),   S = zeros(1,0);     end
end

function y = prod_terms(X,F)
% Y = PROD_TERMS(X,F) returns the polynomial with exponent rows F.deg and
% coefficient rows F.coef at the points X: sum_t prod_k
% X(:,k).^F.deg(t,k) * F.coef(t,:).
% INPUT  X: n x k double; F: 'struct' with fields deg (nt x k double) and
%           coef (nt x p double).
% OUTPUT y: n x p double.
n = size(X,1);  Z = ones(n,size(F.deg,1));
for t = 1:size(F.deg,1)
    for k = 1:size(F.deg,2),    Z(:,t) = Z(:,t).*X(:,k).^F.deg(t,k);    end
end
y = Z*F.coef;
end


% ------------------------------------------------------------------------
function y = indep_apply(b,names,vals,f,S,nq)
% Y = INDEP_APPLY(B,NAMES,VALS,F,S,NQ) returns (P(d) x)(s) for an
% 'sdopvar' block b, with its kernels obtained without the code under
% test: SDVAR2DPVAR gives each gamma cell as a 'dpvar' in (vars.out,
% vars.in_dum); SOSGETSOL substitutes the values by name; the polynomial
% kernels are integrated here. Gamma cell g runs over the SORTED shared
% variables, direction 1 fastest; 1 = multiplier (s' = s), 2 = int_a^s,
% 3 = int_s^b (sopvar document Sec. 4).
% INPUT  b:     'sdopvar' object; names: cellstr, decision variable names;
%        vals:  double, vals(k) the value of names{k};
%        f:     function_handle of TEST_DATA; S: ns x numel(b.vars.out)
%               double, the points; nq: double, Gauss nodes per direction.
% OUTPUT y:     ns x b.dims(1) double.
% Variables, dummy input names, parameters and names of the block; a
% program holding the values, for SOSGETSOL
vin = b.vars.in(:)';    vout = b.vars.out(:)';
vdum = strcat(vin,'_dum');
prm = b.params;     Zd = cellstr(b.Zd(:));
sos = struct('decvartable',{names(:)},'solinfo',struct('RRx',vals(:)));
% Shared (S3, sorted) and input-only (S1) variables, by position
S3 = sort(intersect(vin,vout));     n3 = numel(S3);
[~,i3in] = ismember(S3,vin);        [~,i3out] = ismember(S3,vout);
S1 = setdiff(vin,vout);             [~,i1in] = ismember(S1,vin);
q = b.dims(1);  p = b.dims(2);
[xg,wg] = gauss01(nq);
ns = size(S,1);     y = zeros(ns,q);
for g = 1:numel(prm.A)
    % The kernel of cell g at the values, and its multi-index gamma
    Dg = sdvar2dpvar(struct('A',prm.A{g},'B',prm.B{g}),b.dims, ...
                     struct('out',{vout},'in',{vdum}),b.ZL,b.ZR,Zd);
    K = sosgetsol(sos,Dg);                      % q x p, polynomial or double
    gam = ones(1,n3);
    if n3>0,    c = cell(1,n3);     [c{:}] = ind2sub([3*ones(1,n3),1],g);   gam = [c{1:n3}];    end
    for is = 1:ns
        % Integration limits at the point s: [a, s_k] for gamma_k = 2,
        % [s_k, b] for gamma_k = 3, the domain of an input-only variable
        s = S(is,:);
        lo = [];    hi = [];    col = [];
        for k = 1:n3
            if gam(k)==2,       lo(end+1) = b.dom.in(i3in(k),1); hi(end+1) = s(i3out(k));        col(end+1) = i3in(k); %#ok<AGROW>
            elseif gam(k)==3,   lo(end+1) = s(i3out(k));         hi(end+1) = b.dom.in(i3in(k),2); col(end+1) = i3in(k); %#ok<AGROW>
            end
        end
        for k = 1:numel(S1)
            lo(end+1) = b.dom.in(i1in(k),1);   hi(end+1) = b.dom.in(i1in(k),2);   col(end+1) = i1in(k); %#ok<AGROW>
        end
        % Tensor Gauss points of the integrated variables; a multiplier
        % variable (gamma_k = 1) fixed at s_k
        nd = numel(col);    nx = nq^nd;
        X = zeros(nx,numel(vin));   w = ones(nx,1);
        for k = 1:n3
            if gam(k)==1,   X(:,i3in(k)) = s(i3out(k));   end
        end
        sub = cell(1,nd);
        if nd>0,    [sub{:}] = ind2sub([nq*ones(1,nd),1],(1:nx)');  end
        for t = 1:nd
            X(:,col(t)) = lo(t) + (hi(t)-lo(t))*xg(sub{t});
            w = w.*((hi(t)-lo(t))*wg(sub{t}));
        end
        % Kernel and test function at the points; the weighted sum
        Kv = poly_at(K,[vout,vdum],[repmat(s,nx,1),X]);    % nx x q*p, column-major
        fx = f(X);
        for ii = 1:q
            for jj = 1:p
                y(is,ii) = y(is,ii) + sum(w.*Kv(:,ii+q*(jj-1)).*fx(:,jj));
            end
        end
    end
end
end

function v = poly_at(K,names,Y)
% V = POLY_AT(K,NAMES,Y) returns the values of a matrix polynomial (or
% double) at the rows of Y, whose columns are the variables NAMES; one
% column per entry, column-major.
% INPUT  K: 'polynomial' or 'double' matrix; names: cellstr, containing
%           the variables of K; Y: n x numel(names) double.
% OUTPUT v: n x numel(K) double.
% A double: the same values at every point
n = size(Y,1);
if isa(K,'double')
    v = repmat(reshape(full(K),1,[]),n,1);  return
end
% The monomials of K at the points, times its coefficients
dm = full(K.degmat);    cf = full(K.coefficient);   vn = K.varname;
[tf,loc] = ismember(vn,names);
assert(all(tf),'poly_at: kernel variable not among %s',strjoin(names,','))
Z = ones(n,size(dm,1));
for t = 1:size(dm,1)
    for k = 1:numel(loc),   Z(:,t) = Z(:,t).*Y(:,loc(k)).^dm(t,k);  end
end
v = Z*cf;
end


% ========================================================================
% (ii) a real solved container stability program
% ========================================================================
function out = part_ii()
% OUT = PART_II() runs part (ii): builds and solves the stability LPI of
% CX_PIE2PDESTABILITY on CX_PLANT('rd',0.5) with MOSEK, extracts its
% decision operators and checks them, and errors at the first failure.
% INPUT  none.
% OUTPUT out: 'struct' with fields ndec, build, solve, rel_b, psd_relmin,
%             extract, rel1, rel2, relD, psdP, psdN and controls (R.ii of
%             TEST_GETSOL_SOP).
% Gauss nodes of the kernel integrals and of the output grid; the number
% of test functions
nqi = 10;   nqo = 12;   ntest = 6;
% Build and solve the LPI; print the solve, then require MOSEK's
% certificate, the equality residual and PSD Gram blocks
PIE = cx_plant('rd',0.5);
st = cx_settings('light','mosek');
t0 = tic;   [prog,ops] = cx_PIE2PDEstability(PIE,st);    tb = toc(t0);
r = cx_solve(prog,st.sos_opts);
sol = r.sol;
fprintf('(ii) rd 0.5: build %.1f s, solve %.1f s wall, pinf %g dinf %g numerr %g, rel_b %.2e, psd_min %.2e (rel %.2e)\n',...
    tb,r.wall,r.pinf,r.dinf,r.numerr,r.rel_b,r.psd_min,r.psd_relmin);
assert(r.st==1,'(ii) MOSEK did not certify feasibility (pinf %g, numerr %g)',r.pinf,r.numerr)
assert(r.rel_b<=1e-6,'(ii) equality residual rel_b = %.3g on solinfo.RRx rows',r.rel_b)
assert(r.psd_relmin>=-1e-8,'(ii) Gram blocks not PSD: min eig / max eig = %.3g',r.psd_relmin)
% Print the sizes of decvartable, RRx and x, and the distance of x from
% RRx on their common entries
ntot = numel(sol.decvartable);
nxr = min(ntot,numel(sol.solinfo.x));
fprintf('(ii) %d decision variables, RRx %d, x %d; ||x-RRx|| = %.3g on the first %d, ||RRx|| = %.3g\n',ntot,...
    numel(sol.solinfo.RRx),numel(sol.solinfo.x),norm(sol.solinfo.x(1:nxr)-sol.solinfo.RRx(1:nxr)),nxr,...
    norm(sol.solinfo.RRx(1:ntot)));

% Extract every decision operator, through all three entry points.
nm = {'P','Q','D','N'};     X = struct();
t0 = tic;
for k = 1:numel(nm)
    Xk = getsol_lpivar_sop(sol,ops.(nm{k}));
    assert(isa(Xk,'copvar') && verify(Xk).true,'(ii) %s: extracted operator invalid',nm{k})
    assert(same_ops(Xk,lpigetsol_sop(sol,ops.(nm{k})),0),'(ii) %s: lpigetsol_sop differs',nm{k})
    assert(same_ops(Xk,subs_dvar_sop(ops.(nm{k}),sol.decvartable,sol.solinfo.RRx(1:ntot)),0),...
        '(ii) %s: subs_dvar_sop differs',nm{k})
    X.(nm{k}) = Xk;
end
tx = toc(t0);
% The fixed operators T and A; the decay rate epneg of the settings
T = ops.T;  A = ops.A;
epn = st.epneg;
% Test functions on T's input space (P, Q, N all live there).
F = test_funcs(X.P,ntest);
G = out_grid(X.P,nqo);
% Relation 1: T'*Q - P = 0 (lpi_eq_cdopvar, not symmetric).
L1 = apply_c(T'*X.Q,F,G,nqi);    R1 = apply_c(X.P,F,G,nqi);
% Relation 2: D + N = 0, D = A'*Q + Q'*A + epneg*P ('symmetric').
Dq = A'*X.Q + X.Q'*A;   if epn~=0,  Dq = Dq + epn*X.P;  end
L2 = apply_c(X.D,F,G,nqi);      R2 = apply_c(X.N,F,G,nqi);    L2q = apply_c(Dq,F,G,nqi);
rel1 = rel_res(L1,R1,G);        rel2 = rel_res(L2,neg(R2),G);
relD = rel_res(L2,L2q,G);       % getsol(D) against D rebuilt from getsol(Q)
% PSD on test functions: <x, P x> and <x, N x>, relative to ||x|| ||P x||.
pP = min_quad(F,R1,G);      pN = min_quad(F,R2,G);
% Print, then check the tolerances
fprintf('(ii) extract %.2f s | T''Q-P: rel %.2e | D+N: rel %.2e | D vs A''Q+Q''A: rel %.2e | min <x,Px> %.2e, <x,Nx> %.2e\n',...
    tx,rel1,rel2,relD,pP,pN);
assert(rel1<=1e-6 && rel2<=1e-6,'(ii) extracted operators violate the equality constraints (%.3g, %.3g)',rel1,rel2)
assert(relD<=1e-10,'(ii) getsol(D) differs from A''*getsol(Q)+getsol(Q)''*A by %.3g',relD)
assert(pP>=-1e-8 && pN>=-1e-8,'(ii) extracted P or N not PSD on test functions (%.3g, %.3g)',pP,pN)

% Controls: each must break the relations, else the check above is blind.
ctl = struct('name',{'solinfo.x','RRx+1e-3'},'rel1',{NaN,NaN},'rel2',{NaN,NaN});
% solinfo.x in place of RRx (zero-padded to numel(decvartable)), and RRx
% plus 1e-3*norm(RRx)/sqrt(ntot) times a standard normal vector
xx = zeros(ntot,1);     nx = min(ntot,numel(sol.solinfo.x));    xx(1:nx) = sol.solinfo.x(1:nx);
bad = sol;  bad.solinfo.RRx = xx;
bad2 = sol; bad2.solinfo.RRx = sol.solinfo.RRx + 1e-3*norm(sol.solinfo.RRx)/sqrt(ntot)*randn(size(sol.solinfo.RRx));
bb = {bad,bad2};
% Both relations at each control
for k = 1:2
    Qb = getsol_lpivar_sop(bb{k},ops.Q);    Pb = getsol_lpivar_sop(bb{k},ops.P);
    Nb = getsol_lpivar_sop(bb{k},ops.N);    Db = getsol_lpivar_sop(bb{k},ops.D);
    ctl(k).rel1 = rel_res(apply_c(T'*Qb,F,G,nqi),apply_c(Pb,F,G,nqi),G);
    ctl(k).rel2 = rel_res(apply_c(Db,F,G,nqi),neg(apply_c(Nb,F,G,nqi)),G);
    fprintf('(ii) control %-9s: T''Q-P rel %.2e, D+N rel %.2e\n',ctl(k).name,ctl(k).rel1,ctl(k).rel2);
end
% Only the perturbed RRx is asserted; the solinfo.x control is printed
assert(max(ctl(2).rel1,ctl(2).rel2)>1e3*max(rel1,rel2),'(ii) perturbed RRx not detected')
% Results
out = struct('ndec',ntot,'build',tb,'solve',r.wall,'rel_b',r.rel_b,'psd_relmin',r.psd_relmin,...
    'extract',tx,'rel1',rel1,'rel2',rel2,'relD',relD,'psdP',pP,'psdN',pN,'controls',ctl);
end


% ------------------------------------------------------------------------
function F = test_funcs(P,nt)
% F = TEST_FUNCS(P,NT) returns nt random polynomial test functions on the
% column spaces of container P: F{t}{j} = struct(vars, deg, coef), 5
% monomials of degree <= 3 per variable, N(0,1) coefficients.
% INPUT  P: 'copvar' object; nt: double.
% OUTPUT F: 1 x nt cell; F{t} 1 x N cell of 'struct' with fields vars
%           (cellstr, the variables of input space j), deg (5 x numel(vars)
%           double, the exponents) and coef (5 x P.dim_in(j) double).
[~,N] = size(P.C);
F = cell(1,nt);
for t = 1:nt
    F{t} = cell(1,N);
    for j = 1:N
        vj = P.vars(P.space_in(j,:));
        F{t}{j} = struct('vars',{vj},'deg',randi([0 3],5,numel(vj)),'coef',randn(5,P.dim_in(j)));
    end
end
end

function G = out_grid(P,nq)
% G = OUT_GRID(P,NQ) returns the Gauss grid (points and weights) on each
% row space of P; R^n: one point.
% INPUT  P: 'copvar' object; nq: double, the nodes per direction.
% OUTPUT G: 1 x M cell; G{i} 'struct' with fields vars (cellstr), X
%           (nq^nd x nd double, the points) and w (nq^nd x 1 double, the
%           weights), nd the number of variables of row space i.
[M,~] = size(P.C);
[xg,wg] = gauss01(nq);
G = cell(1,M);
for i = 1:M
    % Variables and domains of space i; the tensor index of every point,
    % then the nodes and weights mapped to each domain
    vi = P.vars(P.space_out(i,:));      di = P.dom(P.space_out(i,:),:);
    nd = numel(vi);     nx = nq^nd;
    Xp = zeros(nx,nd);  w = ones(nx,1);
    sub = cell(1,nd);
    if nd>0,    [sub{:}] = ind2sub([nq*ones(1,nd),1],(1:nx)');  end
    for t = 1:nd
        Xp(:,t) = di(t,1) + (di(t,2)-di(t,1))*xg(sub{t});
        w = w.*((di(t,2)-di(t,1))*wg(sub{t}));
    end
    G{i} = struct('vars',{vi},'X',Xp,'w',w);
end
end

function Y = apply_c(P,F,G,nq)
% Y = APPLY_C(P,F,G,NQ) returns Y{t}{i} = ((P x_t)_i at G{i}.X) =
% sum_j heatNd_apply(P.C{i,j}, x_t,j).
% INPUT  P:  'copvar' object; F: test functions (TEST_FUNCS) on the input
%            spaces of P; G: grids (OUT_GRID) on the output spaces of P;
%        nq: double, the Gauss nodes per direction of the kernel
%            integrals.
% OUTPUT Y:  1 x numel(F) cell; Y{t} 1 x M cell, Y{t}{i} size(G{i}.X,1) x
%            P.dim_out(i) double.
[M,N] = size(P.C);
Y = cell(1,numel(F));
for t = 1:numel(F)
    Y{t} = cell(1,M);
    for i = 1:M
        Y{t}{i} = zeros(size(G{i}.X,1),P.dim_out(i));
        for j = 1:N
            b = P.C{i,j};
            if isempty(b),  continue,   end
            [~,ci] = ismember(b.vars.in,F{t}{j}.vars);      % f columns by name
            [~,co] = ismember(b.vars.out,G{i}.vars);
            f = @(Z) eval_f(Z,ci,F{t}{j});
            Y{t}{i} = Y{t}{i} + heatNd_apply(b,f,G{i}.X(:,co),nq);
        end
    end
end
end

function y = eval_f(Z,ci,Fj)
% Y = EVAL_F(Z,CI,FJ) returns the test function Fj at points Z whose column
% k is variable Fj.vars{ci(k)}; a variable of Fj not in CI is zero.
% INPUT  Z:  n x numel(ci) double; ci: double, positions in Fj.vars;
%        Fj: 'struct' with fields vars, deg and coef (TEST_FUNCS).
% OUTPUT y:  n x size(Fj.coef,2) double.
X = zeros(size(Z,1),numel(Fj.vars));    X(:,ci) = Z;
y = prod_terms(X,Fj);
end

function r = rel_res(L,Rt,G)
% R = REL_RES(L,RT,G) returns max_t ||L_t - R_t|| / max_t max(||L_t||,
% ||R_t||), L2 norms on the grids.
% INPUT  L, Rt: outputs of APPLY_C on the grids G; G: grids (OUT_GRID).
% OUTPUT r:     double.
num = 0;    den = 0;
for t = 1:numel(L)
    num = max(num,l2n(cellfun(@minus,L{t},Rt{t},'uni',0),G));
    den = max([den,l2n(L{t},G),l2n(Rt{t},G)]);
end
r = num/max(den,realmin);
end

function Y = neg(Y)
% Y = NEG(Y) returns the output Y of APPLY_C with every value negated.
% INPUT  Y: 1 x nt cell of 1 x M cells of double arrays.
% OUTPUT Y: the same cells, negated.
for t = 1:numel(Y),     Y{t} = cellfun(@uminus,Y{t},'uni',0);  end
end

function n = l2n(Y,G)
% N = L2N(Y,G) returns the L2 norm of the values Y on the grids G, by
% their quadrature weights, summed over the spaces.
% INPUT  Y: 1 x M cell, Y{i} size(G{i}.X,1) x p double; G: grids
%           (OUT_GRID).
% OUTPUT n: double.
n = 0;
for i = 1:numel(Y),     n = n + sum(G{i}.w.*sum(Y{i}.^2,2));    end
n = sqrt(n);
end

function m = min_quad(F,Y,G)
% M = MIN_QUAD(F,Y,G) returns min_t <x_t, P x_t> / (||x_t|| ||P x_t||),
% x_t evaluated on the grids of the same (square) spaces.
% INPUT  F: test functions (TEST_FUNCS); Y: the values of P x_t on G
%           (APPLY_C); G: grids (OUT_GRID) on the output spaces of P, which
%           equal its input spaces.
% OUTPUT m: double.
m = inf;
for t = 1:numel(F)
    % x_t on the grids, <x_t, P x_t> summed over the spaces, then the ratio
    ip = 0;     xs = cell(1,numel(G));
    for i = 1:numel(G)
        [~,ci] = ismember(G{i}.vars,F{t}{i}.vars);
        xs{i} = eval_f(G{i}.X,ci,F{t}{i});
        ip = ip + sum(G{i}.w.*sum(xs{i}.*Y{t}{i},2));
    end
    m = min(m,ip/max(l2n(xs,G)*l2n(Y{t},G),realmin));
end
end


% ========================================================================
% (iv) a solved program with a free scalar and an inequality: H-infinity
% ========================================================================
function out = part_iv()
% OUT = PART_IV() runs part (iv): builds and solves the 1-D primal
% H-infinity LPI on CX_PLANT('io1') with gamma a decision variable (gam >=
% 0, minimized), extracts its operators and checks them, and errors at the
% first failure. Here solinfo.x holds one more entry than decvartable (the
% slack of the 'ineq'), which is the case where x and RRx differ.
% INPUT  none.
% OUTPUT out: 'struct' with fields gam, rel_b, relK, rel1, maxK, rel_x,
%             ndec and nx (R.iv of TEST_GETSOL_SOP).
% Gauss nodes and test functions; the settings; the PIE and its operators
% as containers
nqi = 10;   nqo = 12;   ntest = 6;
st = cx_settings('light','mosek');
PIE = initialize(cx_plant('io1'));
Tm = opvar2copvar(PIE.T);   Am = opvar2copvar(PIE.A);
Bw = opvar2copvar(PIE.Bw);  Cz = opvar2copvar(PIE.Cz);  Dzw = opvar2copvar(PIE.Dzw);
% Identity operators on the spaces of w and z
% [spw,dmw] = cx_space_list(Bw,'in');     [spz,dmz] = cx_space_list(Cz,'out'); % MMP, 10/06/2026 (was)
[spw,dmw] = copvar_space_list(Bw,'in');     [spz,dmz] = copvar_space_list(Cz,'out'); % MMP, 10/06/2026
Iw = eye_copvar_sop(dmw,spw,PIE.dom);   Iz = eye_copvar_sop(dmz,spz,PIE.dom);
% The program; gamma declared, gamma >= 0, gamma the objective
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
gam = dpvar('gam');
prog = lpidecvar(prog,gam);     prog = lpi_ineq(prog,gam);  prog = lpisetobj(prog,gam);
% R positive on the output spaces of T, Q free with the degrees of R, and
% T'*Q - R = 0
% [prog,Rm] = cx_hinf_lf(prog,Tm,PIE,st);                                   % MMP, 10/06/2026 (was)
% [sp,dm] = cx_space_list(Tm,'out');                                        % MMP, 10/06/2026 (was)
% [prog,Qm] = lpivar_cdopvar(prog,dm,sp,PIE.dom,cx_hinf_qdeg(Rm));          % MMP, 10/06/2026 (was)
[prog,Rm] = poslpivar_settings_sop(prog,Tm,st,'lf','out',PIE.dom);          % MMP, 10/06/2026
[sp,dm] = copvar_space_list(Tm,'out');                                      % MMP, 10/06/2026
[prog,Qm] = lpivar_cdopvar(prog,dm,sp,PIE.dom,get_lpivar_degs_sop(Rm));     % MMP, 10/06/2026
prog = lpi_eq_cdopvar(prog,Tm'*Qm - Rm);
% The KYP operator K, its positive slack N and N + K = 0
Km = [-gam*Iw, Dzw', Bw'*Qm;  Dzw, -gam*Iz, Cz;  Qm'*Bw, Cz', Am'*Qm + Qm'*Am];
% prog = cx_hinf_slack(prog,Km,PIE,st);                                     % MMP, 10/06/2026 (was)
[prog,Nm] = poslpivar_settings_sop(prog,Km,st,'slack','out',PIE.dom);       % MMP, 10/06/2026
prog = lpi_eq_cdopvar(prog,Nm + Km,'symmetric');   % Deop + Dop = 0         % MMP, 10/06/2026
% Solve; print the solve and the sizes, then require the certificate
r = cx_solve(prog,st.sos_opts);     sol = r.sol;
ntot = numel(sol.decvartable);
fprintf('(iv) io1 H-inf: pinf %g numerr %g, rel_b %.2e, psd_relmin %.2e; decvartable %d, RRx %d, x %d\n',...
    r.pinf,r.numerr,r.rel_b,r.psd_relmin,ntot,numel(sol.solinfo.RRx),numel(sol.solinfo.x));
assert(r.numerr==0 && r.pinf==0 && r.rel_b<=1e-6 && r.psd_relmin>=-1e-8,'(iv) solve not certified')
% Extract gamma and the operators K, Q, R
g = double(lpigetsol(sol,gam));             % legacy dpvar path
Ks = getsol_lpivar_sop(sol,Km);     Qs = getsol_lpivar_sop(sol,Qm);     Rs = getsol_lpivar_sop(sol,Rm);
assert(isa(Ks,'copvar') && isa(Qs,'copvar') && isa(Rs,'copvar'),'(iv) not fixed')
% Km at the solution, rebuilt from the legacy gamma and the extracted Q.
Kf = [(-g)*Iw, Dzw', Bw'*Qs;  Dzw, (-g)*Iz, Cz;  Qs'*Bw, Cz', Am'*Qs + Qs'*Am];
% The extracted K against the rebuilt K on test functions
F = test_funcs(Ks,ntest);   G = out_grid(Ks,nqo);
YK = apply_c(Ks,F,G,nqi);
relK = rel_res(YK,apply_c(Kf,F,G,nqi),G);
% T'*Q - R = 0 on the extracted operators; KYP operator <= 0.
FT = test_funcs(Rs,ntest);  GT = out_grid(Rs,nqo);
rel1 = rel_res(apply_c(Tm'*Qs,FT,GT,nqi),apply_c(Rs,FT,GT,nqi),GT);
pK = -min_quad(F,neg(YK),G);                % max <x,Kx>/(||x|| ||Kx||)
% Print, then check the tolerances
fprintf('(iv) gamma %.8g | getsol(Km) vs Km(gamma, getsol(Q)): rel %.2e | T''Q-R: rel %.2e | max <x,Kx> %.2e\n',...
    g,relK,rel1,pK);
assert(relK<=1e-10,'(iv) getsol(Km) differs from Km rebuilt with the legacy gamma (%.3g)',relK)
% 1e-5, not (ii)'s 1e-6: this solve stops at rel_b ~5e-7 (measured 09/29),
% and the operator residual measured 8.4e-7 sits at that level.
assert(rel1<=1e-5,'(iv) extracted operators violate T''Q = R (%.3g)',rel1)
assert(pK<=1e-8,'(iv) extracted KYP operator not negative on test functions (%.3g)',pK)
% Control: solinfo.x in RRx's place (x is the solver's cone vector).
bad = sol;  bad.solinfo.RRx = sol.solinfo.x(1:ntot);
relx = rel_res(apply_c(getsol_lpivar_sop(bad,Km),F,G,nqi),apply_c(Kf,F,G,nqi),G);
fprintf('(iv) control solinfo.x: getsol(Km) vs Km(gamma, getsol(Q)) rel %.2e\n',relx);
% Results
out = struct('gam',g,'rel_b',r.rel_b,'relK',relK,'rel1',rel1,'maxK',pK,'rel_x',relx,...
    'ndec',ntot,'nx',numel(sol.solinfo.x));
end


% ========================================================================
% (iii) dispatch against legacy lpigetsol
% ========================================================================
function out = part_iii()
% OUT = PART_III() runs part (iii): the dispatch of LPIGETSOL_SOP,
% GETSOL_LPIVAR_SOP and SUBS_DVAR_SOP over the legacy and the new classes,
% on programs with a known decision vector, and errors at the first
% failure.
% INPUT  none.
% OUTPUT out: 'struct' with fields nlegacy, nnew and nerr (R.iii of
%             TEST_GETSOL_SOP).
% The counts; a 1-D program on [0,1] with a positive dopvar P1, a dopvar
% Q1, a 2 x 1 dpvar g and a cdopvar Pc on R^1 x L2^1[s]
out = struct('nlegacy',0,'nnew',0,'nerr',0);
pvar s theta
prog = lpiprogram(s,theta,[0 1]);
[prog,P1] = poslpivar(prog,[1 1],2);                % dopvar
[prog,Q1] = lpivar(prog,[1 1;1 1],2);               % dopvar
[prog,g] = lpidecvar(prog,[2 1]);                   % dpvar
[prog,Pc] = lpivar_cdopvar(prog,struct('out',[1;1],'in',[1;1]), ...
            struct('out',{{{},{'s'}}},'in',{{{},{'s'}}}),[0 1],2);   % cdopvar in this program
% A known decision vector and a solver report, so the program reads as
% solved; the fixed opvar of P1
N = numel(prog.decvartable);
prog.solinfo.RRx = randn(N+2,1);                    % 2 entries past decvartable, as a slack
prog.solinfo.info = struct('pinf',0,'dinf',0,'numerr',0);
Psol_leg = lpigetsol(prog,P1);                      % a fixed opvar
% 2-D program for dopvar2d / opvar2d.
pvar s1 s2 t1 t2
prog2 = lpiprogram([s1;s2],[t1;t2],[0 1;0 1]);
[prog2,P2] = lpivar(prog2,[1 1;1 1;1 1;1 1],1);     % dopvar2d
prog2.solinfo.RRx = randn(numel(prog2.decvartable),1);
prog2.solinfo.info = struct('pinf',0);
assert(isa(P2,'dopvar2d'),'(iii) lpivar did not return a dopvar2d')
P2sol = lpigetsol(prog2,P2);

% Legacy objects: lpigetsol_sop and subs_dvar_sop against lpigetsol, bit
% for bit
leg = {  {'[]',[]}, {'double',[1 2;3 4]}, {'polynomial',s^2+3*s}, {'opvar',Psol_leg}, ...
         {'dpvar',g}, {'dpvar expr',2*g(1)+g(2)*s}, {'char',g.dvarname{1}}, {'cellstr',g.dvarname(:)'}, ...
         {'cell of dpvar',{g,2*g}}, {'dopvar',P1}, {'dopvar',Q1} };
for k = 1:numel(leg)
    [lbl,Xk] = leg{k}{:};
    a = lpigetsol(prog,Xk);     b = lpigetsol_sop(prog,Xk);
    c = subs_dvar_sop(Xk,prog.decvartable,prog.solinfo.RRx(1:N));
    assert(bit_equal(a,b),'(iii) %s: lpigetsol_sop differs from lpigetsol',lbl)
    assert(bit_equal(a,c),'(iii) %s: subs_dvar_sop differs from lpigetsol',lbl)
    out.nlegacy = out.nlegacy+1;
end
% The comparison must see a different decision vector (non-vacuous).
alt = prog;     alt.solinfo.RRx = prog.solinfo.RRx + 1;
assert(~bit_equal(lpigetsol(prog,P1),lpigetsol_sop(alt,P1)) && ...
       ~bit_equal(lpigetsol(prog,g),lpigetsol_sop(alt,g)),'(iii) bit_equal is blind to the values')
% Same for the 2-D classes, and getsol_lpivar_sop against getsol_lpivar.
leg2 = {{'dopvar2d',P2},{'opvar2d',P2sol}};
for k = 1:numel(leg2)
    [lbl,Xk] = leg2{k}{:};
    a = lpigetsol(prog2,Xk);    b = lpigetsol_sop(prog2,Xk);
    assert(bit_equal(a,b),'(iii) %s: lpigetsol_sop differs from lpigetsol',lbl)
    out.nlegacy = out.nlegacy+1;
end
assert(bit_equal(getsol_lpivar_sop(prog,P1),getsol_lpivar(prog,P1)),'(iii) getsol_lpivar_sop(dopvar) differs')
assert(bit_equal(getsol_lpivar_sop(prog2,P2),getsol_lpivar(prog2,P2)),'(iii) getsol_lpivar_sop(dopvar2d) differs')
% Names keyed, in any order: a permuted list gives the legacy answer.
pm = randperm(N);
c = subs_dvar_sop(g,prog.decvartable(pm),prog.solinfo.RRx(pm));
assert(bit_equal(c,lpigetsol(prog,g)),'(iii) permuted names change the dpvar result')

% A name listed twice in decvartable (a Gram's (i,j),(j,i)) takes the value
% of its LOWEST row, where getequation puts its At row; sosgetsol agrees.
Zc = cellstr(Pc.Zd(:));     nm1 = Zc{1};
[~,ic] = ismember(Zc,prog.decvartable);     dv = prog.solinfo.RRx(ic);
% The first name of Pc appended (dupA) or prepended (dupB) with value 99
dupA = prog;    dupA.decvartable = [prog.decvartable(:); {nm1}];   dupA.solinfo.RRx = [prog.solinfo.RRx(1:N); 99];
dupB = prog;    dupB.decvartable = [{nm1}; prog.decvartable(:)];   dupB.solinfo.RRx = [99; prog.solinfo.RRx(1:N)];
dvB = dv;       dvB(1) = 99;
assert(same_ops(getsol_lpivar_sop(dupA,Pc),subs_dvar_sop(Pc,Zc,dv),0) && ...
       same_ops(getsol_lpivar_sop(dupB,Pc),subs_dvar_sop(Pc,Zc,dvB),0),'(iii) duplicate names: not the lowest row')
assert(double(sosgetsol(dupA,dpvar({nm1})))==dv(1) && double(sosgetsol(dupB,dpvar({nm1})))==99, ...
       '(iii) sosgetsol reads another row than the lowest')

% New classes.
Pcs = lpigetsol_sop(prog,Pc);
assert(isa(Pcs,'copvar') && same_ops(Pcs,getsol_lpivar_sop(prog,Pc),0),'(iii) cdopvar dispatch')
b12 = Pc.C{2,2};
assert(isa(lpigetsol_sop(prog,b12),'sopvar'),'(iii) sdopvar dispatch')
assert(isequal(lpigetsol_sop(prog,Pcs),Pcs),'(iii) copvar not returned as is')
assert(isequal(lpigetsol_sop(prog,Pcs.C{2,2}),Pcs.C{2,2}),'(iii) sopvar not returned as is')
mix = lpigetsol_sop(prog,{Pc,g;Pcs,[]});
assert(isequal(size(mix),[2 2]) && same_ops(mix{1,1},Pcs,0) && bit_equal(mix{1,2},lpigetsol(prog,g)) ...
    && isequal(mix{2,1},Pcs) && isempty(mix{2,2}),'(iii) mixed cell')
% the substitution needs names: none of Pc's are in prog2, so nothing moves
Pnone = getsol_lpivar_sop(prog2,Pc);
assert(isa(Pnone,'cdopvar') && isequal(cellstr(Pnone.Zd(:)),cellstr(Pc.Zd(:))),'(iii) absent names must stay')
out.nnew = 6;

% Errors.
bads = {struct('a',1),true,int8(3),"coeff_1",{1,struct()}};
if exist('nopvar','class')==8,  bads{end+1} = nopvar();     end
for k = 1:numel(bads)
    try
        lpigetsol_sop(prog,bads{k});
        ok = false;
    catch e
        ok = true;
        if ~iscell(bads{k})     % a cell of legacy content is legacy lpigetsol's error
            ok = strcmp(e.identifier,'lpigetsol_sop:badClass');
        end
    end
    assert(ok,'(iii) class %s not rejected',class(bads{k}))
    out.nerr = out.nerr+1;
end
% An unsolved program (empty solver report), and a program without RRx
unsolved = prog;    unsolved.solinfo.info = [];
try
    lpigetsol_sop(unsolved,Pc);     ok = false;
catch
    ok = true;
end
assert(ok,'(iii) unsolved program accepted');
try
    getsol_lpivar_sop(struct('decvartable',{{'a'}}),Pc);    ok = false;
catch e
    ok = strcmp(e.identifier,'getsol_lpivar_sop:notSolved');
end
assert(ok,'(iii) program without RRx accepted');
fprintf('(iii) %d legacy objects bit-identical to lpigetsol, %d new-class checks, %d rejections\n',...
    out.nlegacy,out.nnew,out.nerr);
end


% ========================================================================
% shared helpers
% ========================================================================
function tf = same_ops(A,B,tol)
% TF = SAME_OPS(A,B,TOL) returns true for the same class, the same block
% grid, the same metadata and parameters within tol (absolute, max entry);
% tol = 0 demands bit identity.
% INPUT  A, B: 'copvar', 'cdopvar', 'sopvar' or 'sdopvar' objects; tol:
%              double.
% OUTPUT tf:   logical.
tf = strcmp(class(A),class(B));
if ~tf,     return,     end
% Containers: the grid, the empty blocks, every block, and for a cdopvar
% its decision-variable list
if isa(A,'copvar') || isa(A,'cdopvar')
    tf = isequal(size(A.C),size(B.C));
    for k = 1:numel(A.C)
        if ~tf,     return,     end
        if isempty(A.C{k}) || isempty(B.C{k})
            tf = isempty(A.C{k}) && isempty(B.C{k});
        else
            tf = same_ops(A.C{k},B.C{k},tol);
        end
    end
    if tf && isa(A,'cdopvar'),  tf = isequal(cellstr(A.Zd(:)),cellstr(B.Zd(:)));   end
    return
end
% Blocks: the metadata, then every parameter cell (A and B cells and the
% decision-variable list of an sdopvar)
tf = isequal(A.vars,B.vars) && isequal(A.dom,B.dom) && isequal(A.dims,B.dims) ...
    && isequal(A.ZL,B.ZL) && isequal(A.ZR,B.ZR);
if ~tf,     return,     end
if isa(A,'sopvar')
    for g = 1:numel(A.params)
        tf = tf && isequal(size(A.params{g}),size(B.params{g})) ...
            && maxabs(A.params{g}-B.params{g})<=tol;
    end
else
    for g = 1:numel(A.params.A)
        tf = tf && isequal(size(A.params.A{g}),size(B.params.A{g})) ...
            && maxabs(A.params.A{g}-B.params.A{g})<=tol ...
            && isequal(size(A.params.B{g}),size(B.params.B{g})) ...
            && maxabs(A.params.B{g}-B.params.B{g})<=tol;
    end
    tf = tf && isequal(cellstr(A.Zd(:)),cellstr(B.Zd(:)));
end
end

function tf = bit_equal(a,b)
% TF = BIT_EQUAL(A,B) returns true for bit-identical data, recursively: the
% same class and size, then a 'polynomial' by its coefficient, degmat,
% varname and matdim, another object by each public property, a struct by
% its field names (in any order) and fields, a cell by its entries, and
% any other value by isequal.
% INPUT  a, b: values of any class.
% OUTPUT tf:   logical.
tf = strcmp(class(a),class(b)) && isequal(size(a),size(b));
if ~tf,     return,     end
% 'isequal' cannot be used directly: the polynomial overload is
% ELEMENTWISE (SOSTOOLS400/multipoly/@polynomial/isequal.m) and returns an
% array.
if isa(a,'polynomial')
    tf = isequal(a.coefficient,b.coefficient) && isequal(a.degmat,b.degmat) ...
        && isequal(a.varname,b.varname) && isequal(a.matdim,b.matdim);
elseif isobject(a)
    pr = properties(a);
    for k = 1:numel(pr)
        tf = tf && bit_equal(a.(pr{k}),b.(pr{k}));
    end
elseif isstruct(a)
    fa = fieldnames(a);
    tf = isequal(sort(fa),sort(fieldnames(b)));
    for i = 1:numel(a)
        for k = 1:numel(fa)
            tf = tf && bit_equal(a(i).(fa{k}),b(i).(fa{k}));
        end
    end
elseif iscell(a)
    for k = 1:numel(a),     tf = tf && bit_equal(a{k},b{k});    end
else
    tf = isequal(a,b);
end
end

function m = maxabs(x)
% M = MAXABS(X) returns the largest absolute entry of X as a full double, 0
% for an empty X.
% INPUT  x: numeric array.
% OUTPUT m: double.
m = 0;  if ~isempty(x),     m = full(max(abs(x(:))));   end
end

function p = derangement(n)
% P = DERANGEMENT(N) returns a random permutation of 1:n with no fixed
% point (n >= 2), by rejection; 1:n for n < 2.
% INPUT  n: double.
% OUTPUT p: 1 x n double.
p = 1:n;
if n<2,     return,     end
while any(p==1:n),  p = randperm(n);    end
end

function [x,w] = gauss01(n)
% [X,W] = GAUSS01(N) returns the Gauss-Legendre nodes and weights on [0,1]
% (Golub-Welsch).
% INPUT  n: double, the number of nodes.
% OUTPUT x: n x 1 double, the nodes in increasing order; w: n x 1 double,
%           the weights.
% Nodes: the eigenvalues of the Jacobi matrix of the Legendre recurrence;
% weights: twice the squared first entries of its eigenvectors; then the
% map from [-1,1] to [0,1]
b = (1:n-1)./sqrt(4*(1:n-1).^2-1);
[Vv,L] = eig(diag(b,1)+diag(b,-1));
[x,ix] = sort(diag(L));     w = 2*Vv(1,ix)'.^2;
x = (x+1)/2;    w = w/2;
end
