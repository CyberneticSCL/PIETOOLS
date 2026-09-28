function prob = Sedumi2Mosek(A,b,c,K)
% Sedumi2Mosek --- Convert Sedumi inputs to the MOSEK prob format.
%
% prob = Sedumi2Mosek(A,b,c,K)
%
% Description: The Sedumi SDP is defined as
% min c^T x :
% Ax=b
% x_i \in K_i
% The elements of K are: K.f=free, K.l=PO, K.q=Lor, K.r=R Lor, K.s=SDP
% We will assume only K.f, K.l and K.s are non-empty. quadratic constraints
% will return an error.
%
% The mosek SDP format is
% min c^Tx + \sum_i <C_i,X_i>
% lbi \le a^Tx + \sum_j <A_(ij),X_j> \le ubi,
%
% which is all contained in the prob structure.
%
% Inputs:
% 
% A: An m x n matrix
% b: A length m column vector
% c: A length n column vector
% K: is the cone structure from the Sedumi input format.
%
% prob: is the input structure accepted by MOSEK
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2021  M. Peet, S. Shivakumar, D. Jagt
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% 04/24/14 - MP - initial coding

% A Sedumi to Mosek converter
% changed sparse construction for faster conversion, SS - 7/20/2021
% 09/26/26 - MMP - bara (the constraint matrices A_ij of Mosek's form) is
%   built per psd block in one pass over that block's nonzeros, instead of
%   one constraint row at a time. The old loop took row i of A on block j,
%   reshaped it with mat() and ran find(tril(X+X')/2), for every row i and
%   every block j. Extracting a row of a sparse (compressed-column) A
%   searches all K.s(j)^2 columns of the block however few nonzeros the row
%   has, so the conversion cost O(m*sum(K.s.^2)), i.e. (equality
%   constraints) x (Gram-matrix entries), plus one cellfun call and one
%   sparse K.s(j) x K.s(j) matrix per row and block. It now costs
%   O(nnz(A) + size(A,2) + m*numel(K.s)). Measured (min of reps, one
%   process): sosdemo11 (m = 5555, K.s = [100 100 150]) 0.86 s -> 0.001 s;
%   sosdemo12 0.73 s -> 0.0015 s; an SDP with m = 13977, K.s = [8 869],
%   nnz(A) = 2.3e6: 72 s -> 0.08 s, peak memory 314 -> 191 MB (86 MB of it
%   prob.bara itself); there the old conversion was 19-24% of the whole
%   Mosek solve.
%   The per-block pieces are joined once after the loop: appending them to
%   the full arrays block by block, as the old code did, costs
%   O(numel(K.s)*nnz(A)); with the new extraction, 37 s against 1.4 s for
%   10^4 blocks of size 3 and m = 10^5.
%   The output is unchanged bit for bit, every prob field, including the
%   order, orientation and empty shapes of the bara triplets: the same
%   operations are applied to the same operands (X+X' with the same
%   operand order, so NaN bits, cancellations and overflow agree; /2;
%   find).
%   Replaces the cellfun construction of SS - 7/20/2021 (kept commented
%   out below).

if (isfield(K,'q') && any(K.q>0)) || (isfield(K,'r') && any(K.r>0))
    error('K.q and K.r are supposed to be empty')
else
    if isfield(K,'q')
        K.q = rmfield(K,'q');
    end
    if isfield(K,'r')
        K.r = rmfield(K,'r');
    end
end

if ~isfield(K,'f')
    K.f=[];
end
if ~isfield(K,'l')
    K.l=[];
end
nf=sum(K.f); nl=sum(K.l);

Kindex = nf+nl+cumsum([0 (K.s).^2])+1; % vector with the starting indices of the SDP variables
n=size(A,2);
m=size(A,1);

nsdpvar=length(K.s);


lpA=A(:,1:(Kindex(1)-1));
lpc=c(1:(Kindex(1)-1));

% The structure is stored in prob
prob.barc.subj=[];
prob.barc.subk=[];
prob.barc.subl=[];
prob.barc.val=[];

prob.bara.subi=[];
prob.bara.subj=[];
prob.bara.subk=[];
prob.bara.subl=[];
prob.bara.val=[];

prob2.bara.subi=[];
prob2.bara.subj=[];
prob2.bara.subk=[];
prob2.bara.subl=[];
prob2.bara.val=[];
% bara pieces per block, joined once after the loop: appending to the       % MMP, 09/26/2026
% full arrays in every block cost O(numel(K.s)*nnz(A)).                     % MMP, 09/26/2026
Bi = cell(1,nsdpvar);  Bj = Bi;  Bk = Bi;  Bl = Bi;  Bv = Bi;               % MMP, 09/26/2026
for j=1:nsdpvar
    temp=mat(c(Kindex(j):(Kindex(j+1)-1)),K.s(j));
    [I,J,V] = find(tril(temp+temp')/2);
    prob.barc.subj=[prob.barc.subj j*ones(1,length(I))];
    prob.barc.subk=[prob.barc.subk I' ] ;
    prob.barc.subl=[prob.barc.subl J' ];
    prob.barc.val=[prob.barc.val V' ];
    clear I J V
    


% The barA matrices are formed from the rows of A
%     for i=1:m
%         barA= mat(A(i,Kindex(j):(Kindex(j+1)-1)),K.s(j));% j references the variable number (K.s(j)) and i references the ith row of A
%         [I,J,V] = find(tril(barA+barA')/2);
%         prob.bara.subi=[prob.bara.subi i*ones(1,length(I))];
%         prob.bara.subj=[prob.bara.subj j*ones(1,length(I))];
%         prob.bara.subk=[prob.bara.subk I' ];
%         prob.bara.subl=[prob.bara.subl J' ];
%         prob.bara.val=[prob.bara.val V' ];
%         clear I J V
%     end

%     cell implementation
    
% BEGIN MMP, 09/26/2026: block j in one pass over its nonzeros, not one
% cellfun call per constraint row (see header).
%   barAcell = num2cell((1:m)');                                            % MMP, 09/26/2026 (was)
%   barAcell = cellfun(@(x) mat(A(x,Kindex(j):(Kindex(j+1)-1)),K.s(j)), barAcell, 'un',0); % MMP, 09/26/2026 (was)
%   [I, J, V] = cellfun(@(x) find(tril(x+x')/2), barAcell, 'un', 0);        % MMP, 09/26/2026 (was)
%   idxI = num2cell((1:m)');                                                % MMP, 09/26/2026 (was)
%   prob2.bara.subi = [prob2.bara.subi (cell2mat(cellfun(@(x,y) y*ones(length(x),1), I, idxI,'un',0)))']; % MMP, 09/26/2026 (was)
%   prob2.bara.subj = [prob2.bara.subj (cell2mat(cellfun(@(x) j*ones(length(x),1), I,'un',0)))']; % MMP, 09/26/2026 (was)
%   prob2.bara.subk = [prob2.bara.subk (cell2mat(I))'];                     % MMP, 09/26/2026 (was)
%   prob2.bara.subl = [prob2.bara.subl (cell2mat(J))'];                     % MMP, 09/26/2026 (was)
%   prob2.bara.val = [prob2.bara.val (cell2mat(V))'];                       % MMP, 09/26/2026 (was)
    % Row i of A on block j is vec(X_i)'; bara holds the lower triangle of  % MMP, 09/26/2026
    % (X_i+X_i')/2, i outer, column-major within X_i. All rows at once, by  % MMP, 09/26/2026
    % the old operations: the same + with the same operand order (so NaN    % MMP, 09/26/2026
    % bits, cancellation, overflow agree), /2 (drops underflow), find.      % MMP, 09/26/2026
    if m==0,  continue,  end        % old appended 0x0 (no rows)            % MMP, 09/26/2026
    s = K.s(j);                                                             % MMP, 09/26/2026
    [kk,ll] = find(tril(true(s)));  % lower triangle k>=l, column-major     % MMP, 09/26/2026
    % Block column (l-1)*s+k is X_i(k,l); Y(i,:) = X_i(k,l)+X_i'(k,l).      % MMP, 09/26/2026
    Y = A(:,Kindex(j)-1+(ll-1)*s+kk) + A(:,Kindex(j)-1+(kk-1)*s+ll);        % MMP, 09/26/2026
    Y = (Y/2).';                    % transposed: find runs i outer, as old % MMP, 09/26/2026
    [p,ii,vv] = find(Y);                                                    % MMP, 09/26/2026
    clear Y                                                                 % MMP, 09/26/2026
    kk = kk(p);  ll = ll(p);                                                % MMP, 09/26/2026
    kk = kk(:)'; ll = ll(:)'; ii = ii(:)'; vv = vv(:)';  % rows, as old     % MMP, 09/26/2026
    Bi{j} = ii;                                                             % MMP, 09/26/2026
    Bj{j} = j*ones(1,numel(ii));                                            % MMP, 09/26/2026
    % Old shape for an empty block with s<=1: find of a 1x1 zero is 0x0,    % MMP, 09/26/2026
    % so cell2mat(I) was 0x0 there (1x0 for s>1, and for subi/subj);        % MMP, 09/26/2026
    % a 0x0 piece is left empty, which the join below skips.                % MMP, 09/26/2026
    if ~isempty(vv) || s>1                                                  % MMP, 09/26/2026
        Bk{j} = kk;                                                         % MMP, 09/26/2026
        Bl{j} = ll;                                                         % MMP, 09/26/2026
        Bv{j} = vv;                                                         % MMP, 09/26/2026
    end                                                                     % MMP, 09/26/2026
% END MMP, 09/26/2026
end
% One horzcat equals the per-block appends it replaces: 0x0 pieces are      % MMP, 09/26/2026
% skipped, 1x0 ones kept, so shapes of empty results are unchanged.         % MMP, 09/26/2026
prob2.bara.subi = [prob2.bara.subi Bi{:}];                                  % MMP, 09/26/2026
prob2.bara.subj = [prob2.bara.subj Bj{:}];                                  % MMP, 09/26/2026
prob2.bara.subk = [prob2.bara.subk Bk{:}];                                  % MMP, 09/26/2026
prob2.bara.subl = [prob2.bara.subl Bl{:}];                                  % MMP, 09/26/2026
prob2.bara.val = [prob2.bara.val Bv{:}];                                    % MMP, 09/26/2026

prob.bara = prob2.bara;

prob.a=sparse(lpA);
prob.c=sparse(lpc);
prob.bardim=K.s;

prob.blc=b; %
prob.buc=b;

idxI = 1:nf; idxJ = ones(nf,1); vals = -inf*ones(nf,1);
M = nf+nl; N = 1;
prob.blx = sparse(idxI,idxJ,vals,M,N);
% prob.blx=sparse([-inf*ones(nf,1);zeros(nl,1)]);
prob.bux=[];

