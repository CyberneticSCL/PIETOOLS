function R = opcheck(prog,H,P,q)
% PROVENANCE.  scratchpad/reach1d/opcheck.m verbatim; gate1d's residual half.
%
% Push a Gram vector q (in NORMALISED-b units) back through the actual PI
% operators and test the LPI directly:
%   Pop = poslpivar(Q1) + eppos*I  (positive whenever Q1 >= 0)
%   Dop = T'PA + A'PT + epneg T'PT ,  Deop = poslpivar(Q2)+poslpivar(Q3)
%   equality to certify:  Dop + Deop = 0
pg = prog;
pg.solinfo.info = struct('verified_by','opcheck');
pg.solinfo.RRx  = full(q(:))*P.nb0;      % back to the original b scaling
Pop  = lpigetsol(pg,H.Pop);
Deop = lpigetsol(pg,H.Deop);
Dop  = H.Top'*Pop*H.Aop + H.Aop'*Pop*H.Top + H.st.epneg*H.Top'*Pop*H.Top;
Res  = Dop + Deop;
R.maxDop = mxop(Dop); R.maxDeop = mxop(Deop); R.maxRes = mxop(Res);
R.rel = R.maxRes/max(R.maxDop,realmin);
R.maxPop = mxop(Pop);
% rel is INVALID wherever the candidate drives Dop -> 0 as an operator: it then
% reads 0/0 and improves for a purely artifactual reason (measured on the
% Example 21 beam, max|Dop| = 2.4e-10 manufactured a clean degree threshold that
% does not exist).  relP divides by max|Pop| instead, which cannot collapse
% because Pop >= eppos2*I by construction.  Accept on relP, never on rel alone.
R.relP = R.maxRes/max(R.maxPop,realmin);
B=numel(P.Ns); R.normQ=zeros(1,B); R.mineig=zeros(1,B); R.rank=zeros(1,B);
for i=1:B
    Q=reshape(q(P.rows{i}),P.Ns(i),P.Ns(i)); Q=(Q+Q')/2; e=sort(eig(Q),'descend');
    R.normQ(i)=norm(Q,'fro')*P.nb0; R.mineig(i)=e(end)*P.nb0;
    R.rank(i)=sum(e>max(e(1),eps)*1e-9);
end
end
function m = mxop(X)
m = max([mxc(X.P) mxc(X.Q1) mxc(X.Q2) mxc(X.R.R0) mxc(X.R.R1) mxc(X.R.R2)]);
end
function m = mxc(X)
if isempty(X), m=0; return; end
if isa(X,'double'), C=X; elseif isa(X,'polynomial'), C=X.coefficient; elseif isa(X,'dpvar'), C=X.C; else, m=NaN; return; end
if isempty(C), m=0; else, m=full(max(abs(C(:)))); end
if isempty(m), m=0; end
end
