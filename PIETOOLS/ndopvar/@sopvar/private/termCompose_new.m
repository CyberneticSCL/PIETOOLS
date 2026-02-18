function [gamma, AB] = termCompose(A,B,alpha,beta,varsin,varsout,varsmid,opAvars,opBvars)
% Input quadPoly parameters have the form:
% A = Z(s4,s2a,s3a)' * CA[alpha] * Z(t2a,t2b,t3a,t3b) : 
% L2[s2a,s2b,s3a,s3b] -> L2[s4,s2a,s3a]

% B = Z(s2a,s2b,s3a,s3b)' * CB[beta] * Z(t1,t3a,t3b)
% L2[s1,s3a,s3b] -> L2[s2a,s2b,s3a,s3b]

% we need to find
% C = A*B
% C = Z(s4,s2a,s3a)'* CC[gamma] * Z(t1,t3b,t3a)
% where we compute CC[gamma] by doing
% C = Z(s4,s2a,s3a)'* CA[alpha] * INT(Z(eta2a,eta2b,eta3a,eta3b)
%                   *Z(eta2a,eta2b,eta3a,eta3b)') *CB[beta] * Z(t1,t3a,t3b)
% L2[s1,s3b,s3a] -> L2[s4,s2a,s3a]

% first, separate out all the necessary metadata about variables
% the below vars are all in s
varsS1 = setdiff(varsin,varsmid);
varsS4 = setdiff(varsout,varsmid);
varsS2a = intersect(setdiff(varsmid,varsin),varsout);
varsS2b = setdiff(setdiff(varsmid,varsin),varsout);
varsS3a = intersect(intersect(varsmid,varsout),varsin);
varsS3b = intersect(setdiff(varsmid,varsout),varsin);

eta2a = strrep(varsS2a,'s','eta');
eta3a = strrep(varsS3a,'s','eta');
eta2b = strrep(varsS2b,'s','eta');
eta3b = strrep(varsS3b,'s','eta');

alph_s3aLoc = find(ismember(opAvars.vars_S3,varsS3a));
beta_s3aLoc = find(ismember(opBvars.vars_S3,varsS3a));

% first, find where the compose terms will be placed
gamma = mapAlphaBetaToGamma(alpha(alph_s3aLoc),beta(beta_s3aLoc));
gammaLinIdx = zeros(size(gamma,1),1);
for i=1:size(gamma,1)
    gammaI = gamma(i,:);
    gammaLinIdx(i) = sub2ind(repmat(3,numel(varsS3a)),gammaI{:});
end

% NOW, we need to find INT(ZA(eta2a,eta2b,eta3a,eta3b)*ZB(eta2a,eta2b,eta3a,eta3b)')
% find which vars are actually present in ZA and ZB
etaMidA = strrep(A.nt,'t','eta'); % changing interior variables in quadPoly to dummies
etaMidB = strrep(B.ns,'s','eta');

% this gives, ZetaA*ZetaB' = S(I\otimes ZetaMid(neta)) = (I\otimesZetaMid'(neta))Sd
[S,Sd,ZetaMid,neta] = monomial_outer_multiple(A.Zt,B.Zs,etaMidA,etaMidB);  

% 2b does not appear in inputs of B or outputs of A, so must be fully
% integrated out
% this gives int_{0}^1 ZetaMid(eta2a,eta2b,eta3a,eta3b) d eta2b = diag(Cint)*ZetaMid2b(eta2a,eta3a,eta3b)
[Cint, ZetaMid2b] = int_definite(ZetaMid,neta,eta2b); 


% we need, int_[2a,2b,3a,3b] ZetaA*ZetaB' = int_[2a,2b,3a,3b] S(I\otimes ZetaMid(neta))
%          = S(I \otimes int_[2a,3a,3b] diag(Cint)*ZetaMid2b(eta2a,eta3a,eta3b))
%          = S(I\otimes diag(Cint)(I\otimes int_[2a,3a,3b] ZetaMid2b(eta2a,eta3a,eta3b)) 
barCA = A.C*S*kron(speye(prod(cellfun(@numel,A.Zt))),spdiags(Cint));


% now, we need to handle int_[2a,3a,3b] ZetaMid2b(eta2a,eta3a,eta3b)
% then can either be, subs eta with s, subs eta with t, multiplier
%                     integrate then, subs eta with s, t, 0, 1



% but ZA and ZB are tensor products. Separate 2a,2b,3a,3b
[SA,ZA2a,ZA2b,ZA3a,ZA3b] = monomial_split(A.Zt,etaMidA,eta2a,eta2b,eta3a,eta3b);  %ZA = SA(ZA2a\otimes ZA2b\otimes ZA3a\otimes ZA3b)
[SB,ZB2a,ZB2b,ZB3a,ZB3b] = monomial_split(B.Zt,etaMidB,eta2a,eta2b,eta3a,eta3b);  %ZB' = (ZB2a\otimes ZB2b\otimes ZB3a\otimes ZB3b)'SB

% now we are at, CA[alpha]*SA*INT((ZA2a\otimes ZA2b\otimes ZA3a\otimes ZA3b)*(ZB2a\otimes ZB2b\otimes ZB3a\otimes ZB3b)')*SB*CB[beta]
% which is same as
% barCA[alpha]*(INT(ZA2a*ZB2a')\otimes INT(ZA2b*ZB2b')\otimes INT(ZA3a*ZB3a')\otimes INT(ZA3b*ZB3b'))*barCB[beta]
% barCA[alpha]*(INT(ZAB2a)\otimes INT(ZAB2b)\otimes INT(ZAB3a)\otimes INT(ZAB3b))*barCB[beta]

barCA = A.C*SA; barCB = SB*B.C;

% next, integrate out eta_2b, because eta_2b appears neither in 
% input nor in output vars
S2b = int(ZA2b,ZB2b,eta2b,lims);

% next, integrate and pull eta_2a to the left, because eta2a only appears
% as s2a at the end
[~,S2a,ZAB2a,sAB2a] = monomial_outer_multiple(ZA2a,ZB2a,etaMidA,etaMidB);   % ZA2a*ZB2a' = (I\otimes ZAB2a)'*S2a
[Cint2a, ZAB2aint] = int_indefinite(ZAB2a); % INT(ZA2a*ZB2a') = (I\otimes ZAB2aint'(sAB2a))*(I\otimes diag(Cint2a))*S2a

% next, integrate and pull eta_3b to the right, because eta3b only appears
% as t3b at the end
[S3b,~,ZAB3b,sAB3b] = monomial_outer_multiple(ZA3b,ZB3b,etaMidA,etaMidB);   % ZA3b*ZB3b' = S3b*(I\otimes ZAB3b)
[Cint3b, ZAB3bint] = int_indefinite(ZAB3b); % INT(ZA3b*ZB3b') = S3b*(I\otimes diag(Cint3b))*(I\otimes ZAB3bint(sAB3b))

% all that is left is the 3a integration


end
