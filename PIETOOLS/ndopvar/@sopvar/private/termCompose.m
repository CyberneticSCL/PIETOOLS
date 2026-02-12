function [gamma, AB] = termCompose(A,B,alpha,beta,varsin,varsout,varsmid)
% Input quadPoly parameters have the form:
% A = Z(s4,s2a,s3a)' * CA[alpha] * Z(s2a',s2b,s3a',s3b) : 
% L2[s2a,s2b,s3a,s3b] -> L2[s4,s2a,s3a)

% B = Z(s2a,s2b,s3a,s3b)' *CB[beta] * Z(s1,s3a',s3b')
% L2[s1,s3a,s3b] -> L2[s2a,s2b,s3a,s3b]

% we need to find
% C = A*B
% C = Z(s4,s2a,s3a)'* CA[alpha] * Z(s2a',s2b,s3a',s3b)
%                   *Z(s2a',s2b,s3a',s3b)' *CB[beta] * Z(s1,s3a'',s3b')

% first find where the compose terms will be placed
gamma = mapAlphaBetaToGamma(alpha,beta);

% now compute the composition of A, B quadPolys
% interior degmats
ZetaA = A.Zt;
ZetaB = B.Zs;

% compute S(I\otimes ZetaMid) = (I\otimes ZetaMid')Sd = ZetaA*ZetaB'
[S,Sd,ZetaMid] = monomial_outer_multiple(ZetaA,ZetaB,A.nt,B.ns);

% next do integration: \int_[s2a',s2b,s3a',s3b] ZetaMid 

end