function [gammaLinIdx, CparamCells] = termCompose(KA, KB, betaIdx, alphaIdx, s2a, s2b, s3a, s3b)
% first pull out K^A and K^B components (coeffs and monomials)

interiorVars = union(s2a,s2b,s3a,s3b);

ZKA_left = KA.Zs;
ZKA_right = KA.Zt;
CKA = KA.C;

ZKB_left = KB.Zs;
ZKB_right = KB.Zt;
CKB = KB.C;

s2aloc = ismember(KB.ns,s2a);
s3bloc = ismember(KA.nt,strrep(s3b,'s','t'));
% get alpha_a, alpha_b, beta_a, beta_b
beta_a = betaIdx(s2aloc); beta_b = betaIdx(~s2aloc);
alpha_a = alphaIdx(~s3bloc); alpha_b = alphaIdx(s3bloc);

%now, let us start to compute K^Z(s2a,s3a|s3a',s3b)
% K^Z(s2a,s3a|s3a',s3b) = int(  int(ZA_right*ZB_left', s2b),   s2a', s3a'', s3b')
% first multiply inner monomials and integrate out s2b, 
% i.e., int(ZA_right*ZB_left', s2b) = (I\otimes Zs2a_p)*Gs3a_pp*(I\otimes Zs3b_p)
[Zs2a_p,Gs3a_pp,Zs3b_p] = int_2b(ZKA_right, ZKB_left, interiorVars, s2a, s2b, s3a, s3b);


% next do the outer integration since variables are separated
% K^Z(s2a,s3a|s3a',s3b) = 
%          =int((I\otimes Zs2a_p), s2a')*int(Gs3a_pp,s3a'')*int(I\otimes Zs3b_p,s3b')
%          =(I\otimes KZhat_left')*KG*(I\otimes KZhat_right)
KZhat_left = int_2a(Zs2a_p,beta_a);
KZhat_right = int_3b(Zs3b_p, alpha_b);
KG = int_3a(Gs3a_pp,beta_b,alpha_a);

%KG has structure KG.gamma, KG.Cgamma, KG.Zleft3a, KG.Zright3a



% we will treat KG last. First KZhat_left,
% (I\otimes ZKA_left')*CKA*(I\otimes KZhat_left') = (I\otimes KC_left')*CL
[KC_left, CL] = leftshiftMonomoials(ZKA_left, CKA, KZhat_left, KA.ns, s2a); 
% likewise, right shift 3b,
% (I\otimes KZhat_right)*CKB*(I\otimes ZKB_right) = CR*(I\otimes KC_right)
[KC_right, CR] = rightshiftMonomials(ZKB_right, CKB, KZhat_right, KB.nt, s3b);
% maybe we can compute CL and CR only once outside???

% next, we have CL*KG*CR for different gamma in KG
% CL*KG(gamma)*CR = (I\otimes Zbar_3a^T)*CG(gamma)*(I\otimes Zbar_3ap)
[Zbar_3a,Zbar_3ap,CG] = leftrightshiftMonomials(CL, KG, CR);

% lastly, we have 
% (I\otimes KC_left')*(I\otimes Zbar_3a^T)*CG(gamma)*(I\otimes Zbar_3ap)*(I\otimes KC_right)
%  = (I\otimes ZC_left)*CC_gamma*(I\otimes ZC_right)
[ZC_left, CC, ZC_right] = combineMonomials(KC_left, KC_right, Zbar_3a,Zbar_3ap,CG);

gammaLinIdx = cell(1,length(CC.gamma));
CparamCells = cell(1,length(CC.gamma));
for i=1:length(CC.gamma)
    gammaLinIdx{i} = CC.gamma;
    CparamCells{i} = quadPoly(CC.C, ZC_left, ZC_right,[KA.dim(1),KB.dim(2)],KA.ns,KB.nt);
end
end

function [Zs2a_p,Gs3a_pp,Zs3b_p] = int_2b(ZKA_right, ZKB_left, interiorVars, s2a, s2b, s3a, s3b)
end
function KZhat_left = int_2a(Zs2a_p,beta_a)
end
function KZhat_right = int_3b(Zs3b_p, alpha_b)
end
function KG = int_3a(Gs3a_pp,beta_b,alpha_a)
end
function [KC_left, CL] = leftshiftMonomoials(ZKA_left, CKA, KZhat_left, KA.ns, s2a)
end
function [KC_right, CR] = rightshiftMonomials(ZKB_right, CKB, KZhat_right, KB.nt, s3b)
end
function [Zbar_3a,Zbar_3ap,CG] = leftrightshiftMonomials(CL, KG, CR)
end
function [ZC_left, CC, ZC_right] = combineMonomials(KC_left, KC_right, Zbar_3a,Zbar_3ap,CG)
end