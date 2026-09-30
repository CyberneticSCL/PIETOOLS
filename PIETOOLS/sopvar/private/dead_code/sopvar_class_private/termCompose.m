function [gammaLinIdx, CparamCells] = termCompose(KA, KB, betaIdx, alphaIdx, s2a, s2b, s3a, s3b)
% first pull out K^A and K^B components (coeffs and monomials)

interiorVars = union(strrep(KA.nt,'t','s'),KB.ns);

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
[Zs2a_p,Gs3a_pp,Zs3b_p] = int_2b(ZKA_right, ZKB_left, strrep(KA.nt,'t','s'), ZKB.ns, s2a, s2b, s3a, s3b);


% next do the outer integration since variables are separated
% K^Z(s2a,s3a|s3a',s3b) = 
%          =int((I\otimes Zs2a_p), s2a')*int(Gs3a_pp,s3a'')*int(I\otimes Zs3b_p,s3b')
%          =(I\otimes KZhat_left')*KG*(I\otimes KZhat_right)
KZhat_left = int_onesided(Zs2a_p,beta_a);  % limits are either delta, 0-to-si, or si-to-1
KZhat_right = int_onesided(Zs3b_p, alpha_b); % limits are either delta, 0-to-si, or si-to-1
KG = int_twosided(Gs3a_pp,beta_b,alpha_a); % limits are in both s3a and s3a'
%KG has structure KG.gamma, KG.Cgamma, KG.Zleft3a, KG.Zright3a
%KZhat_left and _right have structure fields .C, .Z



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


function KZhat = int_onesided(Zs, beta)
n = numel(Zs);
KZhat.Z = cell(1,n);
Ci = cell(1,n);
for i = 1:n
    d = Zs{i}(:);    
    m = numel(d);
    if beta(i) == 1      % no integration
        KZhat.Z{i} = d;
        Ci{i} = speye(m);
    elseif beta(i) == 2  % int_0^{s_i}
        KZhat.Z{i} = d + 1;
        coeff = 1 ./ (d + 1);
        Ci{i} = spdiags(coeff, 0, m, m);   % m x m
    elseif beta(i) == 3  % int_{s_i}^{1}
        KZhat.Z{i} = [0; d + 1];      % length m+1
        coeff = 1 ./ (d + 1);
        % Ci is m x (m+1): first column = +coeff (constant term),
        % remaining columns = -diag(coeff) multiplying s^(d+1)
        Ci{i} = [sparse(coeff), -spdiags(coeff, 0, m, m)];
    else
        error('beta_a(%d) must be 1,2,or 3.', i);
    end
end
KZhat.C = 1;
for i = 1:n
    KZhat.C = kron(KZhat.C, Ci{i});
end
end
function KG = int_twosided(Gs3a_pp,beta_b,alpha_a)
end
function [KC_left, CL] = leftshiftMonomoials(ZKA_left, CKA, KZhat_left, KA.ns, s2a)
end
function [KC_right, CR] = rightshiftMonomials(ZKB_right, CKB, KZhat_right, KB.nt, s3b)
end
function [Zbar_3a,Zbar_3ap,CG] = leftrightshiftMonomials(CL, KG, CR)
end
function [ZC_left, CC, ZC_right] = combineMonomials(KC_left, KC_right, Zbar_3a,Zbar_3ap,CG)
end


function delta_vec = p_multi(gamma, alpha, beta)
    n = numel(gamma);
    delta_vec = zeros(1,n);
    for i = 1:n
        delta_vec(i) = p_scalar(gamma(i), alpha(i), beta(i));
    end
end

function d = p_scalar(g,a,b)
    % p(0,alpha,beta)
    if g == 0
        if a==0 && b==0
            d = 1;
        else
            d = 0;
        end
        return;
    end
    
    % p(1,alpha,beta)  
    if g == 1
        if (b==0 && a==0),  d=0; return; end
        if (b==0 && a==1),  d=1; return; end
        if (b==0 && a==-1), d=0; return; end
        if (b==1 && a==0),  d=2; return; end
        if (b==1 && a==1),  d=7; return; end
        if (b==1 && a==-1), d=4; return; end
        if (b==-1 && a==0), d=0; return; end
        if (b==-1 && a==1), d=5; return; end
        if (b==-1 && a==-1),d=0; return; end
    end
    
    % p(-1,alpha,beta)
    if g == -1
        if (b==0 && a==0),  d=0; return; end
        if (b==0 && a==1),  d=0; return; end
        if (b==0 && a==-1), d=1; return; end
        if (b==1 && a==0),  d=0; return; end
        if (b==1 && a==1),  d=0; return; end
        if (b==1 && a==-1), d=3; return; end
        if (b==-1 && a==0), d=2; return; end
        if (b==-1 && a==1), d=6; return; end
        if (b==-1 && a==-1),d=8; return; end
    end
end