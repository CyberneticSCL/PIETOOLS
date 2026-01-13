function C = plus(A,B)
Cparams = A.params;

Bparams = B.params;
kB = keys(Bparams);

hit = isKey(Cparams, kB);
kBoth = kB(hit);
kNew  = kB(~hit);


if ~isempty(kBoth)
vC = Cparams(kBoth);
vB = Bparams(kBoth);

vSum = vC;
% Add where keys overlap; this can be done in parallel, parfor loops?
for i = 1:numel(kBoth)
    vSum{i} = vSum{i} + vB{i};
end
Cparams(kBoth) = vSum;
end

% Insert new keys
if ~isempty(kNew)
    Cparams{kNew} = Bparams{kNew};
end

C = nDopvar(A.vars_in, A.vars_out, A.dims, Cparams);
end


% naive implementation below

% function C = plus(A,B)
% Aparams = A.params;
% Bparams = B.params;
% 
% keyA = keys(Aparams)';
% keyB = keys(Bparams)';
% keyC = union(keyA,keyB);
% 
% inA = isKey(Aparams, keyC);
% inB = isKey(Bparams, keyC);
% 
% paramsC = cell(1,length(keyC));
% 
% for i=1:length(keyC)
%     if ~ismember(keyC(i),keyA) 
%         paramsC{i} = Bparams{keyC(i)};
%     elseif ~ismember(keyC(i),keyB) 
%         paramsC{i} = Aparams{keyC(i)};
%     else
%         paramsC{i} = Aparams{keyC(i)}+Bparams{keyC(i)}; 
%     end
% end
% 
% C = nDopvar(A.vars_in, A.vars_out, A.dims, dictionary(keyC, paramsC));
% end