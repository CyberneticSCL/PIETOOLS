function C = mtimes(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = mtimes(A,B) composes two sopvar operators C=AB where the output
% domain of B must be the input domain of A
% 
% Date: 1/16/26
% Version: 1.0
% 
% INPUT
% A: sopvar class object L_2^q[Sj,...Sn] to L_2^r[Sk,...,Sm]
% B: sopvar class object L_2^p[Si,...Sn] to L_2^q[Sj,...Sn]
%
% OUTPUT
% C = AB:  L_2^p[Si,...,Sm] to L_2^r[Sk,...Sn]
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - mtimes
%
% Copyright (C)2026  M. Peet, S. Shivakumar
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
% Initial coding MMP, SS  - 1_16_2026
%
% Error handling: Checks to ensure A and B are compatible
if A.dims(1)~=B.dims(2)
    error('number of output components of B is different than number of input components of A');
end
if any(~strcmp(A.vars_in,B.vars_out))
    error('names of B output variables differ from names of A input variables');
end
if any(A.dom_out~=B.dom_in)
    error('output domain of B differs from input domain of A');
end

if numel(A.params)~=numel(B.params)
    error('number of terms in summands is not equal -- one of them is probably malformed');
end

C=A;                % begin with A -- creates correct output vars, dim, dom
C.dims(1)=B.dims(1); 
C.vars_in=B.vars_in;  
C.dom_in=B.dom_in;

% Now, just need to compute the parameters!

for i=1:N
    for j=1:N
        for k=1:N
            quadcompose(A.params{i,j,k},A.params{i,j,k})
        end
    end
end



Aparams = A.params;
Bparams = B.params;
keyA = keys(Aparams);
keyB = keys(Bparams);

valA = Aparams(keyA);   % cell array of polynomials
valB = Bparams(keyB);   % cell array of polynomials

Avars = vars(A);
Bvars = vars(B);
midvars = A.vars_in;

nA = numel(Avars);
nB = numel(Bvars);
nmid = numel(midvars);
for i=1:nmid
Var(i) = pvar(midvars{i});
midVar(i) = pvar([midvars{i},'_mid']);
dumVar(i) = pvar([midvars{i},'_dum']);
end

Cvars = union(A.vars_out,B.vars_in);
nC = numel(Cvars);

maxkeyC = 4^(nC);

keyC = 1:maxkeyC;
paramsC = cell(1,maxkeyC);

keyAidx = nDopvar.keys2index(keyA,nA);
keyBidx = nDopvar.keys2index(keyB,nB);

% Position maps (all computed once)
[~, posAinC] = ismember(Avars, Cvars);
[~, posBinC] = ismember(Bvars, Cvars);

[~, posMidA] = ismember(midvars, Avars);
[~, posMidB] = ismember(midvars, Bvars);
[~, posMidC] = ismember(midvars, Cvars);




for i=1:numel(keyA)
    for j=1:numel(keyB)
        [keysList, paramsList] = ...
            compose(valA{i}, valB{j}, Avars, Bvars, Cvars, A.vars_in, keyAidx(i,:), keyBidx(j,:), nC, nmid,...
            posAinC,posBinC, posMidA,posMidB,posMidC, ...
            Var,midVar,dumVar);
        for k=1:numel(keysList)
            paramsC{keysList(k)} = paramsC{keysList(k)}+paramsList{k};
        end
    end
end

% prune empty parameters
emptyIdx = cellfun('isempty', paramsC);
keyC = keyC(~emptyIdx);
paramsC = paramsC(~emptyIdx);

C = nDopvar(B.vars_in, A.vars_out, [A.dims(1), B.dims(2)], dictionary(keyC,paramsC));
end

function [keys,params] = compose(Aparam, Bparam, Avars, Bvars, Cvars, midvars,...
                            keyAidx, keyBidx, nC, nmid,...
                            posAinC,posBinC,posMidA,posMidB,posMidC, ...
                            Var,midVar,dumVar)
keyCidx(:,posAinC) = keyAidx;
keyCidx(:,posBinC) = keyBidx;

keys = zeros(1,2^nmid);
params = cell(1,2^nmid);

l = 1;



for i=1:nmid
    keyinA = keyAidx(posMidA(i));
    keyinB = keyBidx(posMidB(i));
    keyinC = posMidC(i);
    var = Var(i);
    dumvar = dumVar(i);
    midvar = midVar(i);
    if keyinA==0
        if keyinB==0
            keyCidx(keyinC) = 0;
            keys(l) = nDopvar.index2keys(keyCidx,nC);
            params{l} = Aparam*Bparam;
            l=l+1;
        elseif keyinB==1
            keyCidx(keyinC) = 1;
            keys(l) = nDopvar.index2keys(keyCidx,nC);
            params{l} = Aparam*Bparam;
            l=l+1;
        elseif keyinB==2
            keyCidx(keyinC) = 2;
            keys(l) = nDopvar.index2keys(keyCidx,nC);
            params{l} = Aparam*Bparam;
            l=l+1;
        end
    elseif keyinA==1
        if keyinB==0
            keyCidx(keyinC) = 1;
            keys(l) = nDopvar.index2keys(keyCidx,nC);
            params{l} = Aparam*subs(Bparam,var,dumvar);
            l=l+1;
        elseif keyinB==1
            keyCidx(keyinC) = 1;
            keys(l) = nDopvar.index2keys(keyCidx,nC);
            params{l} = int(subs(Aparam,dumvar,midvar)*subs(Bparam,var,midvar),midvar,dumvar,var);
            l=l+1;
        elseif keyinB==2
            ABparam = subs(Aparam,dumvar,midvar)*subs(Bparam,var,midvar);

            keyCidx(keyinC) = 1;
            keys(l) = nDopvar.index2keys(keyCidx,nC);
            params{l} = int(ABparam, midvar, 0, dumvar);
            l=l+1;

            keyCidx(keyinC) = 2;
            keys(l) = nDopvar.index2keys(keyCidx,nC);
            params{l} = int(ABparam, midvar, 0, var);
            l=l+1;
        end
    elseif keyinA==2
        if keyinB==0
            keyCidx(keyinC) = 2;
            keys(l) = nDopvar.index2keys(keyCidx,nC);
            params{l} = Aparam*subs(Bparam,var,dumvar);
            l=l+1;
        elseif keyinB==1
            ABparam = subs(Aparam,dumvar,midvar)*subs(Bparam,var,midvar);

            keyCidx(keyinC) = 1;
            keys(l) = nDopvar.index2keys(keyCidx,nC);
            params{l} = int(ABparam,midvar, var, 1);
            l=l+1;
            
            keyCidx(keyinC) = 2;
            keys(l) = nDopvar.index2keys(keyCidx,nC);
            params{l} = int(ABparam,midvar, dumvar, 1);
            l=l+1;
        elseif keyinB==2
            keyCidx(keyinC) = 2;
            keys(l) = nDopvar.index2keys(keyCidx,nC);
            params{l} = int(subs(Aparam,dumvar,midvar)*subs(Bparam,var,midvar),midvar,var,dumvar);
            l=l+1;
        end
    elseif keyinA==3
        if keyinB==0
            keyCidx(keyinC) = 3;
            keys(l) = nDopvar.index2keys(keyCidx,nC);
            params{l} = Aparam*Bparam;
            l=l+1;
        elseif keyinB==1
            keyCidx(keyinC) = 3;
            keys(l) = nDopvar.index2keys(keyCidx,nC);
            params{l} = int(subs(Aparam*Bparam,var,midvar),midvar,var,1);
            l=l+1;
        elseif keyinB==2
            keyCidx(keyinC) = 3;
            keys(l) = nDopvar.index2keys(keyCidx,nC);
            params{l} = int(subs(Aparam*Bparam,var,midvar),midvar,0,var);
            l=l+1;
        end
    end
end

% truncate
keys   = keys(1:l-1);
params = params(1:l-1);

% dedupe: group identical keys and sum corresponding params
% if numel(keys) > 1
%     [ks, ord] = sort(keys(:));     % column
%     ps = params(ord);
% 
%     isNew = [true; diff(ks) ~= 0];
%     uKeys = ks(isNew);
% 
%     uParams = cell(size(uKeys));
%     t = 1;
%     for u = 1:numel(uKeys)
%         acc = ps{t};
%         t = t + 1;
%         while t <= numel(ks) && ks(t) == uKeys(u)
%             acc = acc + ps{t};
%             t = t + 1;
%         end
%         uParams{u} = acc;
%     end
% 
%     % return in your original row-ish shape (optional)
%     keys   = uKeys.';      % row
%     params = uParams.';    % row cell
% end

end
