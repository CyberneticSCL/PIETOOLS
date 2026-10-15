% Test_tensopvar_dims.m -- tests C.dims.in and C.dims.out for 'tensopvar'.
%
% For random 'tensopvar' objects C with a known d x 2 array of dimensions of
% the factors, C.dims.in must return dims(:,2) and C.dims.out must return
% dims(:,1), also when followed by an index, e.g. C.dims.out(1). The
% remaining fields, including C.dims itself, are not affected.
%
% DJ (CC), 10/15/2026: Initial coding;

for rep=1:12
    rng(rep);
    d = 1+mod(rep,3);       N = 1+mod(rep,2);       M = 1+mod(floor(rep/2),3);      nterms = 1+mod(rep,2);
    type = [1+mod(rep,2),2];
    dims = randi(3,d,2);    if type(1)==1,  dims(:,1) = dims(1,1);  end
    vars = struct('in',{arrayfun(@(i) sprintf('x%d',i),1:N,'UniformOutput',false)},'out',{arrayfun(@(i) sprintf('y%d',i),1:M,'UniformOutput',false)});
    dom = struct('in',[0,1;-1,2;0,3],'out',[0,2;1,4;-1,0]);        dom.in = dom.in(1:N,:);     dom.out = dom.out(1:M,:);
    depmat = struct('in',rand(d,N)>0.4,'out',rand(d,M)>0.4);
    depmat.in(1,:) = true;      depmat.out(1,:) = true;
    Top = rand_tensopvar(dims,vars,dom,depmat,nterms,type,2,0.5);

    if ~isequal(Top.dims.in,dims(:,2)) || ~isequal(Top.dims.out,dims(:,1))
        error("Test_tensopvar_dims failed: dims.in and dims.out must be the second and first column of dims (rep=%d).",rep)
    elseif ~isequal(size(Top.dims.in),[d,1]) || ~isequal(size(Top.dims.out),[d,1])
        error("Test_tensopvar_dims failed: dims.in and dims.out must be column vectors (rep=%d).",rep)
    end
    % Followed by an index
    k = randi(d);
    if ~isequal(Top.dims.in(k),dims(k,2)) || ~isequal(Top.dims.out(end),dims(end,1))
        error("Test_tensopvar_dims failed: indexing of dims.in and dims.out (rep=%d).",rep)
    end
    % Other fields are unaffected
    if ~isequal(Top.dims,dims) || ~isequal(Top.vars.in,vars.in) || ~isequal(Top.dom.out,dom.out)
        error("Test_tensopvar_dims failed: other fields must be unaffected (rep=%d).",rep)
    end
end

% Other fields cannot be accessed through dims
try
    Top.dims.foo;
    error("Test_tensopvar_dims failed: an invalid field of dims must error.")
catch ME
    if contains(ME.message,"Test_tensopvar_dims failed")
        rethrow(ME)
    end
end

disp("Test_tensopvar_dims passed.")
