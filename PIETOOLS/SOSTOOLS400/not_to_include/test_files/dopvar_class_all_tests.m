function [] = dopvar_class_all_tests(varargin)
testtypes = {'init','plus','mtimes','transpose','horzcat','vertcat'};
ntests = length(testtypes);
testerror = zeros(ntests,1);
tol=1e-13;
if nargin==0
    disp('Running all tests');
    logval = ones(ntests,1);
else
    logval = zeros(ntests,1);
    for i=1:nargin
        if strcmp(varargin{i},'init')
            logval(1) = 1;
        end
        if strcmp(varargin{i},'plus')
            logval(2) = 1;
        end
        if strcmp(varargin{i},'mtimes')
            logval(3) = 1;
        end
        if strcmp(varargin{i},'transpose')
            logval(4) = 1;
        end
        if strcmp(varargin{i},'horzcat')
            logval(5) = 1;
        end
        if strcmp(varargin{i},'vertcat')
            logval(6) = 1;
        end
    end
end

% test creation, create random dpvar
if logval(1)
    disp('Running creation tests');
    disp('----------------------');
    disp('testing no dvars')
    opvardim = randi([0,4],2,2);
    ndvars = 0;
    npvars = 2;
    monodeg = randi([1,20],1,1);
    sparsity = rand(1);
    testdop = random_dopvar(opvardim,monodeg,sparsity,ndvars,npvars-1,npvars);
    err = testdop.chkval;
    disp('testing no pvars')
    opvardim = randi([0,4],2,2);
    ndvars = randi([0,2000],1,1);
    npvars = 0;
    monodeg = 1;
    sparsity = rand(1);
    testdop = random_dopvar(opvardim,monodeg,sparsity,ndvars,0,0);
    err = err&testdop.chkval;
    disp('testing with low density')
    opvardim = randi([0,4],2,2);
    ndvars = randi([0,2000],1,1);
    npvars = 2;
    monodeg = randi([1,20],1,1);
    sparsity = 0.2*rand(1);
    testdop = random_dopvar(opvardim,monodeg,sparsity,ndvars,npvars-1,npvars);
    err = err&testdop.chkval;
    disp('testing with high density')
    opvardim = randi([0,4],2,2);
    ndvars = randi([0,2000],1,1);
    npvars = 2;
    monodeg = randi([1,20],1,1);
    sparsity = 0.8*rand(1)+0.2;
    testdop = random_dopvar(opvardim,monodeg,sparsity,ndvars,npvars-1,npvars);
    err = err&testdop.chkval;
    disp('testing R to R')
    opvardim = [randi([1,4],1,2); 0 0];
    ndvars = randi([0,2000],1,1);
    npvars = 2;
    monodeg = randi([1,20],1,1);
    sparsity = rand(1);
    testdop = random_dopvar(opvardim,monodeg,sparsity,ndvars,npvars-1,npvars);
    err = err&testdop.chkval;
    disp('testing R to L2')
    opvardim = [0, randi([1,4],1,1); randi([1,4],1,1), 0];
    ndvars = randi([0,2000],1,1);
    npvars = 2;
    monodeg = randi([1,20],1,1);
    sparsity = rand(1);
    testdop = random_dopvar(opvardim,monodeg,sparsity,ndvars,npvars-1,npvars);
    err = err&testdop.chkval;
    disp('testing L2 to R')
    opvardim = [randi([1,4],1,1),0;0, randi([1,4],1,1)];
    ndvars = randi([0,2000],1,1);
    npvars = 2;
    monodeg = randi([1,20],1,1);
    sparsity = rand(1);
    testdop = random_dopvar(opvardim,monodeg,sparsity,ndvars,npvars-1,npvars);
    err = err&testdop.chkval;
    disp('testing L2 to L2')
    opvardim = [0,0;randi([1,4],1,2)];
    ndvars = randi([0,2000],1,1);
    npvars = 2;
    monodeg = randi([1,20],1,1);
    sparsity = rand(1);
    testdop = random_dopvar(opvardim,monodeg,sparsity,ndvars,npvars-1,npvars);
    err = err&testdop.chkval;
    disp('testing RL2 to RL2')
    opvardim = randi([1,4],2,2);
    ndvars = randi([0,2000],1,1);
    npvars = 2;
    monodeg = randi([1,20],1,1);
    sparsity = rand(1);
    testdop = random_dopvar(opvardim,monodeg,sparsity,ndvars,npvars-1,npvars);
    err = err&testdop.chkval;
    if ~err
        testerror(1) = 1;
    end
end

if logval(2) % running plus and minus tests
    disp('Running addition tests');
    disp('----------------------');
    disp('no dvars');
    opvardim = randi([0,4],2,2);
    ndvars = 0;
    npvars = 2;
    monodeg = randi([1,20],1,1);
    sparsity = rand(1);
    
    dopvar1 = random_dopvar(opvardim,monodeg,sparsity,ndvars,npvars-1,npvars);
    dopvar2 = random_dopvar(opvardim,monodeg,sparsity,ndvars,npvars-1,npvars);
    dopvarplus = dopvar1+dopvar2;
    opvarplus = dopvar2opvar(dopvar1)+dopvar2opvar(dopvar2);
    err_plus = error_dopvar(dopvarplus-opvarplus);
    dopvarminus = dopvar1-dopvar2;
    opvarminus = dopvar2opvar(dopvar1)-dopvar2opvar(dopvar2);
    err_minus = error_dopvar(dopvarminus-opvarminus);
    disp('no pvars');
    opvardim = randi([0,4],2,2);
    ndvars = randi([0,2000],1,1);
    npvars = 0;
    monodeg = 1;
    sparsity = rand(1);
    
    dopvar1 = random_dopvar(opvardim,monodeg,sparsity,ndvars,0,0);
    dopvar2 = random_dopvar(opvardim,monodeg,sparsity,ndvars,0,0);
    dopvarplus = dopvar1+dopvar2;
    opvarplus = dopvar2opvar(dopvar1)+dopvar2opvar(dopvar2);
    err_plus = err_plus+error_dopvar(dopvarplus-opvarplus);
    dopvarminus = dopvar1-dopvar2;
    opvarminus = dopvar2opvar(dopvar1)-dopvar2opvar(dopvar2);
    err_minus = err_minus+error_dopvar(dopvarminus-opvarminus);
    
    disp('adding RL2 to RL2')
    opvardim = randi([0,4],2,2);
    ndvars = randi([0,2000],1,1);
    npvars = 2;
    monodeg = randi([1,20],1,1);
    sparsity = rand(1);
    
    dopvar1 = random_dopvar(opvardim,monodeg,sparsity,ndvars,npvars-1,npvars);
    dopvar2 = random_dopvar(opvardim,monodeg,sparsity,ndvars,npvars-1,npvars);
    dopvarplus = dopvar1+dopvar2;
    opvarplus = dopvar2opvar(dopvar1)+dopvar2opvar(dopvar2);
    err_plus = err_plus+error_dopvar(dopvarplus-opvarplus);
    dopvarminus = dopvar1-dopvar2;
    opvarminus = dopvar2opvar(dopvar1)-dopvar2opvar(dopvar2);
    err_minus = err_minus+error_dopvar(dopvarminus-opvarminus);
    if (err_plus>tol)||(err_minus>tol)
        testerror(2) = 1;
    end
end

if logval(3) % testing mtimes
    disp('Running mtimes tests');
    disp('----------------------');
    disp('no dvars');
    opvardim = randi([0,4],2,2);
    dopvardim = [opvardim(:,2), randi([0,4],2,1)];
    ndvars = 0;
    npvars = 2;
    monodeg = randi([1,20],1,1);
    sparsity = rand(1);
    
    dopvar1 = random_dopvar(opvardim,monodeg,sparsity,0,npvars-1,npvars);
    opvar1 = dopvar2opvar(dopvar1);
    dopvar2 = random_dopvar(dopvardim,monodeg,sparsity,ndvars,npvars-1,npvars);
    dopvarM1 = opvar1*dopvar2;
    dopvarM2 = dopvar2'*opvar1';
    err_mtimes = error_dopvar(dopvarM1-dopvarM2');
    disp('no pvars');
    opvardim = randi([0,4],2,2);
    dopvardim = [opvardim(:,2), randi([0,4],2,1)];
    ndvars = randi([0,2000],1,1);
    npvars = 0;
    monodeg = 1;
    sparsity = rand(1);
    
    dopvar1 = random_dopvar(opvardim,monodeg,sparsity,0,0,0);
    opvar1 = dopvar2opvar(dopvar1);
    dopvar2 = random_dopvar(dopvardim,monodeg,sparsity,ndvars,0,0);
    dopvarM1 = opvar1*dopvar2;
    dopvarM2 = dopvar2'*opvar1';
    err_mtimes = err_mtimes+error_dopvar(dopvarM1-dopvarM2');
    
    disp('mtimes R to R')
    opvardim = [randi([1,4],1,2); 0 0];
    dopvardim = [opvardim(:,2), [randi([1,4],1,1);0]];
    ndvars = randi([0,2000],1,1);
    npvars = 2;
    monodeg = randi([1,20],1,1);
    sparsity = rand(1);
    
    dopvar1 = random_dopvar(opvardim,monodeg,sparsity,0,npvars-1,npvars);
    opvar1 = dopvar2opvar(dopvar1);
    dopvar2 = random_dopvar(dopvardim,monodeg,sparsity,ndvars,npvars-1,npvars);
    dopvarM1 = opvar1*dopvar2;
    dopvarM2 = dopvar2'*opvar1';
    err_mtimes = err_mtimes+error_dopvar(dopvarM1-dopvarM2');
    
    disp('mtimes L2 to L2')
    opvardim = [0 0;randi([1,4],1,2)];
    dopvardim = [opvardim(:,2), [0;randi([1,4],1,1)]];
    ndvars = randi([0,2000],1,1);
    npvars = 2;
    monodeg = randi([1,20],1,1);
    sparsity = rand(1);
    
    dopvar1 = random_dopvar(opvardim,monodeg,sparsity,0,npvars-1,npvars);
    opvar1 = dopvar2opvar(dopvar1);
    dopvar2 = random_dopvar(dopvardim,monodeg,sparsity,ndvars,npvars-1,npvars);
    dopvarM1 = opvar1*dopvar2;
    dopvarM2 = dopvar2'*opvar1';
    err_mtimes = err_mtimes+error_dopvar(dopvarM1-dopvarM2');
    
    disp('mtimes RL2 to RL2')
    opvardim = randi([1,4],2,2);
    dopvardim = [opvardim(:,2), randi([0,4],2,1)];
    ndvars = randi([0,2000],1,1);
    npvars = 2;
    monodeg = randi([1,20],1,1);
    sparsity = rand(1);
    
    dopvar1 = random_dopvar(opvardim,monodeg,sparsity,0,npvars-1,npvars);
    opvar1 = dopvar2opvar(dopvar1);
    dopvar2 = random_dopvar(dopvardim,monodeg,sparsity,ndvars,npvars-1,npvars);
    dopvarM1 = opvar1*dopvar2;
    dopvarM2 = dopvar2'*opvar1';
    err_mtimes = err_mtimes+error_dopvar(dopvarM1-dopvarM2');
    if (err_mtimes>tol)
        testerror(3) = 1;
    end
end    
if logval(4) % testing transpose
    disp('Running transpose test')
    disp('----------------------')
    disp('no dvars')
    opvardim = randi([0,4],2,2);
    ndvars = 0;
    npvars = 2;
    monodeg = randi([1,20],1,1);
    sparsity = rand(1);
    
    dopvar1 = random_dopvar(opvardim,monodeg,sparsity,ndvars,npvars-1,npvars);
    opvar1 = dopvar2opvar(dopvar1);
    dopvarT = dopvar1';
    opvarT = opvar1';
    err_T = error_dopvar(dopvarT-opvarT);
    disp('no pvars');
    opvardim = randi([0,4],2,2);
    ndvars = randi([0,2000],1,1);
    npvars = 0;
    monodeg = 1;
    sparsity = rand(1);
    
    dopvar1 = random_dopvar(opvardim,monodeg,sparsity,ndvars,0,0);
    opvar1 = dopvar2opvar(dopvar1);
    dopvarT = dopvar1';
    opvarT = opvar1';
    err_T = err_T+error_dopvar(dopvarT-opvarT);
    
    disp('RL2')
    opvardim = randi([0,4],2,2);
    ndvars = randi([0,2000],1,1);
    npvars = 2;
    monodeg = randi([1,20],1,1);
    sparsity = rand(1);
    
    dopvar1 = random_dopvar(opvardim,monodeg,sparsity,ndvars,npvars-1,npvars);
    opvar1 = dopvar2opvar(dopvar1);
    dopvarT = dopvar1';
    opvarT = opvar1';
    err_T = err_T+error_dopvar(dopvarT-opvarT);
    if (err_T>tol)
        testerror(4) = 1;
    end
end
if logval(5) % testing horzcat
    
end
if logval(6) % testing vertcat

end
if ~any(testerror)
    disp('********************');
    disp('No errors were found');
else
    errormethods = testtypes(testerror~=0);
    errormethods = cellfun(@(x) strcat(x,', '),errormethods,'un',0);
    errmsg = strcat('Errors were found in:',errormethods{:});
    disp(errmsg);
end

end

function testDop = random_dopvar(opvardim, monodeg, sparsity, ndvars, npvars1, npvars2) 
dopvar testDop;

testDop.dim = opvardim; testDop.I = [0,1]; testDop.var1 = pvar('pvar_1');
testDop.var2 = pvar('pvar_2');
testDop.P = random_dpvar(opvardim(1,:), 1, sparsity, ndvars, 0);
testDop.Q1 = random_dpvar([opvardim(1,1),opvardim(2,2)], monodeg, sparsity, ndvars, npvars1);
testDop.Q2 = random_dpvar([opvardim(2,1),opvardim(1,2)], monodeg, sparsity, ndvars, npvars1);
testDop.R.R0 = random_dpvar(opvardim(2,:), monodeg, sparsity, ndvars, npvars1);
testDop.R.R1 = random_dpvar(opvardim(2,:), 2*monodeg, sparsity, ndvars, npvars2);
testDop.R.R2 = random_dpvar(opvardim(2,:), 2*monodeg, sparsity, ndvars, npvars2);
end

function testDp = random_dpvar(matdim, monodim, sparsity, ndvars, npvars)
m = matdim(1); n = matdim(2);
dvarname = {}; varname = {};
for i=1:npvars
    varname = [varname, ['pvar_',num2str(i)]];
end
for i=1:ndvars
    dvarname = [dvarname, ['coeff_',num2str(i)]];
end

[~,idxp] = sort(rand(1,npvars));
[~,idxd] = sort(rand(1,ndvars));

if npvars
    degmat = floor(5*sprand(monodim,npvars,sparsity));
else
    degmat = zeros(1,npvars);
end
degmat = unique(degmat,'rows');

C = sprand((ndvars+1)*m,size(degmat,1)*n,sparsity);
testDp = dpvar(C,degmat,varname(idxp),dvarname(idxd),matdim);
end
function err= error_dopvar(dopVar)
P = dopVar.P;
Q1 = dopVar.Q1; Q2 = dopVar.Q2; R0 = dopVar.R.R0; R1 = dopVar.R.R1; R2 = dopVar.R.R2;
err = norm(full(P.C))+ norm(full(Q1.C))+norm(full(Q2.C))+norm(full(R0.C))+norm(full(R1.C))+norm(full(R2.C));
end