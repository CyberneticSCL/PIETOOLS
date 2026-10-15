function [Aop,Bop,info] = bench_build(spec)
% [AOP,BOP,INFO] = BENCH_BUILD(SPEC) builds a random 'intop' object
% Aop (with 'dpvar' kernels in q decision variables) and a random
% 'tensopvar' object Bop, to benchmark mtimes_tensopvar(Aop,Bop).
%
% INPUTS
% - spec:   struct specifying the case, with fields
%               'N':      (required) number of shared spatial variables,
%                         s1,...,sN on [0,1];
%               'd':      (required) number of factors in Bop, each
%                         depending on all N variables;
%               'q':      (required) number of decision variables in the
%                         kernels of Aop (0: polynomial kernels);
%               'nterms': number of terms in Bop. Defaults to 1;
%               'kq':     number of nonzeros per decision variable and
%                         matrix row in each kernel. Defaults to 2;
%               'dnsty':  density of the sparse coefficient matrices of
%                         the factors. Defaults to 0.5;
%               'pdegB':  degree of the kernels of Bop in each variable.
%                         Defaults to 1;
%               'pdegA':  degree of the kernels of Aop in each variable.
%                         Defaults to 1;
%               'seed':   seed for rng. Defaults to 1;
%               'dvmode': 'same' (default), all kernels carry the same
%                         list of decision variables, or 'split', kernel
%                         k carries a different random half of the list;
%               'afam':   (optional) 'full', 'quad' or 'cubic', replaces
%                         the integral structure of Aop ('quad' and
%                         'cubic' require d=2 and d=3, respectively);
%               'kA':     (optional, N=1) keep only the first kA integral
%                         terms of Aop;
%
% OUTPUTS
% - Aop:    'intop' object;
% - Bop:    'tensopvar' object;
% - info:   struct with field 'nkernels', the number of kernels of Aop;
%
% NOTES
% - intop.params are sparse coefficient matrices, not kernels, so
%   kernels(Aop) is modified and Aop declared from the result at the end.
% - Kernels of Aop only involve monomials with a nonzero coefficient, so
%   for q>0 the random draws depend on the monomials of the rand_intop
%   kernels that have a nonzero coefficient.

% DJ (CC), 10/15/2026: Initial coding;

if ~isfield(spec,'nterms'),  spec.nterms = 1;       end
if ~isfield(spec,'kq'),      spec.kq = 2;           end
if ~isfield(spec,'dnsty'),   spec.dnsty = 0.5;      end
if ~isfield(spec,'pdegB'),   spec.pdegB = 1;        end
if ~isfield(spec,'pdegA'),   spec.pdegA = 1;        end
if ~isfield(spec,'seed'),    spec.seed = 1;         end
if ~isfield(spec,'dvmode'),  spec.dvmode = 'same';  end
rng(spec.seed);
N = spec.N;     d = spec.d;
vars = struct('in',{arrayfun(@(i) sprintf('s%d',i),1:N,'UniformOutput',false)},'out',[]);
vars.out = vars.in;
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
depmat = struct('in',true(d,N),'out',true(d,N));
pdegsB = struct('in',spec.pdegB*ones(1,N),'out',spec.pdegB*ones(1,N));
Bop = rand_tensopvar(ones(d,2),vars,dom,depmat,spec.nterms,[2,2],pdegsB,spec.dnsty);
m = Bop.matdim(1);
Aop = rand_intop([1,m],Bop.vars.out,Bop.dom.out,Bop.deps.out,spec.pdegA*ones(1,N),spec.q);
Aparams = kernels(Aop);
Aintmat = Aop.intmat;
if isfield(spec,'afam') && strcmp(spec.afam,'full')
    % Single term: full domain in each variable, first kernel
    intmat_full = cell(1,N);
    for i=1:N
        intmat_full{i} = repmat([0,d+1],1,d);
    end
    Aparams = Aparams(1);       Aintmat = intmat_full;
end
if isfield(spec,'afam') && any(strcmp(spec.afam,{'quad','cubic'}))
    % Same distinct regions in each variable, independent random kernel per combination
    if strcmp(spec.afam,'quad'),    rows = [0 2 0 3; 0 3 0 3; 0 3 1 3];
    else,                           rows = [0 2 0 3 0 4; 3 4 0 3 0 4; 0 4 0 4 0 4];
    end
    assert(size(rows,2)==2*d,'afam quad/cubic needs d=2/3');
    r = size(rows,1);
    names = cell(N*d,1);
    for i=1:N,  for k=1:d,  names{(i-1)*d+k} = sprintf('s%d_dum%d',i,k);   end,  end
    pv = polynomial(names);
    params = cell([r*ones(1,N),1]);
    for k=1:numel(params)
        params{k} = rand_poly([1,m],pv,spec.pdegA,50);
    end
    Aparams = params;       Aintmat = repmat({rows},1,N);
end
if isfield(spec,'kA') && N==1
    kA = min(spec.kA,size(Aintmat{1},1));
    Aparams = Aparams(1:kA);    Aintmat = {Aintmat{1}(1:kA,:)};
end
info = struct('nkernels',numel(Aparams));
if spec.q>0
    dvn = compose('coeff%d',(1:spec.q)');
    for k=1:numel(Aparams)
        P = Aparams{k};
        [mm,nn] = size(P);      nZ = size(P.degmat,1);
        if strcmp(spec.dvmode,'split')
            J = sort(randperm(spec.q,ceil(spec.q/2)))';
        else
            J = (1:spec.q)';
        end
        nJ = numel(J);
        % rows: matrix row outer, (constant term, decision variables) inner
        r = repmat((0:mm-1),nJ*spec.kq,1);
        rows = (r(:)'*(nJ+1)) + 1 + repmat(repelem(1:nJ,spec.kq),1,mm);
        cols = randi(nn*nZ,1,numel(rows));
        vals = randn(1,numel(rows));
        % constant term
        rows = [rows,(0:mm-1)*(nJ+1)+1];
        cols = [cols,randi(nn*nZ,1,mm)];
        vals = [vals,randn(1,mm)];
        C = sparse(rows,cols,vals,mm*(nJ+1),nn*nZ);
        Aparams{k} = dpvar(C,P.degmat,P.varname,dvn(J),[mm,nn]);
    end
end
% Declare the functionals from the kernels
Aop = intop(Aparams,Aintmat,Aop.vars,Aop.dom);
end
