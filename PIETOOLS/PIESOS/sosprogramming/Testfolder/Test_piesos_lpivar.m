function Test_piesos_lpivar()
% TEST_PIESOS_LPIVAR tests, using assertions, piesos_lpivar,
%   [prog,Pop,Pmat,Zop] = ...
%               piesos_lpivar(prog,dims,vars,dom,depmat,degs,opts),
% which declares the indefinite operator variable Pop = Pmat*Zop, returned
% as a struct with fields P=Pmat and Z=Zop.
%
% Ground truths are independent of the function under test:
%  - monomial counts, by hand: for a state and output on one variable s,
%    the 3-PI operator has a multiplier term (kernel in s) and a lower and
%    an upper integral term (kernels in s and the dummy variable). With
%    maximal degree d per variable these have d+1, (d+1)^2 and (d+1)^2
%    monomials, so Zop has 3 rows for d=0 and 2+4+4=10 for d=1; excluding
%    the integral terms leaves the d+1 monomials of the multiplier;
%  - decision variables: Pmat has one decision variable for each of its
%    entries, and as many columns as Zop has rows, one block row per
%    output;
%  - the blocks of Zop of an output are those of monopmat, with the output
%    dependence of that output (stacked over the outputs);
%  - errors.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_piesos_lpivar
%
% Copyright (C) 2026 PIETOOLS Team
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
% DJ (CC), 10/15/2026: Initial coding;

tests = {@case_scalar_counts, @case_defaults, @case_vector_output, @case_blocks_match_monopmat, ...
         @case_exclude, @case_functional, @case_errors};
nfail = 0;
for i=1:numel(tests)
    try
        tests{i}();
        fprintf('PASS  %s\n',func2str(tests{i}));
    catch err
        nfail = nfail+1;
        fprintf('FAIL  %s\n      %s\n',func2str(tests{i}),err.message);
        for k=1:min(3,numel(err.stack))
            fprintf('        in %s, line %d\n',err.stack(k).name,err.stack(k).line);
        end
    end
end
if nfail>0
    error('Test_piesos_lpivar: %d of %d cases failed.',nfail,numel(tests));
end
fprintf('Test_piesos_lpivar: all %d cases passed.\n',numel(tests));
end


%% ---- helpers ----------------------------------------------------------
function [prog,dims,vars,dom,depmat] = setup(dimout,depout)
x1 = polyopvar('x1',polynomial({'s'}),[0,1]);
prog = piesos_program(x1);
dims = struct('in',1,'out',dimout);
vars = struct('in',{{'s'}},'out',{{'s'}});
dom = struct('in',[0 1],'out',[0 1]);
depmat = struct('in',1,'out',depout);
end

function assert_err(f,needle)
try
    f();
catch err
    assert(contains(err.message,needle),'Wrong error message: "%s".',err.message);
    return
end
error('Expected an error containing "%s", none thrown.',needle);
end

function nr = nrows(Zop)
nr = size(Zop,1);
end


%% ---- cases ------------------------------------------------------------
function case_scalar_counts()
[prog,dims,vars,dom,depmat] = setup(1,true);
for d=0:2
    n0 = numel(prog.decvartable);
    [prog1,Pop,Pmat,Zop] = piesos_lpivar(prog,dims,vars,dom,depmat,d);
    nZ = (d+1) + 2*(d+1)^2;
    assert(nrows(Zop)==nZ,'d=%d: expected %d rows of Zop, got %d.',d,nZ,nrows(Zop));
    assert(isequal(size(Pmat),[1,1]) && isequal(size(Pmat{1,1}),[1,nZ]),'d=%d: Pmat must be 1 x %d.',d,nZ);
    assert(numel(prog1.decvartable)-n0==nZ,'d=%d: %d decision variables expected.',d,nZ);
    assert(numel(unique(prog1.decvartable))==numel(prog1.decvartable),'Decision variables must be distinct.');
    assert(isstruct(Pop) && isequal(sort(fieldnames(Pop)),{'P';'Z'}),'Pop must be a struct with fields P and Z.');
    assert(isequal(size(Pop.P),size(Pmat)) && isequal(size(Pop.Z),size(Zop)),'Pop.P and Pop.Z must be Pmat and Zop.');
    assert(isequal(Pop.P{1,1}.dvarname,Pmat{1,1}.dvarname),'Pop.P must hold the decision variables of Pmat.');
end
end

function case_defaults()
% degs omitted = degree 0; opts omitted = default options
[prog,dims,vars,dom,depmat] = setup(1,true);
[~,~,~,Z5] = piesos_lpivar(prog,dims,vars,dom,depmat);
[~,~,~,Z6] = piesos_lpivar(prog,dims,vars,dom,depmat,0);
[~,~,~,Z7] = piesos_lpivar(prog,dims,vars,dom,depmat,0,struct());
assert(nrows(Z5)==3 && nrows(Z6)==3 && nrows(Z7)==3,'Degree 0: 3 rows of Zop expected.');
depmat.out = 1;
[~,~,~,Z8] = piesos_lpivar(prog,dims,vars,dom,depmat,1);
assert(nrows(Z8)==10,'Numeric output dependence 1 must act as true.');
end

function case_vector_output()
% outputs of dimensions 2 and 1, both on s: one block row of Pmat per output, 3 monomials
% (d=0) of Zop per output, 6 rows in all
[prog,dims,vars,dom,depmat] = setup([2;1],[true;true]);
n0 = numel(prog.decvartable);
[prog,Pop,Pmat,Zop] = piesos_lpivar(prog,dims,vars,dom,depmat,0);
assert(nrows(Zop)==6,'Expected 6 rows of Zop, got %d.',nrows(Zop));
assert(isequal(size(Pmat),[2,2]),'Pmat must be a 2 x 2 cell: an output block row, a block of Zop for each output.');
nr = [2,1];
for i=1:2
    for j=1:2
        assert(isequal(size(Pmat{i,j}),[nr(i),3]),'Pmat{%d,%d} must be %d x 3.',i,j,nr(i));
    end
end
assert(numel(prog.decvartable)-n0==2*3+2*3+1*3+1*3,'18 decision variables expected.');
assert(isequal(size(Pop.P),[2,2]) && isequal(size(Pop.P{2,1}),[1,3]) && isequal(size(Pop.Z),size(Zop)),'Pop.P and Pop.Z must be Pmat and Zop.');
names = [Pmat{1,1}.dvarname;Pmat{1,2}.dvarname;Pmat{2,1}.dvarname;Pmat{2,2}.dvarname];
assert(numel(unique(names))==18,'Each entry of Pmat must be a distinct decision variable.');
% output 1 a function of s (3 monomials for d=0), output 2 constant in s (a
% functional of the state: 1 monomial)
[prog,dims,vars,dom,~] = setup([1;1],[true;true]);
depmat = struct('in',1,'out',[true;false]);
[~,~,~,Zop2] = piesos_lpivar(prog,dims,vars,dom,depmat,0);
assert(nrows(Zop2)==4,'Expected 3+1 rows of Zop, got %d.',nrows(Zop2));
end

function case_blocks_match_monopmat()
% Zop stacks monopmat(dims.in,...) for each output, with the depmat of that output
[prog,dims,vars,dom,depmat] = setup([1;1],[true;true]);
[~,~,~,Zop] = piesos_lpivar(prog,dims,vars,dom,depmat,1);
dm = depmat;    dm.out = repmat(depmat.out(1,:),size(depmat.in,1),1);
Zm = monopmat(dims.in,vars,dom,dm,1,struct());
assert(nrows(Zop)==2*nrows(Zm),'Zop must stack two copies of monopmat.');
assert(isequal(Zop.dims.out,[Zm.dims.out;Zm.dims.out]),'Block dimensions must be those of monopmat, stacked.');
assert(isequal(Zop.dims.in,Zm.dims.in),'Input dimensions must be those of monopmat.');
assert(isequal(Zop.vars,Zm.vars) && isequal(Zop.dom,Zm.dom),'Variables and domains must be those of monopmat.');
end

function case_exclude()
% excluding the integral terms leaves the d+1 monomials of the multiplier
[prog,dims,vars,dom,depmat] = setup(1,true);
opts = struct('exclude',{{[false;true;true]}});
[~,~,Pmat,Zop] = piesos_lpivar(prog,dims,vars,dom,depmat,1,opts);
assert(nrows(Zop)==2 && isequal(size(Pmat{1,1}),[1,2]),'Only the multiplier monomials (2) must remain.');
% excluding the multiplier leaves the integrals: 4+4
opts = struct('exclude',{{[true;false;false]}});
[~,~,~,Zop] = piesos_lpivar(prog,dims,vars,dom,depmat,1,opts);
assert(nrows(Zop)==8,'Only the integral monomials (8) must remain.');
end

function case_functional()
% an output without spatial variable: a functional of the state, whose
% kernel (the integral term) has the monomials of degree <=d in s: 2 for d=1
[prog,dims,~,~,~] = setup(1,true);
vars = struct('in',{{'s'}},'out',{cell(1,0)});
dom = struct('in',[0 1],'out',zeros(0,2));
depmat = struct('in',1,'out',false(1,0));
[~,~,Pmat,Zop] = piesos_lpivar(prog,dims,vars,dom,depmat,1);
assert(nrows(Zop)==2,'Expected 2 rows of Zop, got %d.',nrows(Zop));
assert(isequal(size(Pmat{1,1}),[1,2]),'Pmat must be 1 x 2.');
end

function case_errors()
[prog,dims,vars,dom,depmat] = setup(1,true);
depmat.out = 2;
assert_err(@() piesos_lpivar(prog,dims,vars,dom,depmat,1),'cannot be set explicitly');
depmat.out = [0.5];
assert_err(@() piesos_lpivar(prog,dims,vars,dom,depmat,1),'cannot be set explicitly');
end
