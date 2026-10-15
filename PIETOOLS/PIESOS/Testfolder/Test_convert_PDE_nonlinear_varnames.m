% Test_convert_PDE_nonlinear_varnames.m -- names and order of the PIE
% states returned by convert_PIETOOLS_PDE_nonlinear. The kth variable of
% the 'polyopvar' PIE.x (which sorts its variables by name) must be the
% kth state component in x_tab, i.e. the kth column of the operator T and
% the kth state in PIE.varname_pie and PIE.varname_pde, also if the
% converter reorders the state components, and if there are 10 or more
% states. PIE.varname_pie must be sorted, and the same as PIE.x.varname
% and the names in f.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES:
% - 'clear stateNameGenerator' restarts the state numbering in each case,
%   so the PDE states are named x1, x2, ... as the checks (varname_pde =
%   'x<x_tab(i,1)>') assume.

pvar s t
dom = [0,1];
ncases = 0;

% % Reordered components: the second order state u1 is declared before the first order state u2,
% % and the converter puts u2 first
clear stateNameGenerator
u1 = pde_var(s,dom);        u2 = pde_var(s,dom);
PDE = [diff(u1,t)==diff(u1,s,2)-u1+u2;
       diff(u2,t)==-diff(u2,s)-u2;
       subs(u1,s,dom(1))==0;    subs(u1,s,dom(2))==0;
       subs(u2,s,dom(1))==0];
evalc('PIE = convert_PIETOOLS_PDE_nonlinear(PDE,''silent'');');
if isequal(PIE.x_tab(:,1),sort(PIE.x_tab(:,1)))
    error("Test_convert_PDE_nonlinear_varnames: the converter did not reorder the components, the test is vacuous.")
end
check_names(PIE,'reordered');
% The first state is the first order state u2: T = int_0^s
if ~strcmp(PIE.varname_pde{1},'x2')
    error("Test_convert_PDE_nonlinear_varnames failed: the first state is not the PDE state x2.")
end
v = 3*s^2+1;                                        % state value of the first order state
Y = polynomial(apply_op(PIE.T,{v,0*s}));            % T*[v;0]
% expected: u2 = int_0^s v = s^3+s (first row), u1 = 0 (second row)
D1 = Y(1)-(s^3+s);      D2 = Y(2);
if any(abs(full(D1.coefficient))>1e-10) || any(abs(full(D2.coefficient))>1e-10)
    error("Test_convert_PDE_nonlinear_varnames failed: T does not act on the first state as on the first order state.")
end
ncases = ncases+1;

% % More than 9 states: the names are zero-padded, so that sorting preserves the order
clear stateNameGenerator
n = 11;
for i=1:n
    ui = pde_var(s,dom);
    PDEi = [diff(ui,t)==diff(ui,s,2)-i*ui;   subs(ui,s,dom(1))==0;   subs(ui,s,dom(2))==0];
    if i==1,    PDE = PDEi;     else,   PDE = [PDE;PDEi];   end
end
evalc('PIE = convert_PIETOOLS_PDE_nonlinear(PDE,''silent'');');
check_names(PIE,'11 states');
ncases = ncases+1;

fprintf('Test_convert_PDE_nonlinear_varnames passed (%d cases).\n',ncases)


function check_names(PIE,label)
nx = size(PIE.x_tab,1);
names = PIE.varname_pie;
if ~isequal(names,sort(names))
    error("Test_convert_PDE_nonlinear_varnames failed (%s): PIE.varname_pie is not sorted.",label)
elseif ~isequal(reshape(PIE.x.varname,[],1),reshape(names,[],1))
    error("Test_convert_PDE_nonlinear_varnames failed (%s): PIE.x.varname differs from PIE.varname_pie.",label)
elseif numel(PIE.varname_pde)~=nx || numel(names)~=nx
    error("Test_convert_PDE_nonlinear_varnames failed (%s): number of names differs from the number of states.",label)
elseif ~all(ismember(PIE.f.varname,names))
    error("Test_convert_PDE_nonlinear_varnames failed (%s): the names in f are not those in PIE.varname_pie.",label)
end
for i=1:nx
    if ~strcmp(PIE.varname_pde{i},['x',num2str(PIE.x_tab(i,1))])
        error("Test_convert_PDE_nonlinear_varnames failed (%s): varname_pde does not match x_tab.",label)
    end
end
end
