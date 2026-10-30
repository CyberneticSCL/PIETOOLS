% Test_convert_PDE_nonlinear_varnames.m -- names and order of the PIE
% states returned by convert_PIETOOLS_PDE_nonlinear. The kth variable of
% the 'polyopvar' PIE.x (which sorts its variables by name) must be the
% kth state component in x_tab, i.e. the kth column of the operator T and
% the kth state in PIE.space.x.varname and the PDE state with name
% PIE.space_pde.x.varname(PIE.pde2pie.x.idx(k,1)), also if the
% converter reorders the state components, and if there are 10 or more
% states. PIE.space.x.varname must be sorted, and the same as PIE.x.varname
% and the names in f.
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ, 10/29/2026: Use the names in PIE.space.x.varname and
%                   PIE.space_pde.x.varname, instead of PIE.varname_pie and
%                   PIE.varname_pde (which are removed);
%
% NOTES:
% - 'clear stateNameGenerator' restarts the state numbering in each case,
%   so the PDE states are named x1, x2, ... as the checks (space_pde.x.varname =
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
%   if ~strcmp(PIE.varname_pde{1},'x2')                                     % DJ, 10/29/2026 (was)
if ~strcmp(PIE.space_pde.x.varname{PIE.pde2pie.x.idx(1,1)},'x2')            % DJ, 10/29/2026
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
names = PIE.space.x.varname;                                                % DJ, 10/29/2026
names_pde = PIE.space_pde.x.varname(PIE.pde2pie.x.idx(:,1));                % DJ, 10/29/2026
if ~isequal(names,sort(names))
    error("Test_convert_PDE_nonlinear_varnames failed (%s): PIE.space.x.varname is not sorted.",label)
elseif ~isequal(reshape(PIE.x.varname,[],1),reshape(names,[],1))
    error("Test_convert_PDE_nonlinear_varnames failed (%s): PIE.x.varname differs from PIE.space.x.varname.",label)
elseif numel(names_pde)~=nx || numel(names)~=nx
    error("Test_convert_PDE_nonlinear_varnames failed (%s): number of names differs from the number of states.",label)
elseif ~all(ismember(PIE.f.varname,names))
    error("Test_convert_PDE_nonlinear_varnames failed (%s): the names in f are not those in PIE.space.x.varname.",label)
end
for i=1:nx
    if ~strcmp(names_pde{i},['x',num2str(PIE.x_tab(i,1))])
        error("Test_convert_PDE_nonlinear_varnames failed (%s): the PDE state names do not match x_tab.",label)
    end
end
end