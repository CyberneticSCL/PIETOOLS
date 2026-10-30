% Test_Conversion_Linear.m -- checks that the nonlinear converter
% `convert_PIETOOLS_PDE_nonlinear` returns the same PIE as the linear
% converter `convert`, for examples 1 through 32 of the PDE library
% `examples_PDE_library_PIETOOLS_nouinput`.
%
% The operators Tw, Tu, B1, B2 and D11-D22 ('tensopmat' objects, 'opvar'
% objects in the linear PIE), and T, f, g and h ('tensopmat' and
% 'polyopvar' objects) are compared by applying them to random polynomial
% states and inputs, up to the relative tolerance tol. The conversion of a
% linear PIE by `PIE2polyPIE` is tested in the same way. An error is thrown
% if any discrepancy is found.
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/29/2026: Use PIE.space.x.varname for the names of the states;
% DJ (CC), 10/29/2026: Use the polynomials gz and gy, instead of h and g;
% DJ (CC), 10/29/2026: The nonlinear converter now returns an
%                   'npie_struct', with the operators Tw, ..., D22 as
%                   'tensopmat' objects, so these are compared by applying
%                   them to random inputs;
%
% NOTES
% - The state components of the PIE of the nonlinear converter may be
%   ordered differently than the PDE states (e.g. in example 32, x3
%   comes before x2), as in the linear converter, see x_tab. The kth
%   variable of PIE.x, 'xk_pie', is the kth component in x_tab (the kth
%   column of T), and PIE.space_pde.x.varname{PIE.pde2pie.x.idx(k,1)} the associated PDE state.

clear
example_nums = 1:32;
tol = 1e-6;                 % the nonlinear converter cleans coefficients below 1e-8
rng(5);
pvar s
nfail = 0;      nchecked = 0;       nskipped = 0;
for example_num = example_nums
    % Convert to PIE using both functions
    clear stateNameGenerator
    PDE = examples_PDE_library_PIETOOLS_nouinput(example_num);
    evalc('PIE_lin = convert(PDE);');
    clear stateNameGenerator
    PDE = examples_PDE_library_PIETOOLS_nouinput(example_num);
    evalc('PIE_poly = convert_PIETOOLS_PDE_nonlinear(PDE,''silent'');');
    msgs = {};
%   if ~isequal(PIE_poly.varname_pie,PIE_poly.x.varname(:)) || ~isequal(PIE_poly.varname_pie,sort(PIE_poly.varname_pie))    % DJ (CC), 10/29/2026 (was)
%       msgs{end+1} = "PIE.varname_pie is not sorted, or differs from PIE.x.varname";     % DJ (CC), 10/29/2026 (was)
    if ~isequal(PIE_poly.space.x.varname,PIE_poly.x.varname(:)) || ~isequal(PIE_poly.space.x.varname,sort(PIE_poly.space.x.varname))    % DJ (CC), 10/29/2026
        msgs{end+1} = "PIE.space.x.varname is not sorted, or differs from PIE.x.varname";             %#ok<SAGROW>
    end

    % % % Start of changes DJ (CC), 10/29/2026: the operators from the inputs are
    % % % 'tensopmat' objects, which are applied to random inputs, as T below.
    %   ops_old = {'Tw','Tu','B1','B2','D11','D12','D21','D22'};            % DJ (CC), 10/29/2026 (was)
    %   for k=1:numel(ops_old)                                              % DJ (CC), 10/29/2026 (was)
    %       name = ops_old{k};                                              % DJ (CC), 10/29/2026 (was)
    %       Pold = PIE_lin.(name);      Ppoly = PIE_poly.(name);            % DJ (CC), 10/29/2026 (was)
    %       if ~isequal(size(Pold),size(Ppoly)) || ~isequal(Pold.dim,Ppoly.dim)     % DJ (CC), 10/29/2026 (was)
    %           msgs{end+1} = sprintf("%s: dimensions differ",name);        % DJ (CC), 10/29/2026 (was)
    %       elseif ~all(all(Pold==Ppoly))                                   % DJ (CC), 10/29/2026 (was)
    %           msgs{end+1} = sprintf("%s: operators differ",name);         % DJ (CC), 10/29/2026 (was)
    %       end                                                             % DJ (CC), 10/29/2026 (was)
    %       nchecked = nchecked+1;                                          % DJ (CC), 10/29/2026 (was)
    %   end                                                                 % DJ (CC), 10/29/2026 (was)
    ops_in = {'Tw','x_tab','w_tab';  'Tu','x_tab','u_tab';
              'B1','x_tab','w_tab';  'B2','x_tab','u_tab';
              'D11','z_tab','w_tab'; 'D12','z_tab','u_tab';
              'D21','y_tab','w_tab'; 'D22','y_tab','u_tab'};
    for k=1:size(ops_in,1)
        name = ops_in{k,1};
        Pold = PIE_lin.(name);      Ppoly = PIE_poly.(name);
        tab_out = PIE_poly.(ops_in{k,2});       tab_in = PIE_poly.(ops_in{k,3});
        if isempty(tab_out) || isempty(tab_in)
            % Operator maps from or to an empty space
            if size(Pold,1)*size(Pold,2)>0 && ~(isempty(Ppoly))
                msgs{end+1} = sprintf("%s: operator should be empty",name); %#ok<SAGROW>
            end
        else
            [Xin,v_in] = random_states(tab_in,s);
            Y_new = order_as_opvar(apply_op(Ppoly,Xin),tab_out);
            Y_old = apply_opvar(Pold,v_in);
            msgs = [msgs,compare_action(Y_old,Y_new,tol,name)];             %#ok<AGROW>
        end
        nchecked = nchecked+1;
    end
    % % % End of changes DJ (CC), 10/29/2026

    % Check the maps from the fundamental state, applied to random states
    x_tab = PIE_poly.x_tab;
    [X,v_old] = random_states(x_tab,s);
    % T: fundamental to PDE state
    Y_new = apply_op(PIE_poly.T,X);
    Y_new = order_as_opvar(Y_new,x_tab);
    Y_old = apply_opvar(PIE_lin.T,v_old);
    msgs = [msgs,compare_action(Y_old,Y_new,tol,'T')];
    % f and A: state derivative
    Y_new = apply_poly(PIE_poly.f,values_by_name(PIE_poly.f,PIE_poly.x,X));
    Y_new = order_as_opvar(Y_new,x_tab);
    Y_old = apply_opvar(PIE_lin.A,v_old);
    msgs = [msgs,compare_action(Y_old,Y_new,tol,'f=A*x')];
    % h and C1, g and C2: regulated and observed outputs
    outs = {'gz','C1','z_tab';  'gy','C2','y_tab'};
    for k=1:size(outs,1)
        Pold = PIE_lin.(outs{k,2});
        if size(Pold,1)==0 || isempty(PIE_poly.(outs{k,1}).C.ops)
            continue
        end
        tab = PIE_poly.(outs{k,3});
        Y_new = apply_poly(PIE_poly.(outs{k,1}),values_by_name(PIE_poly.(outs{k,1}),PIE_poly.x,X));
        Y_new = order_as_opvar(Y_new,tab);
        Y_old = apply_opvar(Pold,v_old);
        msgs = [msgs,compare_action(Y_old,Y_new,tol,[outs{k,1},'=',outs{k,2},'*x'])]; %#ok<AGROW>
    end

    % Check the conversion of the linear PIE by PIE2polyPIE
    PIE_conv = PIE2polyPIE(PIE_lin);
    [X2,v_old2] = random_states(PIE_conv.x_tab,s);
    Y_new = order_as_opvar(apply_op(PIE_conv.T,X2),PIE_conv.x_tab);
    msgs = [msgs,compare_action(apply_opvar(PIE_lin.T,v_old2),Y_new,tol,'PIE2polyPIE: T')];
    Y_new = apply_poly(PIE_conv.f,values_by_name(PIE_conv.f,PIE_conv.x,X2));
    Y_new = order_as_opvar(Y_new,PIE_conv.x_tab);
    msgs = [msgs,compare_action(apply_opvar(PIE_lin.A,v_old2),Y_new,tol,'PIE2polyPIE: f=A*x')];
    nchecked = nchecked+4;

    if isempty(msgs)
        fprintf("Example %2d: no discrepancies\n",example_num);
    else
        nfail = nfail+1;
        fprintf("Example %2d: %s\n",example_num,strjoin(string(msgs),'; '));
    end
end
if nfail>0
    error("Test_Conversion_Linear: discrepancies in %d of %d examples.",nfail,numel(example_nums))
end
fprintf("\n\n --- No discrepancies were encountered \n\n")


%% ---- local functions ----------------------------------------------------
function [X,v_old] = random_states(x_tab,s)
% Random values for the state components defined by x_tab, as a cell with the
% value of each component (X, in the order of x_tab), and as a single column
% in the order of an 'opvar' (v_old, finite-dimensional components first)
n = size(x_tab,1);
X = cell(1,n);
is_dist = logical(x_tab(:,3));
for i=1:n
    if is_dist(i)
        X{i} = rand_poly([x_tab(i,2),1],s,2,4);
    else
        X{i} = rand(x_tab(i,2),1)-0.5;
    end
end
v_old = [];
for i=old_order(x_tab)'
    v_old = [v_old;X{i}];                                                   %#ok<AGROW>
end
end

function xvals = values_by_name(P,x,X)
% Values of the states of the polynomial P, in the order of P.varname. The
% kth variable in x = PIE.x is the kth state component in x_tab (the columns of
% T, in the order of x_tab), so the value X{k} is assigned to the variable x.varname{k}
xvals = cell(numel(P.varname),1);
for k=1:numel(P.varname)
    idx = find(strcmp(x.varname,P.varname{k}));
    xvals{k} = X{idx};
end
end
function idcs = old_order(tab)
% Order of the state components (rows of tab) in an 'opvar': the finite-
% dimensional components first, then the infinite-dimensional ones, both in the
% order of the rows of tab (as in the PIE of the linear converter)
is_dist = logical(tab(:,3));
idcs = [find(~is_dist);find(is_dist)];
end
function Y = order_as_opvar(Y,tab)
% Reorder the stacked output Y of the components defined by tab, from the
% order of the rows of tab to the order in an 'opvar' (see old_order)
sizes = tab(:,2);
if isempty(Y) || size(Y,1)~=sum(sizes)
    return
end
nnr = [0;cumsum(sizes(:))];
idcs = [];
for i=old_order(tab)'
    idcs = [idcs,nnr(i)+1:nnr(i+1)];                                        %#ok<AGROW>
end
Y = Y(idcs,:);
end
function msgs = compare_action(Y_old,Y_new,tol,label)
% Compare the polynomials Y_old and Y_new, returning a message if they differ
msgs = {};
Y_old = polynomial(Y_old);      Y_new = polynomial(Y_new);
if ~isequal(size(Y_old),size(Y_new))
    msgs = {sprintf("%s: output sizes differ ([%s] vs [%s])",label,num2str(size(Y_old)),num2str(size(Y_new)))};
    return
end
if isempty(Y_old),  return,     end
D = Y_old-Y_new;
err = full(max(abs(D.coefficient(:))));
nrm = full(max([1;abs(Y_old.coefficient(:))]));
if err>tol*nrm
    msgs = {sprintf("%s: outputs differ (by %.2g)",label,err)};
end
end


%{
% % Check if the nonlinear converter returns the same PIE as linear
% % converter for linear examples
clear
example_num = 29;       % 21, 31 and 32 have states that do not contribute

% Convert to PIE using both functions
PDE = examples_PDE_library_PIETOOLS_nouinput(example_num);
PIE_lin = convert(PDE);
PIE_poly = convert_PIETOOLS_PDE_nonlinear(PDE);

% Check that the operators match
% First check the T operator
if ~all(all(PIE_lin.T==PIE_poly.T)) || ~all(all(PIE_lin.Tw==PIE_poly.Tw)) || ~all(all(PIE_lin.Tw==PIE_poly.Tw))
    error("One of the T operators does not match")
end
% Next, check the A operator
Aop_poly = 0*PIE_lin.T;
Aop_poly = Aop_poly(:,end+1:end);
for j=1:numel(PIE_poly.f.C.ops)
    Aop_poly = [Aop_poly,PIE_poly.f.C.ops{j}];
end
if ~all(all(PIE_lin.A==Aop_poly))
    error("A operator does not match")
end
% Check the B operators
if ~all(all(PIE_lin.B1==PIE_poly.B1)) || ~all(all(PIE_lin.B2==PIE_poly.B2))
    error("One of the B operators does not match")
end
% Check the C operators
if size(PIE_lin.C1,1)>0
    C1op_poly = 0*PIE_lin.C1;
    C1op_poly = C1op_poly(:,end+1:end);
    for j=1:numel(PIE_poly.gz.C.ops)
        C1op_poly = [C1op_poly,PIE_poly.gz.C.ops{j}];
    end
    if ~all(all(PIE_lin.C1==C1op_poly))
        error("C1 operator does not match")
    end
end
if size(PIE_lin.C2,1)>0
    C2op_poly = 0*PIE_lin.C2;
    C2op_poly = C2op_poly(:,end+1:end);
    for j=1:numel(PIE_poly.gy.C.ops)
        C2op_poly = [C2op_poly,PIE_poly.gy.C.ops{j}];
    end
    if ~all(all(PIE_lin.C2==C2op_poly))
        error("C2 operator does not match")
    end
end
% Check the D operators
if ~all(all(PIE_lin.D11==PIE_poly.D11)) || ~all(all(PIE_lin.D12==PIE_poly.D12)) ||...
        ~all(all(PIE_lin.D21==PIE_poly.D21)) || ~all(all(PIE_lin.D22==PIE_poly.D22))
    error("One of the D operators does not match")
end

fprintf("\n\n --- No discrepancies were encountered \n\n")
%}
