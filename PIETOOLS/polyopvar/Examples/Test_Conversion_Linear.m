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
if size(PIE_lin.C1)>0
    C1op_poly = 0*PIE_lin.C1;
    C1op_poly = C1op_poly(:,end+1:end);
    for j=1:numel(PIE_poly.h.C.ops)
        C1op_poly = [C1op_poly,PIE_poly.h.C.ops{j}];
    end
    if ~all(all(PIE_lin.C1==C1op_poly))
        error("C1 operator does not match")
    end
end
if size(PIE_lin.C2)>0
    C2op_poly = 0*PIE_lin.C2;
    C2op_poly = C2op_poly(:,end+1:end);
    for j=1:numel(PIE_poly.g.C.ops)
        C2op_poly = [C2op_poly,PIE_poly.g.C.ops{j}];
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