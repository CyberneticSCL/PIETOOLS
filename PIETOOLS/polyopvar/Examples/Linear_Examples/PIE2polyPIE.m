function PIE_out = PIE2polyPIE(PIE_in)
% PIE_OUT = PIE2POLYPIE(PIE_in) takes a linear PIE in 'pie_struct' format
% and returns a structure representing the same PIE but now with the
% right-hand side expressed as a (linear) distributed polynomial
%
% INPUTS
% - PIE_in:     'pie_struct' object representing a linear PIE;
% 
% OUTPUS
% - PIE_out:    'struct' representing the same PIE as the input, but now
%               with the right-hand side of the PIE defined by a
%               'polyopvar' object, PIE_out.f;

Aop = PIE_in.A;
op_dim = Aop.dim(:,1);
if any(op_dim(1:end-1))
    error("PIEs with coupled states on different domains are currently not supported.")
end

% Declare a vector of state variables
nx = size(Aop,1);
x_name = [repmat('x',[nx,1]),num2str((1:nx)')];
x_name = mat2cell(x_name,ones(nx,1),size(x_name,2));

% Declare the right-hand side of the PDE
f = polyopvar(x_name,Aop.var1,Aop.dom);
for i=1:nx
    for j=1:nx
        f.C.ops{i,j} = dopvar2ndopvar(Aop(i,j));
    end
end

% Declare the PIE
PIE_out = struct();
PIE_out.T = dopvar2ndopvar(PIE_in.T);
PIE_out.f = f;
PIE_out.vars = PIE_in.vars;
PIE_out.dom = PIE_in.dom;

end