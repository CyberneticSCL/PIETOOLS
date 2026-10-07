function DMB_test()
    % DMB_test() tests the distributed-monomial basis product implemented
    % by DMB.m for the m=1 case of Def. 4.
    %
    % The test constructs the basis factor U*x used by SOS_DP and
    % V3_Liediff, then verifies that DMB recursively constructs
    %
    %       (U*x)^(otimes 2) and (U*x)^(otimes 3)
    %
    % with the expected distributed degrees, component dimensions,
    % metadata, and tensor-PI coefficient operators. The coefficient
    % operators are compared with independently evaluated otimes calls.
    % Representative input-validation errors are also checked.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - DMB_test
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
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    % GNU General Public License for more details.
    %
    % You should have received a copy of the GNU General Public License
    % along with this program; if not, write to the Free Software
    % Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % If you modify this code, document all changes carefully and include date
    % authorship, and a brief description of modifications
    %
    % CR, 10/06/2026: Initial coding.

    % Add PIETOOLS to the path when this test is called directly from its
    % Tests directory or from another working directory.
    test_dir      = fileparts(mfilename('fullpath'));
    local_dir     = fileparts(test_dir);
    pietools_root = fileparts(fileparts(fileparts(local_dir)));
    addpath(genpath(pietools_root));

    clear stateNameGenerator


    %% Construct a scalar fundamental state.

    % Converting this simple Fisher equation provides the same scalar
    % polyopvar state used by the LocalStability examples.
    pvar s t
    dom = [0,1];
    u   = pde_var(s,dom);
    PDE = [diff(u,t)==diff(u,s,2) + 5*u + u^2;
           subs(u,s,dom(1))==0;
           subs(u,s,dom(2))==0];

    PIE = convert(PDE);
    x   = PIE.f.vartab;


    %% Construct one distributed-monomial basis factor U*x.

    % Use the integral-only basis employed when recursively constructing
    % tensor powers in SOS_DP and V3_Liediff.
    pvar theta
    opdeg    = 1;
    Zmon     = monomials([s,theta],0:opdeg);
    Zop      = opvar();
    Zop.var1 = s;
    Zop.var2 = theta;
    Zop.I    = dom;
    Zop.R.R0 = [0*Zmon;0*Zmon];
    Zop.R.R1 = [Zmon;0*Zmon];
    Zop.R.R2 = [0*Zmon;Zmon];

    U  = dopvar2ndopvar(Zop);
    Ux = U*x;

    assert(isa(Ux,'polyopvar') && isequal(size(Ux.C.ops),[1,1]), ...
        'DMB_test:InvalidBasisFactor', ...
        'The test basis U*x must contain one polyopvar coefficient.');
    assert(isa(Ux.C.ops{1},'tensopvar'), ...
        'DMB_test:InvalidBasisCoefficient', ...
        'The coefficient of U*x must be a tensopvar object.');


    %% Check the degree-two product (U*x) otimes (U*x).

    Ux2 = DMB(Ux,Ux);

    assert(isa(Ux2,'polyopvar'), ...
        'DMB_test:InvalidDegreeTwoClass', ...
        'The degree-two product must be a polyopvar object.');
    assert(isequal(Ux2.degmat,Ux.degmat+Ux.degmat), ...
        'DMB_test:IncorrectDegreeTwoDegree', ...
        'DMB must add the distributed degrees of its two inputs.');
    verify_metadata(Ux2,Ux,'degree-two');

    T1      = Ux.C.ops{1};
    T2      = Ux2.C.ops{1};
    T2_true = otimes(T1,T1,[true,true]);
    k1      = size(T1,1);

    assert(isa(T2,'tensopvar'), ...
        'DMB_test:InvalidDegreeTwoCoefficient', ...
        'The degree-two coefficient must be a tensopvar object.');
    assert(size(T2,1)==k1^2, ...
        'DMB_test:IncorrectDegreeTwoDimension', ...
        'The degree-two output dimension must be k^2.');
    verify_tensopvar(T2,T2_true,'degree-two');


    %% Check the recursive degree-three product.

    Ux3 = DMB(Ux2,Ux);

    assert(isequal(Ux3.degmat,Ux2.degmat+Ux.degmat), ...
        'DMB_test:IncorrectDegreeThreeDegree', ...
        'The recursive DMB product must add the distributed degrees.');
    verify_metadata(Ux3,Ux,'degree-three');

    T3      = Ux3.C.ops{1};
    T3_true = otimes(T2,T1,[true,true]);

    assert(isa(T3,'tensopvar'), ...
        'DMB_test:InvalidDegreeThreeCoefficient', ...
        'The degree-three coefficient must be a tensopvar object.');
    assert(size(T3,1)==k1^3, ...
        'DMB_test:IncorrectDegreeThreeDimension', ...
        'The degree-three output dimension must be k^3.');
    verify_tensopvar(T3,T3_true,'degree-three');


    %% Check representative input validation.

    verify_error(@()DMB(1,Ux),'polyopvar_product:InvalidInput');

    Ux_bad = Ux;
    Ux_bad.varname = {'different_state'};
    verify_error(@()DMB(Ux,Ux_bad), ...
        'polyopvar_product:IncompatibleVariables');

    Ux_bad = Ux;
    Ux_bad.dom = [0,2];
    verify_error(@()DMB(Ux,Ux_bad), ...
        'polyopvar_product:Incompatible_dom');

    Ux_bad = Ux;
    Ux_bad.degmat = [1;1];
    verify_error(@()DMB(Ux,Ux_bad), ...
        'polyopvar_product:Invalid_degmat');

    fprintf(['\n --- DMB tests passed: degree addition, recursive tensor ', ...
             'products, metadata, dimensions, and input validation ---\n\n']);

end


function verify_metadata(C,A,degree_label)
    % Verify metadata which DMB should preserve from its compatible inputs.

    assert(isequal(C.varname,A.varname), ...
        'DMB_test:IncorrectVariableNames', ...
        'DMB did not preserve state names in the %s product.',degree_label);
    assert(isequal(C.varsize,A.varsize(:)), ...
        'DMB_test:IncorrectVariableSize', ...
        'DMB did not preserve state dimensions in the %s product.',degree_label);
    assert(isequal(C.pvarname,A.pvarname), ...
        'DMB_test:IncorrectSpatialVariables', ...
        'DMB did not preserve spatial variables in the %s product.',degree_label);
    assert(isequal(C.dom,A.dom), ...
        'DMB_test:IncorrectDomain', ...
        'DMB did not preserve the spatial domain in the %s product.',degree_label);
    assert(isequal(C.varmat,A.varmat), ...
        'DMB_test:IncorrectVariableMap', ...
        'DMB did not preserve varmat in the %s product.',degree_label);
    assert(isequal(size(C.C.ops),[1,1]), ...
        'DMB_test:IncorrectCoefficientCount', ...
        'DMB must retain one coefficient in the %s product.',degree_label);
end


function verify_tensopvar(actual,expected,degree_label)
    % Compare the DMB coefficient with a direct Def. 4 otimes evaluation.

    % Separate calls to otimes introduce fresh dummy-variable names. Thus,
    % comparing actual.ops or actual.vars with isequal would reject two
    % mathematically identical tensor products merely because their dummy
    % variables have different names. Compare the name-independent tensor
    % invariants instead: total dimensions, factor dimensions/domains,
    % dependency maps, product type, and factor ordering.
    assert(isequal(size(actual),size(expected)), ...
        'DMB_test:IncorrectTensorDimension', ...
        'Incorrect total tensor dimensions in the %s product.',degree_label);
    assert(isequal(size(actual.ops),size(expected.ops)), ...
        'DMB_test:IncorrectTensorFactors', ...
        'Incorrect number of tensor factors in the %s product.',degree_label);
    assert(isequal(actual.dims,expected.dims), ...
        'DMB_test:IncorrectFactorDimensions', ...
        'Incorrect tensor-factor dimensions in the %s product.',degree_label);
    assert(isequal(actual.doms,expected.doms), ...
        'DMB_test:IncorrectFactorDomains', ...
        'Incorrect tensor-factor domains in the %s product.',degree_label);
    assert(isequal(size(actual.vars),size(expected.vars)), ...
        'DMB_test:IncorrectTensorVariableCount', ...
        'Incorrect number of tensor variables in the %s product.',degree_label);
    assert(isequal(actual.dom,expected.dom), ...
        'DMB_test:IncorrectTensorDomain', ...
        'Incorrect tensor domain in the %s product.',degree_label);
    assert(isequal(actual.depmat1,expected.depmat1), ...
        'DMB_test:IncorrectOutputDependencies', ...
        'Incorrect output dependencies in the %s product.',degree_label);
    assert(isequal(actual.depmat2,expected.depmat2), ...
        'DMB_test:IncorrectInputDependencies', ...
        'Incorrect input dependencies in the %s product.',degree_label);
    assert(isequal(actual.type,expected.type), ...
        'DMB_test:IncorrectTensorType', ...
        'Incorrect tensor-product type in the %s product.',degree_label);
    assert(isequal(actual.order,expected.order), ...
        'DMB_test:IncorrectTensorOrder', ...
        'Incorrect tensor-factor order in the %s product.',degree_label);
end


function verify_error(test_fcn,expected_id)
    % Verify that TEST_FCN throws the requested DMB input-validation error.

    did_error = false;
    try
        test_fcn();
    catch ME
        did_error = true;
        assert(strcmp(ME.identifier,expected_id), ...
            'DMB_test:IncorrectErrorIdentifier', ...
            'Expected error "%s", but received "%s".', ...
            expected_id,ME.identifier);
    end

    assert(did_error, ...
        'DMB_test:ExpectedErrorNotThrown', ...
        'Expected error "%s" was not thrown.',expected_id);
end
