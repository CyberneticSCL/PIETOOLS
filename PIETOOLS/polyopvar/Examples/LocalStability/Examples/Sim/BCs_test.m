function [pass, err, tol] = BCs_test(U, BC)
%BCS_TEST Test whether simulation data satisfies specified boundary conditions.
%
%   [pass, err] = BCs_test(U, BC)
%
%   U:
%       Solution matrix with
%           dim 1 = space
%           dim 2 = time
%
%   BC:
%       Structure defining left and right boundary conditions:
%
%       BC.left.type    = 'Dirichlet' or 'Neumann'
%       BC.left.value   = boundary value
%
%       BC.right.type   = 'Dirichlet' or 'Neumann'
%       BC.right.value  = boundary value
%
%       BC.tol          = tolerance (optional, default = 1e-8)
%
%       For Neumann BCs:
%       BC.dx           = spatial grid spacing
%
%   OUTPUTS:
%
%       pass:
%           true if all boundary conditions hold within tolerance
%
%       err:
%           Structure containing the maximum error at each boundary.
%
%       tol:
%           Aceptable tolerance on of test.
%
%   Example:
%
%       BC.left.type   = 'Dirichlet';
%       BC.left.value  = 0;
%       BC.right.type  = 'Dirichlet';
%       BC.right.value = 0;
%       BC.tol         = 1e-8;
%
%       [pass, err] = BCs_test(U, BC);

    % -------------------------------------------------------------
    % Tolerance
    % -------------------------------------------------------------
    if isfield(BC, 'tol')
        tol = BC.tol;
    else
        tol = 1e-8;
    end

    % -------------------------------------------------------------
    % Left boundary
    % -------------------------------------------------------------
    switch lower(BC.left.type)

        case 'dirichlet'

            % u(0,t) = specified value
            left_data = U(1,:);
            left_err  = abs(left_data - BC.left.value);

        case 'neumann'

            if ~isfield(BC, 'dx')
                error('BC.dx must be supplied for Neumann BCs.');
            end

            % First-order approximation to u_x at left boundary
            left_data = (U(2,:) - U(1,:)) / BC.dx;
            left_err  = abs(left_data - BC.left.value);

        otherwise
            error('Unknown left BC type: %s', BC.left.type);
    end

    % -------------------------------------------------------------
    % Right boundary
    % -------------------------------------------------------------
    switch lower(BC.right.type)

        case 'dirichlet'

            % u(L,t) = specified value
            right_data = U(end,:);
            right_err  = abs(right_data - BC.right.value);

        case 'neumann'

            if ~isfield(BC, 'dx')
                error('BC.dx must be supplied for Neumann BCs.');
            end

            % First-order approximation to u_x at right boundary
            right_data = (U(end,:) - U(end-1,:)) / BC.dx;
            right_err  = abs(right_data - BC.right.value);

        otherwise
            error('Unknown right BC type: %s', BC.right.type);
    end

    % -------------------------------------------------------------
    % Maximum errors over all time steps
    % -------------------------------------------------------------
    err.left  = max(left_err);
    err.right = max(right_err);

    % -------------------------------------------------------------
    % Pass/fail
    % -------------------------------------------------------------
    pass = (err.left <= tol) && (err.right <= tol);

end