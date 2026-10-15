function [pass, err, tol] = BCs_test(U, BC)
% [PASS,ERR,TOL] = BCS_TEST(U,BC) tests whether the simulation data U
% satisfies the boundary conditions specified by BC.
%
% INPUTS
% - U:      array with the solution, where dimension 1 is space and
%           dimension 2 is time;
% - BC:     struct defining the left and right boundary conditions, with
%           fields
%               'left':     struct with fields 'type' ('Dirichlet' or
%                           'Neumann') and 'value' (boundary value);
%               'right':    struct with fields 'type' ('Dirichlet' or
%                           'Neumann') and 'value' (boundary value);
%               'tol':      (optional) tolerance. Defaults to 1e-8;
%               'dx':       spatial grid spacing, required for Neumann
%                           BCs;
%
% OUTPUTS
% - pass:   logical, true if all boundary conditions hold within the
%           tolerance;
% - err:    struct with fields 'left' and 'right', specifying the maximum
%           error at each boundary over all time steps;
% - tol:    acceptable tolerance of the test;
%
% NOTES
% - Neumann BCs are tested using a first-order finite-difference
%   approximation of u_x at the boundary.
%
% EXAMPLE
%   BC.left.type   = 'Dirichlet';
%   BC.left.value  = 0;
%   BC.right.type  = 'Dirichlet';
%   BC.right.value = 0;
%   BC.tol         = 1e-8;
%
%   [pass, err] = BCs_test(U, BC);

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