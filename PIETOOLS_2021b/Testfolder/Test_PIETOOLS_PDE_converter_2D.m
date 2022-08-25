

ntests = 100;

dim = 2;
max_diff = 3;
n_comps = [3,2,2,2,2];
max_size = [2,2,2,2,2];
max_terms = 5;
max_deg = 2;
add_BC_terms = 0;

tol = 1e-9;

for test_num=1:ntests

    % Generate a random PDE.
    PDE = PIETOOLS_random_PDE_generator(dim,max_diff,n_comps,max_size,max_terms,max_deg,add_BC_terms);
    %assignin('base','PDE_test',PDE);
    % Convert the PDE to a PIE.
    PIE = convert_PIETOOLS_PDE_2D(PDE);
    %assignin('base','PIE_test',PIE);

    % Check that the PIE is appropriate.
    ismatch_arr = PIETOOLS_compare_PDE_PIE(PDE,PIE,tol);
    % Stop the script if the PDE and PIE do not match.
    if any(~ismatch_arr)
        txt = input('Would you like to continue? (y/n) \n',"s");
        if strcmpi(txt,'n') ||  strcmpi(txt,'no')
            break
        else
            continue
        end
    end


end