

ntests = 100;

dim = 1;
max_diff = 3;
n_comps = randi([1,4],1,5);
max_size = randi([1,4],1,5);
max_terms = 5;
max_deg = 2;
add_BC_terms = 1;

tol = 1e-9;

for test_num=1:ntests

    % Generate a random PDE.
    PDE = PIETOOLS_random_PDE_generator(dim,max_diff,n_comps,max_size,max_terms,max_deg,add_BC_terms);
    %assignin('base','PDE_test',PDE);
    % Convert the PDE to a PIE.
    if ~PDE.dim
        continue;
    end
    PDE = initialize_PIETOOLS_PDE_terms(PDE);
    PIE = convert_PIETOOLS_PDE(PDE);
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