

ntests = 100;

dim = 1;
max_diff = 3;
n_comps = randi([1,4],1,5);
n_comps = [randi(4),1,1,1,1];
max_size = randi([1,4],1,5);
max_size = [1,randi([1,4],1,4)];
max_terms = 10;
max_deg = 2;
add_BC_terms = 0;

tol = 1e-9;

% for test_num=1:ntests
% 
%     % Generate a random PDE in terms format.
%     PDE_t = PIETOOLS_random_PDE_generator(dim,max_diff,n_comps,max_size,max_terms,max_deg,add_BC_terms);
%     % Declare in terms format
%     PDE_p = PIETOOLS_terms2parser(PDE_t);
%     PDE_params = initialize(PDE_p.params,true);
%     % Manual correction: u or w may appear in PDE_t without actually
%     % contributing to any terms
%     if isempty(PDE_params.u) || isempty(PDE_params.u{1})
%         PDE_t.u = PDE_params.u;
%     end
%     if isempty(PDE_params.w) || isempty(PDE_params.w{1})
%         PDE_t.w = PDE_params.w;
%     end
%     
%     % Compare terms to parser format
%     logval = isequal_PDEs_terms(PDE_t,PDE_params);
%     
%     if ~logval
%         txt = input('Would you like to continue? (y/n) \n',"s");
%         if strcmpi(txt,'n') ||  strcmpi(txt,'no')
%             break
%         else
%             continue
%         end
%     end
% 
% end


for test_num=1:29
    % Generate a random PDE in terms format.
    PDE_t = examples_PDE_library_PIETOOLS(test_num,'terms');
    % Declare in terms format
    PDE_p = PIETOOLS_terms2parser(PDE_t);
    PDE_params = initialize(PDE_p.params,true);
    
    % Compare terms to parser format
    logval = isequal_PDEs_terms(PDE_t,PDE_params);
    
    if ~logval
        txt = input('Would you like to continue? (y/n) \n',"s");
        if strcmpi(txt,'n') ||  strcmpi(txt,'no')
            break
        else
            continue
        end
    end

end