function b = subsasgn(a,L,RHS)



switch L(1).type
    case '.'
        if length(L) == 1
            if strcmp(L(1).subs,'vars') || strcmp(L(1).subs,'dom')
                b = set(a,L(1).subs,RHS);
            elseif isa(RHS,'cell')
                b = set(a,L(1).subs,RHS);
            elseif isa(RHS,'struct')
                Lnew = L;
                Lnew(2).type = '{}';
                Lnew(2).subs = numel(a.(L.subs))+1;
                b = subsasgn(a,Lnew,RHS);
            else
                error(['The field ',L(1).subs,' must be a cell of structs.'])
            end
            return
        else
            % Peform all subsasgn but L(1)
            if strcmp(L(1).subs,'Axx')
                L(1).subs = 'A';
            elseif strcmp(L(1).subs,'Bxu')
                L(1).subs = 'Bu';
            elseif strcmp(L(1).subs,'Bxw')
                L(1).subs = 'Bw';
            elseif strcmp(L(1).subs,'Cyx')
                L(1).subs = 'Cy';
            elseif strcmp(L(1).subs,'Czx')
                L(1).subs = 'Cz';
            end
            param = subsref(a,L(1));
            if ismember(L(1).subs,{'x';'u';'w';'y';'z'})
                tab_name = [L(1).subs,'_tab'];
                state_num = L(2).subs;
                tab_new = a.(tab_name);
                if state_num > size(a.(tab_name),1)
                    tab_new = [tab_new; zeros(1,size(tab_new,2))];
                end
                
                if length(L)==2
                    if ~isfield(RHS,'vars')
                        if strcmp(L(1).subs,'x')
                            varnames = a.vars(:,1).varname;
                        else
                            varnames = {};
                        end
                    end
                    if numel(RHS.vars)>2
                        error('At most 2 spatial variables are supported.')
                    elseif iscellstr(RHS.vars)
                        varnames = RHS.vars(:);
                    elseif ispvar(RHS.vars)
                        varnames = RHS.vars.varname;
                    else
                        error('Variables must be specified as pvar (''polynomial'' class) objects')
                    end
                    
                    [is_oldvar, var_indx] = ismember(varnames,a.vars(:,1).varname);
                    if any(~is_oldvar)
                        error('The variables in the proposed field do not match the current spatial variables.')
                    end
                    tab_new = a.(tab_name);
                    tab_new(state_num,2+var_indx(is_oldvar)) = is_oldvar;
                    
                    if ~isfield(RHS,'diff')
                        % Assume no differentiability if not specified
                        diff = zeros(1,size(varnames));
                    end
                    tab_new(state_num,4+var_indx(is_oldvar)) = diff(is_oldvar);
                    
                end
                    
                a = set(a,tab_name,[a.(tab_name); tab_new]);
                        
                
            elseif ismember(L(1).subs,{'A';'Bu';'Bw';'Cy';'Cz';'Dyu';'Dyw';'Dzu';'Dzw'})
                
                
            elseif ~ismember(L(1).subs,{'dom';'vars'})
                error(['Proposed field ',L(1).subs,' is not settable.']);
            end
                
                
            switch L(2).type
                case '{}'
                    if length(L(2).subs)==1
                        [nr,nc] = size(param);
                        [rindx,cindx] = ind2sub([nr,nc],L(2).subs{1});
                    elseif length(L(2).subs)>=2
                        rindx = L(2).subs{1};   cindx = L(2).subs{2};
                    end
                    if length(L)==2
                        param{rindx,cindx} = RHS;
                    else
                        param{rindx,cindx} = subsasgn(param{rindx,cindx},L(3:end),RHS);
                    end
                otherwise
                    param = subsasgn(param,L(2:end),RHS);
            end
        end
        b = set(a,L(1).subs,param);
        
    case '()'
        error('()- like subsassign is currently not supported for opvar2d objects.');
        
%         if length(L)==1
%             temp = polynomial(RHS);
%         else
%             % Peform all subsasgn but L(1)
%             temp = subsref(a,L(1));
%             temp = subsasgn(temp,L(2:end),RHS);
%         end
%         
%         %  Three '()'-subsasgn cases
%         if length(L(1).subs)==1 && strcmp(L(1).subs{1},':')
%             b = PVsubsasgn_colon(a,L,temp);
%         elseif length(L(1).subs)==1
%             b = PVsubsasgn_1idx(a,L,temp);
%         else
%             if strcmp(L(1).subs{1},':')
%                 % To handle the case a(:,1)=x1 when a is undefined or 0-by-0
%                 if all(sza==[0 0])
%                     sza(1) = 1;
%                 end
%                 L(1).subs{1} = 1:sza(1);
%             end
%             if strcmp(L(1).subs{2},':')
%                 % To handle the case a(1,:)=x1 when a is undefined or 0-by-0
%                 if all(sza==[0 0])
%                     sza(2) = 1;
%                 end
%                 L(1).subs{2} = 1:sza(2);
%             end
%             b = PVsubsasgn_2idx(a,L,temp);
%         end
        
    case '{}'
        error('{}- like subsassign is not supported for opvar2d objects.');
end


end

