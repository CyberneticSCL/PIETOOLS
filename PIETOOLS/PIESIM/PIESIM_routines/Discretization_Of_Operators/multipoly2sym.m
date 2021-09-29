function Rop=multipoly2sym(R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multipoly2sym.m     PIETOOLS 2021a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function converts a pvar object in s and theta to a symbolic object
% in sym_s and sym_theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
% 6/7/2021 - SS removed an unnecessary summation in two pvars case

syms sym_s sym_theta;
assume(sym_s, 'real');
assume(sym_theta, 'real');

if isa(R,'polynomial')
    coeff = R.coefficient;
    dmat = R.degmat;
    switch(R.nvars)
        case 2
        % Two variables present
            symlist = [sym_s.^dmat(:,1) sym_theta.^dmat(:,2)];
            symlist = symlist(:,1).*symlist(:,2);
            Rop = symlist'*coeff;
%             Rop = sum(Rop,1);
            Rop=reshape(Rop,size(R));
        case 1
        % One variable present
            var=cell2mat(R.varname);
            if (var=='s')
                symlist = sym_s.^dmat;
                Rop=symlist'*coeff;
            else
                symlist = sym_theta.^dmat;
                Rop=symlist'*coeff;
            end
            Rop = reshape(Rop,size(R));
        case 0
            % Scalar value
            Rop=reshape(coeff(1,:),size(R));
    end
else
    Rop = R;
end
   