function [Pcell_filtered_full, Pcell_filtered_nonzero] = filter_Pcell(Pcell, threshold)
    % [...] = filter_Pcell(...) removes small
    % coefficients from the solved Gram-matrix blocks in Pcell.
    %
    % Each polynomial coefficient whose absolute value is strictly less
    % than threshold is set equal to zero. Polynomial terms which become
    % zero are then removed using cleanpoly. The dimensions and cell-array
    % structure of Pcell are preserved.
    %
    % INPUTS
    % - Pcell       Cell array containing solved polynomial Gram-matrix
    %               blocks, such as the solved Pcell returned by
    %               V3_Liediff. Numeric blocks are also supported.
    % - threshold   Finite, nonnegative scalar specifying the absolute
    %               coefficient threshold. Coefficients with magnitude
    %               equal to threshold are retained.
    %
    % OUTPUTS
    % - Pcell_filtered_full  Filtered copy of Pcell with all coefficients of
    %                   magnitude below threshold set to zero.
    % - Pcell_filtered_nonzero Structure containing the size, values and positions
    %                   of the non-zero values in Pcell_filtered_full.
    %
    % NOTES
    % - Pcell should first be evaluated at the PIESOS solution using
    %
    %       Pcell_sol = piesos_getsol(prog_sol,Pcell);
    %
    %   Filtering an unsolved dpvar block is not supported because its
    %   entries still depend on unknown decision variables.
    % - Filtering is intended only to improve readability. A threshold
    %   that is too large may materially alter the recovered certificate.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - filter_Pcell.m
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
    % authorship, and a brief description of modifications.
    %
    % CR, 10/07/2026: Initial coding.

    narginchk(2,2);

    if ~iscell(Pcell)
        error('filter_Pcell:InvalidPcell', ...
            'Pcell must be provided as a cell array.');
    end

    if ~isnumeric(threshold) || ~isscalar(threshold) || ...
            ~isreal(threshold) || ~isfinite(threshold) || threshold < 0
        error('filter_Pcell:InvalidThreshold', ...
            'threshold must be a finite, nonnegative, real scalar.');
    end

    Pcell_filtered_full = Pcell;

    for block_idx = 1:numel(Pcell)
        block = Pcell{block_idx};

        if isempty(block)
            continue
        elseif isa(block,'polynomial')
            % cleanpoly filters every scalar coefficient in the polynomial
            % matrix and removes monomials whose coefficients become zero.
            Pcell_filtered_full{block_idx} = cleanpoly(block,threshold);
        elseif isnumeric(block)
            block(abs(block)<threshold) = 0;
            Pcell_filtered_full{block_idx} = block;
        elseif isa(block,'dpvar')
            error('filter_Pcell:UnsolvedDecisionVariable', ...
                ['Pcell{%d} is a dpvar block. Evaluate Pcell using ', ...
                 'piesos_getsol before applying filter_Pcell.'],block_idx);
        else
            error('filter_Pcell:UnsupportedBlock', ...
                'Pcell{%d} has unsupported class "%s".', ...
                block_idx,class(block));
        end
    end

    % Extract and store non-zero elements of Pcell_filtered and there positions. 
    Pcell_filtered = Pcell_filtered_full{1,1};

    positions = zeros(0,2);
    values    = {};
    
    for row_idx = 1:size(Pcell_filtered,1)
        for col_idx = 1:size(Pcell_filtered,2)
            entry = Pcell_filtered(row_idx,col_idx);
    
            % The entry is nonzero if at least one polynomial coefficient
            % remains nonzero after filtering.
            is_nonzero = any(entry.coefficient(:) ~= 0);
    
            if is_nonzero
                positions(end+1,:) = [row_idx,col_idx];
                values{end+1,1}    = entry;
            end
        end
    end
    
    Pcell_filtered_nonzero.size      = size(Pcell_filtered);
    Pcell_filtered_nonzero.positions = positions;
    Pcell_filtered_nonzero.values    = values;


end
