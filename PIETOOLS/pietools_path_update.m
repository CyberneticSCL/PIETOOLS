% This script updates the MATLAB path to include all PIETOOLS directories
% and subdirectories. Save your path after running to keep directories on
% path for future sessions. Don't move this file!
%
% INPUT
% - none; the script adds the folder of this file and all its subfolders,
%   whatever the current folder.
%
% OUTPUT
% - none; the MATLAB path of the current session is extended.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding DJ - 01/06/2025
% MMP, 11/19/2026: Release comment pass (comment-only).
% MMP, 11/21/2026: Add the folder of this file, not pwd: called by name from
%                   another folder, the script added that folder tree and no
%                   PIETOOLS folder.

% Add the folder of this file and all its subfolders to the path.           % MMP, 11/21/2026
% current_folder = pwd;                                                     % MMP, 11/21/2026 (was)
current_folder = fileparts(mfilename('fullpath'));                          % MMP, 11/21/2026
addpath(genpath(current_folder));