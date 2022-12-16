% Installation/Update script for PIETOOLS
% This script downloads the latest version of PIETOOLS, with SOSTOOLS, 
% Multipoly and dpvar, and adds to the path, overriding
% previous installations. 
% The script also installs latest stable version of SeDuMi overriding
% previous installations.
% Updated - 12/15/2022

clc;
disp('----------------------------------------------');
disp('Update PIETOOLS to latest version on GitHub.');
disp('----------------------------------------------');

current_dir = pwd;

% % % Install latest version of PIETOOLS (with SOSTOOLS included)
if exist('PIETOOLS_PDE','file')
    % % PIETOOLS has already been downloaded.

    % Move to directory containing current version of PIETOOLS.
    old_dir = fileparts(which('PIETOOLS_PDE'));
    cd(old_dir);
    cd ..

    if exist('PIETOOLS','dir')
        fprintf('\n Your current version of PIETOOLS will be overwritten. \n');
        msg = ['    Would you like to proceed? [y/n] \n ---> '];
        resp = input(msg,'s');
        resp = strrep(resp,'''','');
        if ~(strcmpi(resp,'y') || strcmpi(resp,'yes'))
            fprintf('\n Update aborted. Rename your current PIETOOLS folder and try updating again. \n');
            return
        end
    end
    download_dir = pwd;

    % Download PIETOOLS in same directory as previous install.
    disp(' ')
    disp('Downloading PIETOOLS from GitHub...');
    f = websave('PIETOOLS.zip','https://github.com/CyberneticSCL/PIETOOLS/zipball/master/');
    [~,~,f3] = fileparts(f);
    if strcmp(f3,'.html')
        error('Could not download PIETOOLS from the internet. Check internet access or website status.');
    end
    rehash;
    
    % Unzip PIETOOLS.
    disp("Extracting PIETOOLS files..");
    filenames = unzip('PIETOOLS.zip');
    delete('PIETOOLS.zip');
    
    % Enter PIETOOLS directory.
    cd(filenames{1});

    % Move files from "shell" folder into PIETOOLS code folder.
    movefile('LICENSE',['PIETOOLS',filesep,'LICENSE'])
    movefile('Readme.m',['PIETOOLS',filesep,'Readme.m'])
    movefile('README.md',['PIETOOLS',filesep,'README.md'])

    % Extract PIETOOLS code folder from "shell" folder.
    cd ..
    is_moved = movefile([filenames{1},filesep,'PIETOOLS']);
    if ~is_moved
        error(['The folder',[filenames{1},filesep,'PIETOOLS'],' containing the PIETOOLS code could not be moved.',...
               ' Please manually move this folder to the desired location.'])
    end

    % Remove the "shell" folder
    is_removed = rmdir(filenames{1},'s');
    if ~is_removed
        rm_path = [pwd,filesep,filenames{1}];
        error(['The folder "',rm_path,'" was not removed. Please manually remove this folder.'])
    end

else
    % % Fresh install of PIETOOLS 
    disp(' ');
    disp('No version of PIETOOLS was encountered; downloading the latest version...')

    % Allow user to select download directory.
    disp(' ');
    fprintf(['Choose the installation directory.\n',...
      'A folder "PIETOOLS" is going to be created in the specified location.\n',...
      'If not specified, the current directory will be chosen as installation directory.\n']);

    default_dir = pwd;
    download_dir = uigetdir(pwd);
    if isequal(download_dir,0)
        fprintf(['No directory selected.\n',... 
            'Installing in the current directory "%s"?\n'],default_dir);
        download_dir = default_dir;
    end

    % Move to download directory and download.
    cd(download_dir)    
    disp('Downloading PIETOOLS from GitHub...');
    f = websave('PIETOOLS.zip','https://github.com/CyberneticSCL/PIETOOLS/zipball/master/');
    [~,~,f3] = fileparts(f);
    if strcmp(f3,'.html')
        error('Could not download PIETOOLS from the internet. Check internet access or website status.');
    end
    rehash;

    % Unzip PIETOOLS, and rename folder to PIETOOLS
    disp("Extracting PIETOOLS files...");
    filenames = unzip('PIETOOLS.zip');
    delete('PIETOOLS.zip');
    movefile(filenames{1},'PIETOOLS');
    cd('PIETOOLS');
end

% Determine path at which PIETOOLS is located.
toolbox_path = [pwd,filesep,'PIETOOLS'];
cd(toolbox_path)

  
% % % Install latest version of SOSTOOLS

% Check if tbxmanager is already installed
if ~exist('tbxmanager')
    % Create a new directory in PIETOOLS folder
    d = [toolbox_path,filesep,'tbxmanager'];
    if isequal(exist(d,'dir'),7)
        error('The installation directory "%s" already exists.\nPlease, remove or rename the folder or change the installation path.',d);
    end
    disp('Creating the directory "tbxmanager".');
    out = mkdir(d);
    if ~out
        error(['An error appeared when trying to create the folder "%s".\n',...
              'Please, check administrative access to modify the destination directory.'],c); 
    end
    
    % Enter tbxmanager directory
    cd(d);
    
    % Download tbxmanager
    disp(' ');
    disp('Downloading the Toolbox Manager from the internet...');
    [~, c] = urlwrite('http://www.tbxmanager.com/tbxmanager.m', 'tbxmanager.m');
    if ~c
        fprintf(2,['\n Could not download tbxmanager from the internet.\n',...
                    ' Please manually install an SDP solver (e.g. SeDuMi) if not already installed.\n']);
        return
    end
    rehash;
end

% Create the initialization file to set the path 
disp(' ');
disp('Creating the initialization file "startup.m"...');
p = which('startup.m');
d = which('tbxmanager'); d = strrep(d,'tbxmanager.m','');
if isempty(p)
    p = [d,'startup.m'];
end
fid = fopen(p,'a');
if isequal(fid,-1)
    error(['Could not modify the initialization file "startup.m".',...
           'Edit this file in the folder "%s" manually and insert there the line:  tbxmanager restorepath.'],p);
end
filestrings = fileread(p);
pathexist = contains(filestrings,'tbxmanager restorepath');
if ~pathexist
    fprintf(fid,'tbxmanager restorepath\n');
    fclose(fid);
    disp('File has been created.');
end

% Download SeDuMi
disp(' ');
disp('Downloading and installing SeDuMi...');
tbxmanager install sedumi 
rehash;


% % % Update the path.

% Add path to PIETOOLS + SeDuMi
disp(' ');
disp('Adding PIETOOLS path to MATLAB...');
addpath(genpath(toolbox_path));

% Save path for future sessions
disp(' ');
disp('Saving path for future sessions...');
status = savepath;

% Check if the path was properly saved.
if status
    fprintf('Could not save the path to a default location,\nplease provide a location where you want to save the path.');
    cn = uigetdir(pwd);
    if isequal(cn,0)
        disp(' ');
        fprintf('No directory specified, saving the path to the current directory "%s".\n\n',default_dir);
        cn = default_dir;
    end
    sn = savepath([cn,filesep,'pathdef.m']);
    if sn
        error(['Could not save the path automatically.\n',...
            'Please, open the "Set Path" button in the Matlab menu and save the path manually to some location.']);
    end
end

disp(' ');
disp('Installation finished.');
cd(current_dir) % Return to starting directory