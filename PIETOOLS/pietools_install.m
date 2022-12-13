% Installation Script for PIETOOLS
% This script adds SOSTOOLS, MULTIPOLY and PIETOOLS to the path overriding
% previous installations. 
% The script also installs latest stable version of SeDuMi overriding
% previous installations.
% The script does not delete any previous versions of sostools and sedumi, 
% but just changes path priority.
% Updated - 12/30/2021

clc;
disp('----------------------------------------------');
disp('Installation of PIETOOLS.');
disp('----------------------------------------------');
disp(' ');
fprintf(['Choose the installation directory.\n',...
      'A folder "PIETOOLS" is going to be created in the specified location.\n',...
      'If not specified current directory will be chosen as the installation directory/\n']);
  

% get the installation folder
default_dir = pwd;
c = uigetdir(pwd);
if isequal(c,0)
    fprintf(['No directory selected.\n',... 
        'Installing in the current directory "%s"?\n'],default_dir);
    c = default_dir;
end

c_p = [c,filesep,'PIETOOLS_2021b'];
out = mkdir(c_p);
if ~out
    error(['An error appeared when trying to create the folder "%s".\n',...
          'Please, check administrative access to modify the destination directory.'],c_p); 
end


%check if tbxmanager is already installed
if ~exist('tbxmanager')
% create a new directory in that folder
d = [c_p,filesep,'tbxmanager'];
if isequal(exist(d,'dir'),7)
    error('The installation directory "%s" already exists.\nPlease, remove or rename the folder or change the installation path.',d);
end
disp('Creating the directory "tbxmanager".');
out = mkdir(d);
if ~out
    error(['An error appeared when trying to create the folder "%s".\n',...
          'Please, check administrative access to modify the destination directory.'],c); 
end

% enter that directory
cd(d);

%download tbxmanager
disp(' ');
disp('Downloading the Toolbox manager from the internet.');
[f, c] = urlwrite('http://www.tbxmanager.com/tbxmanager.m', 'tbxmanager.m');
rehash;
end

% create the initialization file to set the path 
disp(' ');
disp('Creating the initialization file "startup.m".');
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

% download SeDuMi
disp(' ');
disp('Downloading and installing SeDuMi');
tbxmanager install sedumi 
rehash;



% get back to the original directory
cd(default_dir);

% % add path to tbxmanager
% disp(' ');
% disp('Adding path to Matlab.');
% addpath(c_p);


cd(c_p);

% Downloading SOSTOOLS, MULTIPOlY and PIETOOLS  %%%Change Path to PIETOOLS zip%%%

disp(' ');
disp('Downloading the PIETOOLS from control.asu.edu.');
[f, c] = urlwrite('http://control.asu.edu/pietools/PIETOOLS_2021b.zip', 'PIETOOLS_2021b.zip');
rehash;



if isequal(c,0)
    error('Could not download PIETOOLS from the internet. Check internet access or website status.');
end
% unzip PIETOOLS
disp("Extracting PIETOOLS files");
unzip('PIETOOLS_2021b.zip');
delete PIETOOLS_2021b.zip;

% get back to the original directory
cd(default_dir);

% add path to PIETOOLS+SeDuMi
disp(' ');
disp('Adding PIETOOLS path to MATLAB.');
addpath(genpath(c_p));

% save path for future
disp(' ');
disp('Saving path for future sessions.');
status = savepath;

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