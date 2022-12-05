% example library batch and terms comparison test
close all; clc; clear;
pvar s theta;
examp_list = [1:30];

for i=examp_list
PDEb = examples_PDE_library_PIETOOLS(i,'batch');
PDEt = examples_PDE_library_PIETOOLS(i,'terms');
PIEb = convert_PIETOOLS_PDE(PDEb);
PIEt = convert_PIETOOLS_PDE(PDEt);
if ~(PIEb==PIEt)
    logvals(end+1) = i;
end
end
