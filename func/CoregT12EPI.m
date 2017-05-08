function [] = CoregT12EPI(t1file,epifile)

% This function will perform the following spatial pre-processing steps.
% 1 - coregistration (rigid body) of EPI to T1 image
%
% -------
% INPUTS:
% -------
%
% t1file 	- a string containing the full path and filename of a t1 image. 
% 			e.g., '/path/to/t1/t1.nii'
% 
% epifile   - a string containing the full path and filename of a 4D input EPI file. 
%           e.g., '/path/to/epi/4d.nii'
%
% -------
% OUTPUTS:
% -------
%
% coregistered epi file with prefix 'c'
%
% =========================================================================

spm('defaults','fmri');
spm_jobman('initcfg')

%-----------------------------------------------------------------------
% Job configuration created by cfg_util
%-----------------------------------------------------------------------
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[epifile,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[t1file,',1']};
matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',matlabbatch);
