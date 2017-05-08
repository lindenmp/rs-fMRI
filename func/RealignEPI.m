function [] = RealignEPI(epifile, N)

% This function will perform the following spatial pre-processing steps.
% 1 - realignment of fMRI time series
%
% -------
% INPUTS:
% -------
%
% epifile   - a string containing the full path and filename of a 4D input EPI file. 
% 			e.g., '/path/to/epi/4d.nii'
%
% N       - Length of time series (no of volumes in EPI 4D file).
%
% -------
% OUTPUTS:
% -------
%
% realigned epi file with prefix 'r'
%
% =========================================================================

spm('defaults','fmri');
spm_jobman('initcfg')

%-----------------------------------------------------------------------
% Job configuration created by cfg_util
%-----------------------------------------------------------------------
%%


% define epi list
for i = 1:N
    matlabbatch{1}.spm.spatial.realign.estwrite.data{1}{i,1} = [epifile,',',num2str(i)];
end

% Realignment
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = {''};
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

spm_jobman('run',matlabbatch);
