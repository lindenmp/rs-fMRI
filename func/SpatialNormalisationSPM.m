function [] = SpatialNormalisationSPM(epifile, t1file, mni_template, voxdim, N)

% This function will perform the following spatial pre-processing steps using SPM's methods
% 1 - realignment of fMRI time series
% 2 - coregistration of EPI to T1
% 3 - spatial normalization of T1 to MNI
% 4 - Application of the normalization parameters to coregistered EPI
%
% -------
% INPUTS:
% -------
%
% epifile   - a string containing the full path and filename of a 4D input EPI file. 
%           e.g., '/path/to/epi/4d.nii'
%
% t1file   - path and filename of t1 to coregister to.
%
% mni_template - path and filename of template used for spatial
%               normalization
% 
% voxdim  - scalar indicating the desired voxel dimensions (in mm) 
%	    of the normalised output. e.g, voxdim = 2 will result in output 
%	    with 2 mm^3 voxels.
%
% N       - Length of time series (no of volumes in EPI 4D file).
%
% -------
% OUTPUTS:
% -------
%
% realigned files will be output with the prefix 'r', coregistered with
% 'cr' and spatially normalized with 'w'.
%
% =========================================================================
spm('defaults','fmri');
spm_jobman('initcfg')

%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
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

% Coregistration of EPI to T1
matlabbatch{2}.spm.spatial.coreg.estimate.ref = {[t1file,',1']};
matlabbatch{2}.spm.spatial.coreg.estimate.source(1) = cfg_dep;
matlabbatch{2}.spm.spatial.coreg.estimate.source(1).tname = 'Source Image';
matlabbatch{2}.spm.spatial.coreg.estimate.source(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{2}.spm.spatial.coreg.estimate.source(1).tgt_spec{1}(1).value = 'image';
matlabbatch{2}.spm.spatial.coreg.estimate.source(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{2}.spm.spatial.coreg.estimate.source(1).tgt_spec{1}(2).value = 'e';
matlabbatch{2}.spm.spatial.coreg.estimate.source(1).sname = 'Realign: Estimate & Reslice: Mean Image';
matlabbatch{2}.spm.spatial.coreg.estimate.source(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.spm.spatial.coreg.estimate.source(1).src_output = substruct('.','rmean');
matlabbatch{2}.spm.spatial.coreg.estimate.other(1) = cfg_dep;
matlabbatch{2}.spm.spatial.coreg.estimate.other(1).tname = 'Other Images';
matlabbatch{2}.spm.spatial.coreg.estimate.other(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{2}.spm.spatial.coreg.estimate.other(1).tgt_spec{1}(1).value = 'image';
matlabbatch{2}.spm.spatial.coreg.estimate.other(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{2}.spm.spatial.coreg.estimate.other(1).tgt_spec{1}(2).value = 'e';
matlabbatch{2}.spm.spatial.coreg.estimate.other(1).sname = 'Realign: Estimate & Reslice: Realigned Images (Sess 1)';
matlabbatch{2}.spm.spatial.coreg.estimate.other(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.spm.spatial.coreg.estimate.other(1).src_output = substruct('.','sess', '()',{1}, '.','cfiles');
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

% Normalization of T1 to MNI and apply parameters to co-registered EPI
matlabbatch{3}.spm.spatial.normalise.estwrite.subj.source = {[t1file,',1']};
matlabbatch{3}.spm.spatial.normalise.estwrite.subj.wtsrc = '';
matlabbatch{3}.spm.spatial.normalise.estwrite.subj.resample = {[t1file,',1']};
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.template = {[mni_template,',1']};
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.weight = '';
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.smosrc = 8;
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.smoref = 0;
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.regtype = 'mni';
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.cutoff = 25;
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.nits = 16;
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.reg = 1;
matlabbatch{3}.spm.spatial.normalise.estwrite.roptions.preserve = 0;
matlabbatch{3}.spm.spatial.normalise.estwrite.roptions.bb = [-90 -126 -72
                                                             90 90 108];
matlabbatch{3}.spm.spatial.normalise.estwrite.roptions.vox = repmat(voxdim,1,3);
matlabbatch{3}.spm.spatial.normalise.estwrite.roptions.interp = 1;
matlabbatch{3}.spm.spatial.normalise.estwrite.roptions.wrap = [0 0 0];
matlabbatch{3}.spm.spatial.normalise.estwrite.roptions.prefix = 'w';

matlabbatch{4}.spm.spatial.normalise.write.subj.matname(1) = cfg_dep;
matlabbatch{4}.spm.spatial.normalise.write.subj.matname(1).tname = 'Parameter File';
matlabbatch{4}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{4}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{4}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{4}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(2).value = 'e';
matlabbatch{4}.spm.spatial.normalise.write.subj.matname(1).sname = 'Normalise: Estimate & Write: Norm Params File (Subj 1)';
matlabbatch{4}.spm.spatial.normalise.write.subj.matname(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{4}.spm.spatial.normalise.write.subj.matname(1).src_output = substruct('()',{1}, '.','params');
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep;
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1).tname = 'Images to Write';
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(1).value = 'image';
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(2).value = 'e';
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1).sname = 'Coreg: Estimate: Coregistered Images';
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1).src_output = substruct('.','cfiles');
matlabbatch{4}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{4}.spm.spatial.normalise.write.roptions.bb = [-90 -126 -72
                                                             90 90 108];
matlabbatch{4}.spm.spatial.normalise.write.roptions.vox = repmat(voxdim,1,3);
matlabbatch{4}.spm.spatial.normalise.write.roptions.interp = 1;
matlabbatch{4}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{4}.spm.spatial.normalise.write.roptions.prefix = 'w';

spm_jobman('run',matlabbatch);