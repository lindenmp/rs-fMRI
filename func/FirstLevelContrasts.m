function [] = FirstLevelContrasts(WhichSeed,matfile)
%
%
% ------
% INPUTS
% ------
%
% matfile 	- path and file name of first level estimation output .mat file (typically SPM.mat)
%
% -------
% OUTPUTS
% -------
%
% SPM.mat file containing contrats for ROIs for 1st level analysis. located in spmdir.
%
% =========================================================================

spm('defaults','fmri');
spm_jobman('initcfg')
spm_get_defaults('mask.thresh',-Inf)

% run contrasts
switch WhichSeed
	case 'TriStri'
		matlabbatch{1}.spm.stats.con.spmmat = {matfile};
		matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Sensorimotor';
		matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = [1 0 0];
		matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
		matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Associative';
		matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = [0 1 0];
		matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
		matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'Ventral';
		matlabbatch{1}.spm.stats.con.consess{3}.tcon.convec = [0 0 1];
		matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
		matlabbatch{1}.spm.stats.con.delete = 0;
	case 'DiMartino'
		matlabbatch{1}.spm.stats.con.spmmat = {matfile};
		matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'VSi';
		matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = [1 0 0 0 0 0];
		matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
		matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'VSs';
		matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = [0 1 0 0 0 0];
		matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
		matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'DC';
		matlabbatch{1}.spm.stats.con.consess{3}.tcon.convec = [0 0 1 0 0 0];
		matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
		matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'DCP';
		matlabbatch{1}.spm.stats.con.consess{4}.tcon.convec = [0 0 0 1 0 0];
		matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
		matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'DRP';
		matlabbatch{1}.spm.stats.con.consess{5}.tcon.convec = [0 0 0 0 1 0];
		matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
		matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'VRP';
		matlabbatch{1}.spm.stats.con.consess{6}.tcon.convec = [0 0 0 0 0 1];
		matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
		matlabbatch{1}.spm.stats.con.delete = 0;
end

spm_jobman('run',matlabbatch);

clear matlabbatch
