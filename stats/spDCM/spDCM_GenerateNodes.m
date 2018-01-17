% ------------------------------------------------------------------------------
% This script produces a set of anatomical masks (MNI space) to use with masking
% the main effect of seed GLM to create nodes for spDCM.
% There is a lot of hard coding in this script, thus it is bespoke to my project.
% Linden Parkes, 2017
% ------------------------------------------------------------------------------
clear all; clc
% fsldir = '/usr/local/fsl/bin/';
fsldir = '/usr/local/fsl/5.0.9/fsl/bin/';
setenv('FSLDIR',fsldir(1:end-4));
setenv('FSLOUTPUTTYPE','NIFTI');
fslstandard = '/usr/local/fsl/5.0.9/fsl/data/standard/';
% fslstandard = '/usr/local/fsl/data/standard/';

projdir = '/scratch/kg98/Linden/ResProjects/rfMRI_DCM/OCDPG/';
% projdir = '~/Dropbox/Work/ResProjects/rfMRI_DCM/OCDPG/';

ROIDir = '/projects/kg98/Linden/ROIs/';
% ROIDir = '~/Dropbox/Work/ROIs/';

% ------------------------------------------------------------------------------
% Load in TriStri ROIs and create distance penalty masks
% ------------------------------------------------------------------------------
fprintf(1, 'Making penalty masks\n');

% Load TriStri parc
[hdr, Parc] = read([ROIDir,'TriStri/TriStri.nii']);

outDir = [projdir,'ROIs/TriStri/'];
if exist(outDir) == 0
	mkdir(outDir)
elseif exist(outDir) == 7
	rmdir(outDir,'s')
	mkdir(outDir)
end
cd(outDir)

i = 3;
Mask = Parc == i;
write(hdr,Mask,[num2str(i),'.nii'])

[Vox_coords, ~] = GetROICoords([num2str(i),'.nii']);
CoM = round(mean(Vox_coords));

system([fsldir,'fslmaths ',fslstandard,'MNI152_T1_2mm.nii.gz -mul 0 -add 1 -roi ',num2str(CoM(1)),' 1 ',num2str(CoM(2)),' 1 ',num2str(CoM(3)),' 1 0 1 PenaltyPoint -odt float']);
system([fsldir,'fslmaths PenaltyPoint.nii -kernel box 40 -fmean Penalty -odt float']);
system([fsldir,'fslmaths Penalty.nii -bin Penalty']);

[~,Penalty] = read('Penalty.nii');

% ------------------------------------------------------------------------------
% Make cortical masks
% ------------------------------------------------------------------------------
fprintf(1, 'Making cortical masks\n');

Atlases = {'AAL'};
WhichAtlas = Atlases{1};

switch WhichAtlas
	case 'AAL'
		% Load parcellation
		[hdr, Parc] = read([ROIDir,'AAL/AAL_2mm/raal.nii']);

		% Output directory for create masks
		maskDir = [projdir,'ROIs/AAL_Masks/'];
		if exist(maskDir) == 0
			mkdir(maskDir)
		elseif exist(maskDir) == 7
			rmdir(maskDir,'s')
			mkdir(maskDir)
		end
		cd(maskDir)

		% List of masks to create
		theMasks = {'L_ACC','L_lOFC','L_mOFC','L_Thal'};
		% List of AAL ROIs that constitute the above masks
		ROIs2Use = {31,[9,15],[5,27],77};

		for i = 1:length(theMasks)
			idx = ROIs2Use{i};
			% initialise mask
			Mask = zeros(size(Parc));

			for j = 1:length(idx)
				Mask = Mask + double(Parc == idx(j));
			end

			write(hdr,Mask,[theMasks{i},'.nii'])
		end
end

% ------------------------------------------------------------------------------
% Get peaks for cortical masks
% ------------------------------------------------------------------------------
fprintf(1, 'Getting peak t-value from cortical masks\n');

nodeDir = [projdir,'nodes/'];
if exist(nodeDir) == 0
	mkdir(nodeDir)
elseif exist(nodeDir) == 7
	rmdir(nodeDir,'s')
	mkdir(nodeDir)
end
cd(nodeDir)

spmDir = '/scratch/kg98/Linden/ResProjects/rfMRI_DCM/OCDPG/SecondLevelSPM/Factorial/ICA-AROMA+2P/HC+OCD/TriStri/ROI_3/';
% spmDir = '~/Dropbox/Work/ResProjects/rfMRI_DCM/OCDPG/SecondLevelSPM/Factorial/ICA-AROMA+2P/HC+OCD/TriStri/ROI_3/';
conNum = 3;

[hdr,SPM] = read([spmDir,'spmT_000',num2str(conNum),'.nii']);

for i = 1:length(theMasks)
	% load mask
	[~,Mask] = read([maskDir,theMasks{i},'.nii']);

	% subtract penalty
	Mask = double(Mask) - Penalty;

	% mask SPM
	SPMm = SPM;
	SPMm(Mask ~= 1) = NaN;

	% find max stat
	[M,I] = max(SPMm(:));

	% got voxel coords of max stat
	[I1,I2,I3] = ind2sub(size(SPMm),I);

	% create spherical node around vox coords
	system([fsldir,'fslmaths ',fslstandard,'MNI152_T1_2mm.nii.gz -mul 0 -add 1 -roi ',num2str(I1),' 1 ',num2str(I2),' 1 ',num2str(I3),' 1 0 1 NodePoint -odt float']);
	system([fsldir,'fslmaths NodePoint.nii -kernel sphere 4 -fmean Node -odt float']);
	delete('NodePoint.nii')
	system([fsldir,'fslmaths Node.nii -bin ',theMasks{i}]);
	delete('Node.nii')
end

