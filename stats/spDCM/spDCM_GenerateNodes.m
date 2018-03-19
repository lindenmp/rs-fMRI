% ------------------------------------------------------------------------------
% This script produces a set of anatomical masks (MNI space) to use with masking
% the main effect of seed GLM to create nodes for spDCM.
% There is a lot of hard coding in this script, thus it is bespoke to my project.
% 
% This script is based on the method of defining DCM nodes outlined in Heim et al.
% 2009, HBM.
% Specifically, it generates nodes on a subject-specific basis constrained by
% anatomy (a priori) and function (second level SPM) using the following steps:
% 1) Define anatomical regions of interest using AAL atlas
% 2) Find peak functional FC between TriStri seed region and AAL mask at the 
% second level (second level SPM needs to be run prior - see FactorialSPM.m)
% 3) Generate 16mm sphere around this peak and mask out any voxels that do not 
% overlap with anatomical prior
% 4) Find subject-specific peak FC within this 16mm mask
% This method ensures the comparability between DCM VOIs across subjects by
% combining anatomical and functional constraints, while still allowing for
% individual variability in FC profiles.
% Linden Parkes, 2018
% ------------------------------------------------------------------------------
clear all; clc

% ------------------------------------------------------------------------------
% Parent dir
% ------------------------------------------------------------------------------
parentdir = '/home/lindenmp/kg98/Linden/';
parentdir_scratch = '/home/lindenmp/kg98_scratch/Linden/';

% ------------------------------------------------------------------------------
% Add paths - edit this section
% ------------------------------------------------------------------------------
% where the prepro scripts are
funcdir = [parentdir,'Scripts/rs-fMRI/func/'];
addpath(funcdir)

fsldir = '/usr/local/fsl/5.0.9/fsl/bin/';
setenv('FSLDIR',fsldir(1:end-4));
setenv('FSLOUTPUTTYPE','NIFTI');
fslstandard = '/usr/local/fsl/5.0.9/fsl/data/standard/';

% Where my collections of ROI atlases are (read only)
ROIDir = [parentdir,'ROIs/'];

% ------------------------------------------------------------------------------
% Set options
% ------------------------------------------------------------------------------
Projects = {'OCDPG_DCM'};
WhichProject = Projects{1}

WhichNoise = 'ICA-AROMA+2P/'
% WhichNoise = 'ICA-AROMA+2P+GSR/'

% TriStri parcel
WhichSeed = 'TriStri'
WhichStri = 2 % dorsal
% WhichStri = 3 % ventral

% ------------------------------------------------------------------------------
% Set project variables
% ------------------------------------------------------------------------------
switch WhichProject
    case 'OCDPG_DCM'
        projdir = [parentdir_scratch,'ResProjects/rfMRI_DCM/OCDPG/'];
        datadir = [projdir,'data/'];
        preprostr = '/func/prepro/';
end

% Where all the nodes/masks etc will be stored
searchVolDir = [projdir,WhichSeed,'_',num2str(WhichStri),'_SearchVolumes/'];
if exist(searchVolDir) == 0
	mkdir(searchVolDir)
elseif exist(searchVolDir) == 7
	rmdir(searchVolDir,'s')
	mkdir(searchVolDir)
end
cd(searchVolDir)

% ------------------------------------------------------------------------------
% Load in TriStri ROIs and create distance penalty masks
% ------------------------------------------------------------------------------
fprintf(1, 'Making penalty masks\n');

% Load TriStri parc
[hdr, parc] = read([ROIDir,'TriStri/TriStri.nii']);

penaltyDir = [searchVolDir,'/TriStri_Penalty/'];
mkdir(penaltyDir)

% this pulls out left only atm
mask = parc == WhichStri;
write(hdr,mask,[penaltyDir,WhichSeed,'_',num2str(WhichStri),'.nii'])

[Vox_coords, ~] = GetROICoords([penaltyDir,WhichSeed,'_',num2str(WhichStri),'.nii']);
CoM = round(mean(Vox_coords));

cd(penaltyDir)
system([fsldir,'fslmaths ',fslstandard,'MNI152_T1_2mm.nii.gz -mul 0 -add 1 -roi ',num2str(CoM(1)),' 1 ',num2str(CoM(2)),' 1 ',num2str(CoM(3)),' 1 0 1 penaltyPoint -odt float']);
system([fsldir,'fslmaths penaltyPoint.nii -kernel sphere 20 -fmean penalty -odt float']);
system([fsldir,'fslmaths penalty.nii -mul 2 -thr `',fsldir,'fslstats penalty.nii -p 100` -bin penalty -odt float']);
[~,penalty] = read([penaltyDir,'penalty.nii']);
penalty = logical(penalty);

% ------------------------------------------------------------------------------
% Make cortical masks
% ------------------------------------------------------------------------------
fprintf(1, 'Making cortical masks\n');

% Load parcellation
[hdr, parc] = read([ROIDir,'AAL/AAL_2mm/raal.nii']);

% Output directory for create masks
maskDir = [searchVolDir,'AAL_Masks/'];
mkdir(maskDir)

% List of masks to create & list of AAL ROIs that constitute the above masks
if WhichStri == 2 
	% DORSAL
	theMasks = {'L_Thal','L_ACC','L_mOFC'};
	ROIs2Use = {77,31,[5,27]};
elseif WhichStri == 3
	% VENTRAL
	theMasks = {'L_Thal','L_ACC','L_mOFC','L_dlPFC'};
	ROIs2Use = {77,31,[5,27],[3,7]};

	% theMasks = {'L_ACC','L_mOFC','L_lOFC','L_dlPFC'};
	% ROIs2Use = {31,[5,27],[9,15],[3,7]};
end

for i = 1:length(theMasks)
	idx = ROIs2Use{i};
	% initialise mask
	mask = zeros(size(parc));

	for j = 1:length(idx)
		mask = mask + double(parc == idx(j));
	end

	write(hdr,mask,[maskDir,theMasks{i},'.nii'])
end

% ------------------------------------------------------------------------------
% Create search volume using second level main effect of seed
% ------------------------------------------------------------------------------
fprintf(1, 'Getting second-level peak t-value from cortical masks\n');

spmDir = [projdir,'SecondLevel/SPM/Factorial/',WhichNoise,WhichSeed,'/'];
if WhichStri == 2
	con1Num = 2; % first level contrast number - e.g., the parcellation in TriStri
elseif WhichStri == 3
	con1Num = 3; % first level contrast number - e.g., the parcellation in TriStri
end

con2Num = 1; % second level contrast number - i.e., whichever contrast denotes the main effect of seed at second level for whole sample
			 % i typicially try to keep this set to the first contrast no matter what I'm doing. 
[hdr,SPM] = read([spmDir,'ROI_',num2str(con1Num),'/spmT_000',num2str(con2Num),'.nii']);

% store search volumes in a cell
searchVol = cell(1,length(theMasks));

for i = 1:length(theMasks)
	% load mask
	[~,mask] = read([maskDir,theMasks{i},'.nii']);
	mask = logical(mask);

	% apply distance penalty
	mask(penalty) = 0;

	% mask SPM
	SPMm = SPM;
	SPMm(~mask) = NaN;

	% find max stat
	[M,I] = max(SPMm(:));

	% got voxel coords of max stat
	[I1,I2,I3] = ind2sub(size(SPMm),I);

	% create search volume of 16mm radius around vox coords
	system([fsldir,'fslmaths ',fslstandard,'MNI152_T1_2mm.nii.gz -mul 0 -add 1 -roi ',num2str(I1),' 1 ',num2str(I2),' 1 ',num2str(I3),' 1 0 1 searchVolPoint -odt float']);
	system([fsldir,'fslmaths searchVolPoint.nii -kernel sphere 16 -fmean searchVol.nii -odt float']);
	system([fsldir,'fslmaths searchVol.nii -mul 2 -thr `',fsldir,'fslstats searchVol.nii -p 100` -bin searchVol -odt float']);
	[hdr_searchVol,searchVol{i}] = read('searchVol.nii');

	% mask out any voxels in search volume that are not within original anatomical mask and rename
	searchVol{i}(~mask) = 0;

	% writing out search volume
	write(hdr_searchVol,searchVol{i},[searchVolDir,'search_',theMasks{i},'.nii'])

	delete('searchVolPoint.nii')
	delete('searchVol.nii')
end

% ------------------------------------------------------------------------------
% Generate subject-specific nodes using first level main effect of seed
% ------------------------------------------------------------------------------
load([spmDir,'metadata.mat'])
numSubs = size(metadata,1);

% Focusing on left hemisphere right now
firstdir = [preprostr,'/',WhichNoise,'FirstLevel_L_',WhichSeed,'/'];

storeNode = cell(length(theMasks),1);

for i = 1:numSubs
	fprintf(1, 'Generating DCM nodes for %s\n', metadata.ParticipantID{i});
	nodeDir = [datadir,metadata.ParticipantID{i},preprostr,WhichNoise,WhichSeed,'_',num2str(WhichStri),'_nodes_wThal/'];
	if exist(nodeDir) == 0
		mkdir(nodeDir)
	elseif exist(nodeDir) == 7
		rmdir(nodeDir,'s')
		mkdir(nodeDir)
	end
	cd(nodeDir)

	% first, just duplicate the TriStri ROI that is needed into each subjects directory
	copyfile([penaltyDir,WhichSeed,'_',num2str(WhichStri),'.nii'],nodeDir)

	% load first level SPM
	[hdr,SPM] = read([datadir,'/',metadata.ParticipantID{i},firstdir,'/con_000',num2str(con1Num),'.nii,1']);

	for j = 1:length(theMasks)
		% load search volume
		% [~,searchVol] = read([searchVolDir,'search_',theMasks{j},'.nii']);

		% mask SPM
		SPMm = SPM;
		SPMm(~searchVol{j}) = NaN;

		% find max stat
		[M,I] = max(SPMm(:));

		% got voxel coords of max stat
		[I1,I2,I3] = ind2sub(size(SPMm),I);

		% create node of 3mm radius around vox coords
		system([fsldir,'fslmaths ',fslstandard,'MNI152_T1_2mm.nii.gz -mul 0 -add 1 -roi ',num2str(I1),' 1 ',num2str(I2),' 1 ',num2str(I3),' 1 0 1 -bin nodePoint -odt float']);
		system([fsldir,'fslmaths nodePoint.nii -kernel sphere 3 -fmean node.nii -odt float']);
		system([fsldir,'fslmaths node.nii -mul 2 -thr `',fsldir,'fslstats node.nii -p 100` -bin node -odt float']);
		[hdr_node,node] = read('node.nii');

		% writing out search volume
		write(hdr_node,node,[nodeDir,theMasks{j},'.nii'])

		% Store node point in matrix for plotting
		[~,nodePoint] = read('nodePoint.nii');
		storeNode{j}(:,:,:,i) = nodePoint;

		delete('nodePoint.nii')
		delete('node.nii')
	end

	% duplicate without Thal
	nodeDir2 = [datadir,metadata.ParticipantID{i},preprostr,WhichNoise,WhichSeed,'_',num2str(WhichStri),'_nodes/'];
	copyfile(nodeDir,nodeDir2)
	cd(nodeDir2)
	delete('L_Thal.nii')
end

% ------------------------------------------------------------------------------
% Produce probabilistic node maps for each search vol
% ------------------------------------------------------------------------------
cd(searchVolDir)

for i = 1:length(theMasks)
	nodeProb{i} = sum(storeNode{i},4)./numSubs;
	write(hdr_searchVol,nodeProb{i},[searchVolDir,'prob_',theMasks{i},'.nii'])
end
