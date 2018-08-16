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
clear all; close all; clc
rng('default')

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
Projects = {'OCDPG_DCM','GenCog'};
WhichProject = Projects{1}

% TriStri parcel
WhichSeed = 'TriStri'
Stris = [2,3]; % 2 = dorsal, 3 = ventral
offset = 3;

% DiMartino
% WhichSeed = 'DiMartino'
% Stris = [3,1]; % 3 = DC, 1 = VSi
% offset = 6;

projection = 0 % ipsilateral projections
% projection = 1 % contralateral projections
% projection = 2 % bilateral projections

penaltySize = 20 % for ocdpg
% penaltySize = 30 % for gencog

searchSize = 16
% searchSize = 14
% searchSize = 12
% searchSize = 10

searchSizeThal = 12

extraStr = '';

% ------------------------------------------------------------------------------
% Set project variables
% ------------------------------------------------------------------------------
switch WhichProject
    case 'OCDPG_DCM'
        projdir = [parentdir_scratch,'ResProjects/rfMRI_DCM/OCDPG/'];
        datadir = [projdir,'data/'];
        preprostr = '/func/prepro/';
		
		WhichNoise = 'ICA-AROMA+2P/'
		% WhichNoise = 'ICA-AROMA+2P+GSR/'
    case 'GenCog'
        projdir = [parentdir_scratch,'ResProjects/rfMRI_DCM/GenCog/'];
        datadir = [projdir,'data/'];
        preprostr = '/func/prepro/';

		WhichNoise = 'ICA-FIX/'
		% WhichNoise = 'ICA-FIX+GSR/'
end

% ------------------------------------------------------------------------------
% Generate search vols
% ------------------------------------------------------------------------------
for WhichStri = Stris
	% create string label
	if WhichStri == 2 & ismember('TriStri',WhichSeed,'rows')
		StriLabel = 'dorsal'; % first level contrast number - e.g., the parcellation in TriStri
	elseif WhichStri == 3 & ismember('TriStri',WhichSeed,'rows')
		StriLabel = 'ventral'; % first level contrast number - e.g., the parcellation in TriStri
	elseif WhichStri == 3 & ismember('DiMartino',WhichSeed,'rows')
		StriLabel = 'dorsal'; % first level contrast number
	elseif WhichStri == 1 & ismember('DiMartino',WhichSeed,'rows')
		StriLabel = 'ventral'; % first level contrast number
	end	
	
	WhichStri
	StriLabel
	
	for h = 1:2
		if h == 1
			hemi = 'L';
			StriIdx = WhichStri;
		elseif h == 2;
			hemi = 'R';
			StriIdx = WhichStri + offset;
		end

		fprintf(1, 'Making search volumes for %s hemisphere \n',hemi);

		% Where all the nodes/masks etc will be stored
		if projection == 0
			searchVolDir = [projdir,'SecondLevel/spDCM/SearchVolumes',extraStr,'/',WhichNoise,WhichSeed,'_',num2str(WhichStri),'_',hemi,'/'];
		elseif projection == 1
			searchVolDir = [projdir,'SecondLevel/spDCM/SearchVolumes',extraStr,'/',WhichNoise,WhichSeed,'_',num2str(WhichStri),'_',hemi,'c/'];
		elseif projection == 2
			searchVolDir = [projdir,'SecondLevel/spDCM/SearchVolumes',extraStr,'/',WhichNoise,WhichSeed,'_',num2str(WhichStri),'_',hemi,'b/'];
		end

		if exist(searchVolDir) == 0
			mkdir(searchVolDir)
		elseif exist(searchVolDir) == 7
			rmdir(searchVolDir,'s')
			mkdir(searchVolDir)
		end
		cd(searchVolDir)

		% ------------------------------------------------------------------------------
		% Write params out for future reference
		% ------------------------------------------------------------------------------
		csvwrite('projection.txt',projection)
		csvwrite('penaltySize.txt',penaltySize)
		csvwrite('searchSize.txt',searchSize)
		csvwrite('searchSizeThal.txt',searchSizeThal)

		% ------------------------------------------------------------------------------
		% Load in ROIs and create distance penalty masks
		% ------------------------------------------------------------------------------
		fprintf(1, '\tMaking penalty mask\n');

		% Load striatum parc
		if ismember('TriStri',WhichSeed,'rows')
			[hdr, parc] = read([ROIDir,'TriStri/TriStri.nii']);
		elseif ismember('DiMartino',WhichSeed,'rows')
			[hdr, parc] = read([ROIDir,'DiMartino/SphereParc02.nii']);
		end

		penaltyDir = [searchVolDir,'/Penalty/'];
		mkdir(penaltyDir)

		mask = parc == StriIdx;
		write(hdr,mask,[penaltyDir,WhichSeed,'_',num2str(WhichStri),'_',hemi,'.nii'])
		[Vox_coords, ~] = GetROICoords([penaltyDir,WhichSeed,'_',num2str(WhichStri),'_',hemi,'.nii']);

		cd(penaltyDir)
		CoM = round(mean(Vox_coords));
		system([fsldir,'fslmaths ',fslstandard,'MNI152_T1_2mm.nii.gz -mul 0 -add 1 -roi ',num2str(CoM(1)),' 1 ',num2str(CoM(2)),' 1 ',num2str(CoM(3)),' 1 0 1 penaltyPoint -odt float']);
		system([fsldir,'fslmaths penaltyPoint.nii -kernel sphere ',num2str(penaltySize),' -fmean penalty -odt float']);
		system([fsldir,'fslmaths penalty.nii -mul 2 -thr `',fsldir,'fslstats penalty.nii -p 100` -bin penalty -odt float']);
		[~,penalty] = read([penaltyDir,'penalty.nii']);
		penalty = logical(penalty);

		% ------------------------------------------------------------------------------
		% Make cortical masks
		% ------------------------------------------------------------------------------
		fprintf(1, '\tMaking cortical masks\n');

		% Load parcellation
		[hdr, parc] = read([ROIDir,'AAL/AAL_2mm/raal.nii']);

		% Output directory for create masks
		maskDir = [searchVolDir,'AAL_Masks/'];
		mkdir(maskDir)

		% List of masks to create & list of AAL ROIs that constitute the above masks
		if ismember('dorsal',StriLabel,'rows') & h == 1 & projection == 0 | ismember('dorsal',StriLabel,'rows') & h == 2 & projection == 1
			% DORSAL LEFT
			theMasks = {'Thal_L','ACC_L','mOFC_L'};
			ROIs2Use = {77,31,[5,27]};
		elseif ismember('dorsal',StriLabel,'rows') & h == 2 & projection == 0 | ismember('dorsal',StriLabel,'rows') & h == 1 & projection == 1
			% DORSAL RIGHT
			theMasks = {'Thal_R','ACC_R','mOFC_R'};
			ROIs2Use = {78,32,[6,28]};
		elseif ismember('ventral',StriLabel,'rows') & h == 1 & projection == 0 | ismember('ventral',StriLabel,'rows') & h == 2 & projection == 1
			% VENTRAL LEFT
			theMasks = {'Thal_L','ACC_L','mOFC_L','dlPFC_L'};
			ROIs2Use = {77,31,[5,27],[3,7]};
		elseif ismember('ventral',StriLabel,'rows') & h == 2 & projection == 0 | ismember('ventral',StriLabel,'rows') & h == 1 & projection == 1
			% VENTRAL RIGHT
			theMasks = {'Thal_R','ACC_R','mOFC_R','dlPFC_R'};
			ROIs2Use = {78,32,[6,28],[4,8]};
		elseif ismember('dorsal',StriLabel,'rows') & projection == 2
			% DORSAL BILAT
			theMasks = {'Thal_L','ACC_L','mOFC_L','Thal_R','ACC_R','mOFC_R'};
			ROIs2Use = {77,31,[5,27],78,32,[6,28]};
		elseif ismember('ventral',StriLabel,'rows') & projection == 2
			% VENTRAL BILAT
			theMasks = {'Thal_L','ACC_L','mOFC_L','dlPFC_L','Thal_R','ACC_R','mOFC_R','dlPFC_R'};
			ROIs2Use = {77,31,[5,27],[3,7],78,32,[6,28],[4,8]};
		end

		for i = 1:length(theMasks)
			idx = ROIs2Use{i};
			% initialise mask
			mask = zeros(size(parc));

			for j = 1:length(idx)
				mask = mask + double(parc == idx(j));
			end

			write(hdr,mask,[maskDir,theMasks{i},'.nii'])

			% perform single erosing of mask
			% system([fsldir,'fslmaths ',maskDir,theMasks{i},'.nii -eroF -bin ',maskDir,theMasks{i},'.nii']);
		end

		% ------------------------------------------------------------------------------
		% Create search volume using second level main effect of seed
		% ------------------------------------------------------------------------------
		fprintf(1, '\tGetting second-level peak t-value from cortical masks\n');

		spmDir = [projdir,'SecondLevel/SPM/Factorial/',WhichNoise,WhichSeed,extraStr,'/'];

		if h == 1
			con2Num = 2; % second level contrast number - i.e., whichever contrast denotes the main effect of seed at second level for whole sample
		elseif h == 2
			con2Num = 3; % second level contrast number - i.e., whichever contrast denotes the main effect of seed at second level for whole sample
		end
		% con2Num = 1; % second level contrast number - i.e., whichever contrast denotes the main effect of seed at second level for whole sample

		[hdr,SPM] = read([spmDir,'ROI_',num2str(WhichStri),'/spmT_000',num2str(con2Num),'.nii']);

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

			% create search volume around vox coords
			system([fsldir,'fslmaths ',fslstandard,'MNI152_T1_2mm.nii.gz -mul 0 -add 1 -roi ',num2str(I1),' 1 ',num2str(I2),' 1 ',num2str(I3),' 1 0 1 searchVolPoint -odt float']);
			if contains(theMasks{i},'Thal')
				system([fsldir,'fslmaths searchVolPoint.nii -kernel sphere ',num2str(searchSizeThal),' -fmean searchVol.nii -odt float']);
			elseif ~contains(theMasks{i},'Thal')
				system([fsldir,'fslmaths searchVolPoint.nii -kernel sphere ',num2str(searchSize),' -fmean searchVol.nii -odt float']);
			end

			system([fsldir,'fslmaths searchVol.nii -mul 2 -thr `',fsldir,'fslstats searchVol.nii -p 100` -bin searchVol -odt float']);
			[hdr_searchVol,searchVol] = read('searchVol.nii');

			% mask out any voxels in search volume that are not within original anatomical mask and rename
			searchVol(~mask) = 0;

			% writing out search volume
			write(hdr_searchVol,searchVol,[searchVolDir,WhichSeed,'_',num2str(WhichStri),'_',hemi,'-',theMasks{i},'.nii'])
			% write out centroids for future refernece
			csvwrite([searchVolDir,WhichSeed,'_',num2str(WhichStri),'_',hemi,'-',theMasks{i},'.txt'],[I1,I2,I3])

			delete('searchVolPoint.nii')
			delete('searchVol.nii')
		end
	end
end

% ------------------------------------------------------------------------------
% Check for overlap in duplicates
% ------------------------------------------------------------------------------
% dorsal
theMasks = {'ACC_L.nii','ACC_R.nii';'mOFC_L.nii','mOFC_R.nii';'Thal_L.nii','Thal_R.nii'}'; % '

% compile into parcellation files
numMaskPairs = size(theMasks,2);

for h = 1:2
	if h == 1
		hemi = 'L';
	elseif h == 2;
		hemi = 'R';
	end
	for i = 1:numMaskPairs
		% DORSAL
		searchVolDir = [projdir,'SecondLevel/spDCM/SearchVolumes',extraStr,'/',WhichNoise,WhichSeed,'_',num2str(Stris(1)),'_',hemi,'/'];
		[hdr_dors,dors] = read([searchVolDir,WhichSeed,'_',num2str(Stris(1)),'_',hemi,'-',theMasks{h,i}]);
		
		% VENTRAL
		searchVolDir = [projdir,'SecondLevel/spDCM/SearchVolumes',extraStr,'/',WhichNoise,WhichSeed,'_',num2str(Stris(2)),'_',hemi,'/'];
		[hdr_vent,vent] = read([searchVolDir,WhichSeed,'_',num2str(Stris(2)),'_',hemi,'-',theMasks{h,i}]);

		% multiply
		overlap = dors.*vent;

		if numel(unique(overlap)) > 1
			fprintf(1, 'There are overlapping voxels in %s \n', theMasks{h,i});

			% subtract out the overlap
			dors = dors - overlap;
			vent = vent - overlap;

			% write out new search vols
			% DORSAL
			searchVolDir = [projdir,'SecondLevel/spDCM/SearchVolumes',extraStr,'/',WhichNoise,WhichSeed,'_',num2str(Stris(1)),'_',hemi,'/'];
			write(hdr_dors,dors,[searchVolDir,WhichSeed,'_',num2str(Stris(1)),'_',hemi,'-',theMasks{h,i}])
			
			% VENTRAL
			searchVolDir = [projdir,'SecondLevel/spDCM/SearchVolumes',extraStr,'/',WhichNoise,WhichSeed,'_',num2str(Stris(2)),'_',hemi,'/'];
			write(hdr_vent,vent,[searchVolDir,WhichSeed,'_',num2str(Stris(2)),'_',hemi,'-',theMasks{h,i}])
		elseif numel(unique(overlap)) == 1
			fprintf(1, 'There are NO overlapping voxels in %s \n', theMasks{h,i});
		end
	end
end
