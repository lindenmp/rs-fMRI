% ------------------------------------------------------------------------------
% 
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

% ------------------------------------------------------------------------------
% Set project variables
% ------------------------------------------------------------------------------
switch WhichProject
    case 'OCDPG_DCM'
        projdir = [parentdir_scratch,'ResProjects/rfMRI_DCM/OCDPG/'];
        datadir = [projdir,'data/'];
        preprostr = '/func/prepro/';

        EPI = 'epi_prepro.nii';
end

% Where all the nodes/masks etc will be stored
corVolDir = [projdir,'CorVolumes/'];
if exist(corVolDir) == 0
	mkdir(corVolDir)
elseif exist(corVolDir) == 7
	rmdir(corVolDir,'s')
	mkdir(corVolDir)
end
cd(corVolDir)

% ------------------------------------------------------------------------------
% Make cortical masks
% ------------------------------------------------------------------------------
fprintf(1, 'Making cortical masks\n');

% Load parcellation
[hdr, parc] = read([ROIDir,'AAL/AAL_2mm/raal.nii']);

% Output directory for create masks
maskDir = [corVolDir,'AAL_Masks/'];
mkdir(maskDir)

% List of masks to create & list of AAL ROIs that constitute the above masks
theMasks = {'L_Thal','R_Thal';'L_ACC','R_ACC';'L_mOFC','R_mOFC';'L_dlPFC','R_dlPFC'}';
ROIs2Use = {77,78;31,32;[5,27],[6,28];[3,7],[4,8]}';

for i = 1:numel(theMasks)
	idx = ROIs2Use{i};
	% initialise mask
	mask = zeros(size(parc));

	for j = 1:length(idx)
		mask = mask + double(parc == idx(j));
	end

	write(hdr,mask,[maskDir,theMasks{i},'.nii'])
end

% also add in the TriStri parcellation
[hdr, parc] = read([ROIDir,'TriStri/TriStri.nii']);
theMasks_TriStri = {'L_Caud','R_Caud';'L_Dors','R_Dors';'L_Vent','R_Vent'}';
ROIs2Use_TriStri = {1,4;2,5;3,6}';

for i = 1:numel(theMasks_TriStri)
	idx = ROIs2Use_TriStri{i};
	% initialise mask
	mask = zeros(size(parc));

	for j = 1:length(idx)
		mask = mask + double(parc == idx(j));
	end

	write(hdr,mask,[maskDir,theMasks_TriStri{i},'.nii'])
end

theMasks = [theMasks,theMasks_TriStri];
ROIs2Use = [ROIs2Use,ROIs2Use_TriStri];

% compile into parcellation files
numMaskPairs = size(theMasks,2);

for i = 1:2
	mask = zeros(size(parc));
	for j = 1:numMaskPairs
		[~,temp] = read([maskDir,theMasks{i,j}]);
		temp = double(temp);
		temp(temp == 1) = j;
		mask = mask + temp;
	end

	% remove any overlap
	mask(mask > numMaskPairs) = 0;

	if i == 1
		write(hdr,mask,[maskDir,'lh.nii'])
	elseif i == 2
		write(hdr,mask,[maskDir,'rh.nii'])
	end
end

% ------------------------------------------------------------------------------
% Get time series and run correlation
% ------------------------------------------------------------------------------
spmDir = [projdir,'SecondLevel/SPM/Factorial/',WhichNoise,WhichSeed,'/'];
load([spmDir,'metadata.mat'])
numSubs = size(metadata,1);

R = zeros(numSubs,numMaskPairs);

for i = 1:numSubs
	fprintf(1, 'Correlating hemis for %s\n', metadata.ParticipantID{i});
	workDir = [datadir,metadata.ParticipantID{i},preprostr,WhichNoise,'/'];

	% left hemi
	system([fsldir,'fslmeants -i ',[workDir,EPI],' --label=',maskDir,'lh.nii -o lh.txt']);
	lh = dlmread('lh.txt');
	delete('lh.txt')

	% right hemi
	system([fsldir,'fslmeants -i ',[workDir,EPI],' --label=',maskDir,'rh.nii -o rh.txt']);
	rh = dlmread('rh.txt');
	delete('rh.txt')

	r = corr(lh,rh,'type','Pearson');
	R(i,:) = diag(r);
end

% ------------------------------------------------------------------------------
% Plot
% ------------------------------------------------------------------------------
figure('color','w'); box('on'); hold on;
for i = 1:numMaskPairs
	histogram(R(:,i),10)
end
maskNames = {theMasks{1,:}};
maskNames = strrep(maskNames,'L_','');
legend(maskNames)

save('R.mat','R','theMasks','ROIs2Use','metadata')
