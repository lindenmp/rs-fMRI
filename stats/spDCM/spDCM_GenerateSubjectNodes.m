%% spDCM_GenerateSubjectNodes: 
function [] = spDCM_GenerateSubjectNodes(searchVolDir,con1Name,seedName,nodeDir)
	% ------------------------------------------------------------------------------
	% searchVolDir			- string directing to search volumes generated using 
	% 						spDCM_GenerateSearchVolumes.m
	% con1Name				- string of the directory and name contrast image to use.
	% 						e.g, /path/to/first/level/SPM/con_0001.nii
	% seedName				- string of the seed image used to generate con image.
	% 						This will be stored in searchVolDir
	% nodeDir 				- string for output directory, this is where the nodes
	% 						will be saved.
	% ------------------------------------------------------------------------------
	% Generate subject-specific nodes using first level main effect of seed
	% ------------------------------------------------------------------------------
	fprintf(1, 'Generating DCM nodes \n');
	rng('default')

    % ------------------------------------------------------------------------------
    % Parent dir
    % ------------------------------------------------------------------------------
    parentdir = '/home/lindenmp/kg98/Linden/';

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

	% ------------------------------------------------------------------------------
	% Find search volumes
	% ------------------------------------------------------------------------------
	searchVols = dir(searchVolDir);
	searchVols = {searchVols.name};
	idx = find(contains(searchVols,'.nii'));
	searchVols = searchVols(idx);
	
	theMasks = searchVols;
	% theMasks = strrep(searchVols,'search_','');

	storeNode = cell(length(theMasks),1);

	if exist(nodeDir) == 0
		mkdir(nodeDir)
	elseif exist(nodeDir) == 7
		rmdir(nodeDir,'s')
		mkdir(nodeDir)
	end
	cd(nodeDir)

	% first, just duplicate the TriStri ROI that is needed into each subjects directory
	penaltyDir = [searchVolDir,'/Penalty/'];
	copyfile([penaltyDir,seedName],nodeDir)

	% load first level SPM
	[hdr,SPM] = read([con1Name,',1']);

	for j = 1:length(theMasks)
		% load search volume
		[~,searchVol] = read([searchVolDir,searchVols{j}]);

		% mask SPM
		SPMm = SPM;
		SPMm(~searchVol) = NaN;

		% find max stat
		[M,I] = max(SPMm(:));

		% got voxel coords of max stat
		[I1,I2,I3] = ind2sub(size(SPMm),I);

		% create node of 3mm radius around vox coords
		system([fsldir,'fslmaths ',fslstandard,'MNI152_T1_2mm.nii.gz -mul 0 -add 1 -roi ',num2str(I1),' 1 ',num2str(I2),' 1 ',num2str(I3),' 1 0 1 -bin nodePoint -odt float']);
		system([fsldir,'fslmaths nodePoint.nii -kernel sphere 3 -fmean node.nii -odt float']);
		system([fsldir,'fslmaths node.nii -mul 2 -thr `',fsldir,'fslstats node.nii -p 100` -bin node -odt float']);
		[hdr_node,node] = read('node.nii');

		% writing out node
		write(hdr_node,node,[nodeDir,theMasks{j}])

		% Store node point in matrix for plotting
		[~,nodePoint] = read('nodePoint.nii');
		storeNode{j} = nodePoint;

		delete('nodePoint.nii')
		delete('node.nii')
	end
	fprintf(1, 'Done \n');
end
