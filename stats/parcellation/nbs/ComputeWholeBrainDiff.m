function [] = ComputeWholeBrainDiff(WhichProject,WhichSplit,WhichParc,WhichNoise,P,outDir,runCensor)
	% ------------------------------------------------------------------------------
	% ComputeWholeBrainDiff.m
	% ------------------------------------------------------------------------------
	% This functon will run nbs contrasts the OCDPG or UCLA datasets.
	% contrasts can be performed either 1) case-control or 2) high v low motion
	% 
	% ------------------------------------------------------------------------------
	% Linden Parkes, Brain & Mental Health Laboratory, 2016
	% ------------------------------------------------------------------------------
	if nargin < 5
		P = 0.05;
	end

	if nargin < 6
		outDir = 'NBS';
	end

	if nargin < 7
		runCensor = 'false';
	end

	fprintf(1, 'Dataset: %s\n Split: %s\n Parcellation: %s\n Noise correction: %s\n', WhichProject,WhichSplit,WhichParc,WhichNoise);

	switch WhichParc
		case 'Gordon'
			Parc = 1;
		case 'Power'
			Parc = 2;
	end

	% ------------------------------------------------------------------------------
	% Add paths - edit this section
	% ------------------------------------------------------------------------------
	funcdir = '/gpfs/M2Home/projects/Monash076/Linden/scripts/rs-fMRI/func/';
	addpath(funcdir)
	funcdir1 = '/gpfs/M2Home/projects/Monash076/Linden/scripts/Func/';
	addpath(funcdir1)

	addpath(genpath('/gpfs/M2Home/projects/Monash076/Linden/scripts/Software/NBS1.2'))

	% ------------------------------------------------------------------------------
	% Select project
	% ------------------------------------------------------------------------------
	switch WhichProject
		case 'OCDPG'
			sublist = '/gpfs/M2Home/projects/Monash076/Linden/Sublists/OCDPGe.csv';
			projdir = '/gpfs/M2Home/projects/Monash076/Linden/OCDPG/';
			datadir = [projdir,'data/'];
			
			preprostr = '/rfMRI/prepro/';
			TR = 2.5;
		case 'UCLA'
			sublist = '/gpfs/M2Home/projects/Monash076/Linden/Sublists/UCLA.csv';
			projdir = '/gpfs/M2Home/projects/Monash076/Linden/UCLA/';
			datadir = [projdir,'data/'];
	        
	        preprostr = '/rfMRI/prepro/';
			TR = 2;
	end

	% ------------------------------------------------------------------------------
	% Setup output directory
	% ------------------------------------------------------------------------------
	cd(projdir)

	if exist(outDir) == 0
		fprintf(1,'Initialising outDir\n')
		mkdir(outDir)
	end

	cd([projdir,outDir])

	% ------------------------------------------------------------------------------
	% Subject list
	% ------------------------------------------------------------------------------
	fileID = fopen(sublist);
	switch WhichProject
		case 'OCDPG'
			metadata = textscan(fileID, '%s %u %u %s %s','HeaderLines',1, 'delimiter',',');
			metadata{2} = double(metadata{2}); metadata{3} = double(metadata{3});

			Cov = [];
			% age
			Cov(:,1) = metadata{3};
			% mean center age
			Cov(:,1) = Cov(:,1) - mean(Cov(:,1));
			% gender
			Cov(strmatch('F',metadata{4}),2) = 1; Cov(strmatch('M',metadata{4}),2) = 0;
			zeroPad = ',0,0';
		case 'UCLA'
			metadata = textscan(fileID, '%s %u %u %s %s %u','HeaderLines',1, 'delimiter',',');
			metadata{2} = double(metadata{2}); metadata{3} = double(metadata{3});
			
			Cov = [];
			% age
			Cov(:,1) = metadata{3};
			% mean center age
			Cov(:,1) = Cov(:,1) - mean(Cov(:,1));
			% gender
			Cov(strmatch('F',metadata{4}),2) = 1; Cov(strmatch('M',metadata{4}),2) = 0;
			% scanner
			Cov(:,3) = metadata{6};
			zeroPad = ',0,0,0';
	end

	ParticipantIDs = metadata{1};
	Group_Diagnostic = metadata{2};

	% ------------------------------------------------------------------------------
	% Jenkinson's mean FD
	% ------------------------------------------------------------------------------
	fprintf(1, 'Loading Jenkinson''s mean FD metric\n');
	[exclude,~,fdJenk,fdJenk_m] = GetExcludeForSample(datadir,ParticipantIDs,preprostr);
	fprintf(1, 'done\n');

	% compute number of volumes using the length of fdJenk
	% note, this is assumed to be same for all subjects!
	numVols = length(fdJenk{1});

	% ------------------------------------------------------------------------------
	% Perform initial exclusion based on gross movement
	% ------------------------------------------------------------------------------
	ParticipantIDs(exclude(:,1)) = [];
	Group_Diagnostic(exclude(:,1)) = [];
	fdJenk_m(exclude(:,1)) = [];
	fdJenk(exclude(:,1)) = [];
	Cov(exclude(:,1),:) = [];

	% compute numsubs
	numSubs = length(ParticipantIDs);

	% ------------------------------------------------------------------------------
	% Power's scrubbing
	% ------------------------------------------------------------------------------
	if ismember('scrub',runCensor,'rows') == 1
		fprintf(1, 'Loading scrubbing mask\n');
		[ScrubMask,~,~,~] = GetScrubbingForSample(datadir,ParticipantIDs,preprostr);
		fprintf(1, 'done\n');
	end

	% ------------------------------------------------------------------------------
	% Create exclusion vector based on <4 minutes of post-scrub BOLD data
	% Note, subjects not excluded yet. that is done during the preprocessing loop below
	% ------------------------------------------------------------------------------
	if ismember('false',runCensor,'rows') == 0
		% exclusion vector
		excludeCensor = zeros(numSubs,1);
		% threshold for exclusion in minutes
		thresh = 4;

		% Loop over subjects and check if censoring would reduce data to <4 minutes
		for i = 1:numSubs
			% Find number of vols left after censoring
			if ismember('scrub',runCensor,'rows') == 1
				% number of volumes - number of scrubbed volumes
				numCVols = numVols - sum(ScrubMask{i});
			elseif ismember('spikeReg',runCensor,'rows') == 1
				% Get spike regressors for subject i
				spikereg = GetSpikeRegressors(fdJenk{i},0.25);
				% number of volumes - number of spike regressors (columns)
				numCVols = numVols - size(spikereg,2);
			end	

			% Compute length, in minutes, of time series data left after censoring
			NTime = (numCVols * TR)/60;
			% if less than threshold, mark for exclusion
			if NTime < thresh
				excludeCensor(i) = 1;
			end
		end
		% convert to logical
		excludeCensor = logical(excludeCensor);

		ParticipantIDs(excludeCensor) = [];
		Group_Diagnostic(excludeCensor) = [];
		fdJenk_m(excludeCensor) = [];
		fdJenk(excludeCensor) = [];
		Cov(excludeCensor,:) = [];

		% compute numsubs
		numSubs = length(ParticipantIDs);

		fprintf(1, '%u subjects excluded based on 4 minute censoring criteria \n', sum(excludeCensor));
	end

	% ------------------------------------------------------------------------------
	% Split subject's
	% ------------------------------------------------------------------------------
	switch WhichSplit
		case 'Motion'
		% 1) based on motion

			% retain only one group (HCs assumed to be denoted by value 1)
			ParticipantIDs = ParticipantIDs(Group_Diagnostic == 1);
			fdJenk = fdJenk(Group_Diagnostic == 1);
			fdJenk_m = fdJenk_m(Group_Diagnostic == 1);
			Cov = Cov(Group_Diagnostic == 1,:);

			% Sort movement
			[~,idx] = sort(fdJenk_m,'ascend');

			% rearrange variables according to motion
			ParticipantIDs = ParticipantIDs(idx);
			fdJenk = fdJenk(idx);
			fdJenk_m = fdJenk_m(idx);
			Cov = Cov(idx);

			% recompute numsubs
			numSubs = length(ParticipantIDs);

			% create grouping variable
			% this will ensure that group 1 and 3 are the same size
			% group 2 might be a slightly different size, but contrasts are run on 1 and 3
			thirdPoint = floor(numSubs/3);
			Group = zeros(numSubs,1);
			Group(1:thirdPoint) = 1;
			Group(end-thirdPoint+1:end) = 3;
			Group(Group == 0) = 2;

			% Select which groups to perform t-tests on
			g1 = 1; % Low motion
			g2 = 3; % High motion
		case 'Diagnostic'
		% 2) Case vs HCs
			Group = Group_Diagnostic;
			% Select which groups to perform t-tests on
			g1 = 1; % HCs
			g2 = 2; % non-HCs
	end

	% ------------------------------------------------------------------------------
	% Load in time series data
	% ------------------------------------------------------------------------------
	fprintf(1, 'Loading time series data: %s\n',WhichNoise);
	cfg = [];
	for j = 1:numSubs
	    tsdir = [datadir,ParticipantIDs{j},preprostr,WhichNoise,'/'];
	    
	    clear temp
	    temp = load([tsdir,'cfg.mat']);
	    
	    cfg = [cfg temp.cfg];
	end
	fprintf(1, 'done\n');

	% ------------------------------------------------------------------------------
	% Compute correlations
	% ------------------------------------------------------------------------------
	numROIs = size(cfg(1).roiTS{Parc},2);

	FC = zeros(numROIs,numROIs,numSubs); 
	for j = 1:numSubs
		FC(:,:,j) = corr(cfg(j).roiTS{Parc});
		% Perform fisher z transform
		FC(:,:,j) = atanh(FC(:,:,j));
	end

	% ------------------------------------------------------------------------------
	% String for file names
	% ------------------------------------------------------------------------------
	nameStr = [WhichSplit,'_',WhichParc,'_',WhichNoise];

	% ------------------------------------------------------------------------------
	% N per group
	% ------------------------------------------------------------------------------
	numG1 = sum(Group == g1);
	numG2 = sum(Group == g2);

	% ------------------------------------------------------------------------------
	% Calculate primary threshold
	% ------------------------------------------------------------------------------
	df = numG1 + numG2 - 2;
	Tval = tinv(1 - P,df);
	% Tval = tinv(1 - 0.05,df);
	% Tval = tinv(1 - 0.01,df);
	% Tval = tinv(1 - 0.001,df);

	% ------------------------------------------------------------------------------
	% File names
	% ------------------------------------------------------------------------------
	matricesName = [nameStr,'_matrices_t',num2str(Tval),'.mat'];
	designName = [nameStr,'_designMatrix_t',num2str(Tval),'.mat'];
	outName = [nameStr,'_NBS_t',num2str(Tval),'.mat'];

	% ------------------------------------------------------------------------------
	% Generate 3D FC matrix for HC vs OCD
	% ------------------------------------------------------------------------------
	Mat = cat(3,FC(:,:,Group == g1),FC(:,:,Group == g2));
	save(matricesName,'Mat')

	% ------------------------------------------------------------------------------
	% Add tDOF as an additional covariate for ICA-AROMA or aCC50 pipeline
	% ------------------------------------------------------------------------------
	WhichNoiseSplit = strsplit(WhichNoise,'+');

	if any(strmatch('aCC50',WhichNoiseSplit,'exact')) == 1 | ...
		any(strmatch('sICA-AROMA',WhichNoiseSplit,'exact')) == 1
		fprintf(1, 'Computing tDOF: %s\n',WhichNoise);
		
		tDOF = zeros(numSubs,1);
		for j = 1:numSubs
			% Get noiseTS
			x = dlmread([datadir,ParticipantIDs{j},preprostr,WhichNoise,'/noiseTS.txt']);
			tDOF(j) = size(x,2);
				        		
			% If ICA-AROMA, get that too.
	        if any(strmatch('sICA-AROMA',WhichNoiseSplit,'exact')) == 1
				y = dlmread([datadir,ParticipantIDs{j},preprostr,WhichNoise,'/classified_motion_ICs.txt']);
				tDOF(j) = tDOF(j) + size(y,2);
			end
		end

		% mean center
		tDOF = tDOF - mean(tDOF);
		% add to Cov
		Cov = [Cov tDOF];
		% update zero pad
		zeroPad = [zeroPad,',0'];
	end

	% ------------------------------------------------------------------------------
	% Generate design matrix
	% ------------------------------------------------------------------------------
	designMatrix = [[ones(numG1,1); zeros(numG2,1)], [zeros(numG1,1); ones(numG2,1)], [Cov(Group == g1,:); Cov(Group == g2,:)]];
	save(designName,'designMatrix')

	% ------------------------------------------------------------------------------
	% Loop over NBS twice, once for each t contrast
	% ------------------------------------------------------------------------------
	nbsOut = cell(1,2);
	for j = 1:2
		% for reproducibility of t-contrasts in isolationg
		rng('default')

		fprintf(1, 'Running NBS for contrast %u\n', j);
		% ------------------------------------------------------------------------------
		% Run NBS with default settings (based on SCZ example in manual)
		% ------------------------------------------------------------------------------
		clear UI
		UI.method.ui = 'Run NBS'; 
		UI.test.ui = 't-test';
		UI.size.ui = 'Extent';
		UI.thresh.ui = num2str(Tval);
		UI.perms.ui = '5000';
		UI.alpha.ui = '0.05';
		if j == 1
			% Group 1 > Group 2
			UI.contrast.ui = ['[1,-1',zeroPad,']'];
		elseif j == 2
			% Group 1 < Group 2
			UI.contrast.ui = ['[-1,1',zeroPad,']'];
		end
		UI.matrices.ui = matricesName;
		UI.design.ui = designName;
		UI.exchange.ui = ''; 
		UI.node_label.ui = '';
		switch WhichParc
			case 'Gordon'
				UI.node_coor.ui = '/gpfs/M2Home/projects/Monash076/Linden/ROIs/Gordon/Gordon_Centroids.txt';
				% UI.node_coor.ui = '~/Dropbox/Work/ROIs/Gordon/Gordon_Centroids.txt';
			case 'Power'
				UI.node_coor.ui = '/gpfs/M2Home/projects/Monash076/Linden/ROIs/Power/Power2011_xyz_MNI.txt';
				% UI.node_coor.ui = '~/Dropbox/Work/ROIs/Power/Power2011_xyz_MNI.txt';
			end

		S = [];
		NBSrun(UI,S)

		% ------------------------------------------------------------------------------
		% Store and clean up
		% ------------------------------------------------------------------------------
		global nbs

		nbsOut{j} = nbs;

		clear nbs
		clear global

	end

	% ------------------------------------------------------------------------------
	% Save
	% ------------------------------------------------------------------------------
	save(outName,'nbsOut','P','Tval','df','Mat','designMatrix')

	% ------------------------------------------------------------------------------
	% Clean up
	% ------------------------------------------------------------------------------
	delete(matricesName,designName)
end
