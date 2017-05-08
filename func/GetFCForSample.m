%% GetFCForSample: 
function [cfg,FC,FC_vec,VarCovar,Var,GCOR] = GetFCForSample(datadir,ParticipantIDs,str,removeNoise,cfgFile,Parc,numROIs,numConnections)
	% ------------------------------------------------------------------------------
	% 
	% Linden Parkes, Brain & Mental Health Laboratory, 2016
	% ------------------------------------------------------------------------------
	numSubs = length(ParticipantIDs);

	% ------------------------------------------------------------------------------
	% Containers
	% ------------------------------------------------------------------------------
    cfg = [];

	FC = zeros(numROIs,numROIs,numSubs);
	
	FC_vec = zeros(numSubs,numConnections);

	VarCovar = zeros(numROIs,numROIs,numSubs);

	Var = zeros(numROIs,numROIs,numSubs);

	GCOR = zeros(numSubs,1);

	% ------------------------------------------------------------------------------
	% Loop
	% ------------------------------------------------------------------------------
	for i = 1:numSubs
		% msg = sprintf('\tProcessing subject %u/%u: %s\n',i,numSubs,ParticipantIDs{i});
	    % fprintf(msg);

		% ------------------------------------------------------------------------------
		% Load in time series data
		% ------------------------------------------------------------------------------
		workdir = [datadir,ParticipantIDs{i},str,removeNoise,'/'];

		clear temp
		temp = load([workdir,cfgFile]);

		cfg = [cfg temp.cfg];

	    % ------------------------------------------------------------------------------
	    % Compute correlations
	    % ------------------------------------------------------------------------------
		FC(:,:,i) = corr(cfg(i).roiTS{Parc});
		% Perform fisher z transform
		FC(:,:,i) = fisherz(FC(:,:,i));
		% get GCOR
		GCOR(i) = mean(LP_FlatMat(FC(:,:,i)));
		% Store
		FC_vec(i,:) = LP_FlatMat(FC(:,:,i));

		% ------------------------------------------------------------------------------
		% Compute var/covariance matrix
		% ------------------------------------------------------------------------------
		VarCovar(:,:,i) = cov(cfg(i).roiTS{Parc});

		% ------------------------------------------------------------------------------
		% Compute correlation square root of variance product
		% ------------------------------------------------------------------------------
		X = var(cfg(i).roiTS{Parc});
		for j = 1:numROIs
			for jj = 1:numROIs
				Var(j,jj,i) = sqrt(X(j)*X(jj));
			end
		end

	    % n = numel(msg);
	    % fprintf(repmat('\b',1,n));
	end