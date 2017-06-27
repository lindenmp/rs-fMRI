%% QC_SaveDVARS: 
function [] = QC_SaveDVARS(WhichMASSIVE,WhichProject)
	fprintf(1, 'Saving DVARS for %s.\n',WhichProject);

	% ------------------------------------------------------------------------------
	% Linden Parkes, Brain & Mental Health Laboratory, 2016
	% ------------------------------------------------------------------------------
	% clear all; close all; clc
	% WhichProject = 'OCDPG' % 'OCDPG' 'UCLA'

    % ------------------------------------------------------------------------------
    % Add paths - edit this section
    % ------------------------------------------------------------------------------
    switch WhichMASSIVE
        case 'M2'
            % where the prepro scripts are
            scriptdir = '/gpfs/M2Home/projects/Monash076/Linden/scripts/rs-fMRI/prepro/';
            addpath(scriptdir)
            funcdir = '/gpfs/M2Home/projects/Monash076/Linden/scripts/rs-fMRI/func/';
            addpath(funcdir)
        case 'M3'
            % where the prepro scripts are
            scriptdir = '/home/lindenmp/kg98/Linden/Scripts/rs-fMRI/prepro/';
            addpath(scriptdir)
            funcdir = '/home/lindenmp/kg98/Linden/Scripts/rs-fMRI/func/';
            addpath(funcdir)
    end

	% ------------------------------------------------------------------------------
	% Select project
	% ------------------------------------------------------------------------------
	switch WhichProject
		case 'OCDPG'
			projdir = '/gpfs/M2Home/projects/Monash076/Linden/OCDPG/';
			sublist = '/gpfs/M2Home/projects/Monash076/Linden/Sublists/OCDPG.txt';
			datadir = [projdir,'data/'];
			
			preprostr = '/rfMRI/prepro/';
	        EPI = 'epi.nii';
		case 'UCLA'
			projdir = '/gpfs/M2Home/projects/Monash076/Linden/UCLA/';
			sublist = '/gpfs/M2Home/projects/Monash076/Linden/Sublists/UCLA.txt';
			datadir = [projdir,'data/'];
	        
	        preprostr = '/rfMRI/prepro/';
	        EPI = 'epi.nii';
	    case 'NYU_2'
			projdir = '/gpfs/M2Home/projects/Monash076/Linden/NYU_2/';
			sublist = '/gpfs/M2Home/projects/Monash076/Linden/Sublists/NYU_2.txt';
			datadir = [projdir,'data/'];
	        
	        % preprostr = '/session_1/rest_1/prepro/';
	        preprostr = '/session_1/rest_2/prepro/';
	        % preprostr = '/session_2/rest_1/prepro/';
	        EPI = 'rest.nii';
		case 'M3_COBRE'
			projdir = '/home/lindenmp/kg98/Linden/ResProjects/SCZ_HCTSA/COBRE/';
			sublist = [projdir,'COBRE_SubjectIDs.txt'];
			datadir = [projdir,'data/'];
			
			preprostr = '/session_1/rest_1/prepro/';
	        EPI = 'rest.nii';
		case 'M3_UCLA'
			projdir = '/home/lindenmp/kg98/Linden/ResProjects/SCZ_HCTSA/UCLA/';
			sublist = [projdir,'UCLA_SubjectIDs.txt'];
			datadir = [projdir,'data/'];
			
			preprostr = '/func/prepro/';
		case 'M3_NAMIC'
			projdir = '/home/lindenmp/kg98/Linden/ResProjects/SCZ_HCTSA/NAMIC/';
			sublist = [projdir,'NAMIC_SubjectIDs.txt'];
			datadir = [projdir,'data/'];
			
			preprostr = '/func/prepro/';
	end

	% ------------------------------------------------------------------------------
	% Subject list
	% ------------------------------------------------------------------------------
	fileID = fopen(sublist);
	ParticipantIDs = textscan(fileID,'%s');
	ParticipantIDs = ParticipantIDs{1};

	numSubs = length(ParticipantIDs);

	% ------------------------------------------------------------------------------
	% Compute DVARS
	% NOTE: Important that the epi used for this is the one BEFORE major noise correction
	% Here I use the smoothed detrended epi output from core image preprocessing
	% ------------------------------------------------------------------------------
	for i = 1:numSubs
		fprintf(1,'Processing subject %u/%u: %s\n',i,numSubs,ParticipantIDs{i})

		epidir = [datadir,ParticipantIDs{i},preprostr];
		cd(epidir)

		if ismember('M3_UCLA',WhichProject,'rows') | ismember('M3_NAMIC',WhichProject,'rows')
			EPI = [ParticipantIDs{i},'_task-rest_bold.nii'];
		end
		
		dvarsExtract = ['idbwrdat',EPI];
		dvarsMask = 'epi_brain_mask.nii';

		% ------------------------------------------------------------------------------
		% Decompress if .gz
		% ------------------------------------------------------------------------------
	    % Decompress files
	    filesToCheck = {dvarsExtract, ...
	                    dvarsMask};

	    dirsOfFiles = {epidir, ...
	                    epidir};

	    for j = 1:length(filesToCheck)
	        if exist([dirsOfFiles{j},filesToCheck{j},'.gz']) == 2
	            fprintf(1, '\t\t Decompressing %s\n',filesToCheck{j});
	            gunzip([dirsOfFiles{j},filesToCheck{j},'.gz'],dirsOfFiles{j})
	            delete([dirsOfFiles{j},filesToCheck{j},'.gz'])
	        elseif exist([dirsOfFiles{j},filesToCheck{j},'.gz']) == 0
	            fprintf(1, '\t\t No need to decompress %s\n',filesToCheck{j});
	        end
	    end

        % ------------------------------------------------------------------------------
	    % Get dvars
        % ------------------------------------------------------------------------------
	    dvars = GetDVARS(dvarsExtract,dvarsMask);

	    % Save dvars out
	    save('dvars.mat','dvars','dvarsExtract','dvarsMask')

	    clear dvars

	    % ------------------------------------------------------------------------------
	    % Clean up files
	    % ------------------------------------------------------------------------------
		% Recompress files
	    for j = 1:length(filesToCheck)
	        fprintf(1, '\t\t Recompressing %s\n',filesToCheck{j});
	        gzip([dirsOfFiles{j},filesToCheck{j}],dirsOfFiles{j})
	        delete([dirsOfFiles{j},filesToCheck{j}])
	    end

	end

	fprintf(1, 'Finished.\n');