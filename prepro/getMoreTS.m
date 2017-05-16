% ------------------------------------------------------------------------------
% Set project variables
% ------------------------------------------------------------------------------
projdir = '/gpfs/M2Home/projects/Monash076/Linden/OCDPG/';
sublist = '/gpfs/M2Home/projects/Monash076/Linden/Sublists/OCDPG.txt';
datadir = [projdir,'data/'];

fileID = fopen(sublist);
ParticipantIDs = textscan(fileID,'%s');
ParticipantIDs = ParticipantIDs{1};
numSubs = length(ParticipantIDs);

% All pipelines
noiseOptions = {'6P','6P+2P','6P+2P+GSR','24P','24P+8P','24P+8P+4GSR','24P+8P+SpikeReg','24P+8P+4GSR+SpikeReg','12P+aCC','24P+aCC','12P+aCC50','24P+aCC50','24P+aCC+4GSR','24P+aCC50+4GSR','24P+aCC+SpikeReg','24P+aCC+4GSR+SpikeReg','sICA-AROMA+2P','sICA-AROMA+2P+SpikeReg','sICA-AROMA+GSR','sICA-AROMA+2P+GSR','sICA-AROMA+8P','sICA-AROMA+4GSR','sICA-AROMA+8P+4GSR'};
numPrePro = length(noiseOptions);

% Loop over noise correction options
for i = 1:numSubs
    subject = ParticipantIDs{i};
    preprodir = [datadir,subject,'/rfMRI/prepro/'];

    for j = 1:numPrePro
        removeNoise = noiseOptions{j};

        cd([preprodir,removeNoise])
        load('cfg.mat')
        
        % Parcellation file for time series extraction
        cfg.new_parcFiles = {'/gpfs/M2Home/projects/Monash076/Linden/ROIs/TriStri/Striatum.nii'};
        cfg.parcFiles{end+1} = cfg.new_parcFiles{1};

        cfg.new_parcWeightGM = {'no'};
        cfg.parcWeightGM{end+1} = cfg.new_parcWeightGM{1};

        % Loop over parcellation files
        for k = 1:length(cfg.new_parcFiles)
            % Set parcellation file
            cfg.parcFile = cfg.new_parcFiles{k};
            % set GM weight
            cfg.weightGM = cfg.new_parcWeightGM{k};
            % extract time series 
            cfg.roiTS{end+1} = prepro_extractTS_FSL(cfg);
        end

        cfg = rmfield(cfg,{'parcFile','weightGM','new_parcFiles','new_parcWeightGM'});
    
        % Save data
        save('cfg.mat','cfg')
        clear cfg
    end
end