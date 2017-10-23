%% ThePlot:
function [] = ThePlot(plotTitle,mov,fdPower,ScrubMask,fdJenk,dvars,ts_compartment,key_compartment,TR,movThr,fdPowerThr,fdJenkThr,dvarsThr)

	% This function plots a series of movement traces over top of a plot of timeseries as
	% as in Power (2016), called 'ThePlot'
	% It will overlay dashed horizontal lines on the movement parameter traces and fdPower traces
	% according to movThr and fdPowerThr, respectively
	%
	% ------
	% INPUTS
	% ------
	% plotTitle			- a string.
	% 					Note, this is just for naming the output .png file
	% mov 				- an numVols*6 matrix containing movement parameters extracted from SPM8's realignment
	% fdPower			- an numVols*1 vector containing Power's framewise displacement (see GetFDPower.m)
	% fdJenk			- an numVols*1 vector containing Jenkinson's framewise displacement (see GetFDJenk.m)
	% dvars 			- an numVols*1 vector containg dvars (see GetDVARS.m)
	% 
	% ts_compartment 	- a numVols * numVoxels timeseries matrix compartmentalised by
	% 					grey/white/csf (see GetTSCompartment.m)
	% key_compartment 	- a vector denoting which compartment a voxel belongs to (see GetTSCompartment.m)
	% 
	% TR 				- Integer representing EPI acquisition time in seconds. This is used to estimate the amount
	% 					time, in minutes, would remain following volume censoring and to mark plotTitles for exclusion
	% 					if <4 minutes remain
	% 
	% movThr 			- cut off for movement parameters. default = 2
	% 					Note, only for visualisation
	% fdPowerThr 		- cut off for fdPower trace. default = 0.2
	% fdJenkThr 		- cut off for fdJenk trace. default = 0.25
	% -------
	% OUTPUTS
	% -------
	% A pretty plot...
	%
	% Linden Parkes, Brain & Mental Health Laboratory, 2016
	% ------------------------------------------------------------------------------	

	if nargin < 10
		movThr = 2;
	end

	if nargin < 11
		fdPowerThr = 0.2;
	end

	if nargin < 12
		fdJenkThr = 0.25;
	end

	if nargin < 13
		dvarsThr = 20; % 2% signal change
	end

	% threshold for exclusion in minutes
	thresh = 4;

    numVols = size(ts_compartment,1);

	% ------------------------------------------------------------------------------
	% Plot
	% ------------------------------------------------------------------------------
	FSize = 10;
	h1 = figure('color','w', 'units', 'centimeters', 'pos', [0 0 21 29.7], 'name',['ThePlot: ',plotTitle]); box('on'); movegui(h1,'center');
	h2 = suptitle(plotTitle);
    pos = get(h2,'Position');
    set(h2,'Position',[pos(1)*1, pos(2)*0.5, pos(3)*1]);

	% ------------------------------------------------------------------------------
	% Movement: translation
	% ------------------------------------------------------------------------------
	sp1 = subplot(6,2,1);
    pos1 = get(sp1,'Position');
	plot(mov(:,1),'LineWidth',1.5);
	hold on
	plot(mov(:,2),'g','LineWidth',1.5);
	plot(mov(:,3),'r','LineWidth',1.5);
	title('translation','fontweight','bold')
	ylabel('mm')
	legend({'x','y','z'},'Orientation','horizontal','Location','best')
	legend('boxoff')
	xlim([1 numVols])
    set(sp1,'XTickLabel','');
	set(gca,'FontSize',FSize)

	% overlay threshold line
	if any(max(abs(mov(:,1:3))) > movThr)
		line([0 numVols],[movThr movThr],'LineStyle','--','Color','k')
		line([0 numVols],[-movThr -movThr],'LineStyle','--','Color','k')
		ylim([-movThr-1 movThr+1])
	end
	
	% ------------------------------------------------------------------------------
	% Movement: rotation
	% ------------------------------------------------------------------------------
	% First convert to mm using Power's approach
	mov(:,4:6) = 50*pi/180*mov(:,4:6);

	sp2 = subplot(6,2,3);
    pos2 = get(sp2,'Position');
	plot(mov(:,4),'LineWidth',1.5);
	hold on
	plot(mov(:,5),'g','LineWidth',1.5);
	plot(mov(:,6),'r','LineWidth',1.5);
	title('rotation','fontweight','bold');
	ylabel('mm')
	legend({'pitch','roll','yaw'},'Orientation','horizontal','Location','best')
	legend('boxoff')
	xlim([1 numVols])
    set(sp2,'XTickLabel','');
	set(gca,'FontSize',FSize)

	% overlay threshold line
	if any(max(abs(mov(:,4:6))) > movThr)
		line([0 numVols],[movThr movThr],'LineStyle','--','Color','k')
		line([0 numVols],[-movThr -movThr],'LineStyle','--','Color','k')
		ylim([-movThr-1 movThr+1])
	end

	% ------------------------------------------------------------------------------
	% FD Jenk
	% ------------------------------------------------------------------------------
	sp3 = subplot(6,2,5);
    pos3 = get(sp3,'Position');
	plot(fdJenk,'LineWidth',1.5);
	hold on
	title('fdJenk','fontweight','bold')
	ylabel('mm')
	xlim([1 numVols])
	ylim([0 max(fdJenk)+(max(fdJenk)*.25)])
    set(sp3,'XTickLabel','');
	set(gca,'FontSize',FSize)

	% Add text to plot
	fdJenk_m = mean(fdJenk); % mean FD
	fdJenkPerc = sum(fdJenk > fdJenkThr) / numVols * 100; % spike percentage
	spikereg = GetSpikeRegressors(fdJenk,fdJenkThr); % Spike regression exclusion
	numCVols = numVols - size(spikereg,2); % number of volumes - number of spike regressors (columns)
	NTime = (numCVols * TR)/60; % Compute length, in minutes, of time series data left after censoring
	
	str1 = ['Mean = ',num2str(fdJenk_m,'%0.2f'),'mm. Supra threshold spikes = ',num2str(fdJenkPerc,'%0.1f'),'%. '];
	str2 = [num2str(NTime,'%0.0f'),' mins of uncensored data.'];
	
	yLimits = ylim;
	text(round(numVols*.05),yLimits(2) - yLimits(2)*.15,[str1,str2],... 
					'HorizontalAlignment','left',...
					'VerticalAlignment','middle',...
					'Color','black',...
					'FontSize', FSize)

	% overlay threshold line
	if max(fdJenk) > fdJenkThr
		line([0 numVols],[fdJenkThr fdJenkThr],'LineStyle','--','Color','k','LineWidth',1.5);
	end

	% ------------------------------------------------------------------------------
	% FD Power
	% ------------------------------------------------------------------------------
	sp4 = subplot(6,2,7);
    pos4 = get(sp4,'Position');
	plot(fdPower,'LineWidth',1.5);
	hold on
	title('fdPower','fontweight','bold')
	ylabel('mm')
	xlim([1 numVols])
	ylim([0 max(fdPower)+(max(fdPower)*.25)])
    set(sp4,'XTickLabel','');
	set(gca,'FontSize',FSize)

	% Add text to plot
	fdPower_m = mean(fdPower); % mean FD
	fdPowerPerc = sum(fdPower > fdPowerThr) / numVols * 100; % spike percentage
	numCVols = numVols - sum(ScrubMask); % number of volumes - number of spike regressors (columns)
	NTime = (numCVols * TR)/60; % Compute length, in minutes, of time series data left after censoring

	str1 = ['Mean = ',num2str(fdPower_m,'%0.2f'),'mm. Supra threshold spikes = ',num2str(fdPowerPerc,'%0.1f'),'%. '];
	str2 = [num2str(NTime,'%0.0f'),' mins of uncensored data.'];
	
	yLimits = ylim;
	text(round(numVols*.05),yLimits(2) - yLimits(2)*.15,[str1,str2],... 
					'HorizontalAlignment','left',...
					'VerticalAlignment','middle',...
					'Color','black',...
					'FontSize', FSize)

	% overlay threshold line
	if max(fdPower) > fdPowerThr
		line([0 numVols],[fdPowerThr fdPowerThr],'LineStyle','--','Color','k','LineWidth',1.5);
	end

	% ------------------------------------------------------------------------------
	% DVARS
	% ------------------------------------------------------------------------------
	dvars = dvars/10;
	dvarsThr = dvarsThr/10;

	sp5 = subplot(6,2,9);
    pos5 = get(sp5,'Position');
	plot(dvars,'LineWidth',1.5);
	title('dvars','fontweight','bold')
	ylabel('rms signal change')
	xlim([1 numVols])
	ylim([0 max(dvars)+(max(dvars)*.10)])
    set(sp5,'XTickLabel','');
	set(gca,'FontSize',FSize)

	% overlay threshold line
	if max(dvars) > dvarsThr
		line([0 numVols],[dvarsThr dvarsThr],'LineStyle','--','Color','k','LineWidth',1.5);
	end

	% ------------------------------------------------------------------------------
	% Time series
	% ------------------------------------------------------------------------------
	sp6 = subplot(6,2,11);
    pos6 = get(sp6,'Position');
	imagesc(ts_compartment')
	colormap(gray)
	caxis([0 1])
	title('Time Series','fontweight','bold')
	ylabel('WHITE                            GREY') % Laziness for the win...
	xlabel('time (volumes)')
	sp6_1 = colorbar;
    pos6_1 = get(sp6_1,'Position');
	set(gca,'FontSize',FSize)

	sp7 = subplot(6,2,12);
    pos7 = get(sp7,'Position');
	imagesc(key_compartment')
    set(sp7,'YTickLabel','','XTickLabel','','TickLength',[0,0]);
	set(gca,'FontSize',FSize)
	colormap(gray)
	caxis([0 1])

    % ------------------------------------------------------------------------------
    % Sizing
    % ------------------------------------------------------------------------------
	% [left bottom width height]
    set(sp1,'Position',[pos1(1)*.6, pos1(2)*1.02 pos1(3)*2.5, pos1(4)*1]);
    set(sp2,'Position',[pos2(1)*.6, pos2(2)*1.035, pos2(3)*2.5, pos2(4)*1]);
    set(sp3,'Position',[pos3(1)*.6, pos3(2)*1.05, pos3(3)*2.5, pos3(4)*1]);
    set(sp4,'Position',[pos4(1)*.6, pos4(2)*1.08, pos4(3)*2.5, pos4(4)*1]);
    set(sp5,'Position',[pos5(1)*.6, pos5(2)*1.125, pos5(3)*2.5, pos5(4)*1]);
    
    set(sp6,'Position',[pos6(1)*.6, pos6(2)*0.35, pos6(3)*2.5, pos6(4)*2]);
    set(sp6_1,'Position',[pos6_1(1)*2.325, pos6_1(2)*.85, pos6_1(3)*1, pos6_1(4)*1]);
    set(sp7,'Position',[pos7(1)*.08, pos7(2)*0.35, pos7(3)*.07, pos7(4)*2]);
end
