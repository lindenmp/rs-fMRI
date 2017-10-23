function [] = TheBarChart(data,data_std,makeFigure,extraParams)
% Linden Parkes, BMH, 2016

	% ------------------------------------------------------------------------------
	% set default args

	if nargin < 2
		for i = 1:length(data)
			data_std{i} = zeros(size(data{i}));
		end
	end

	if nargin < 3
	    makeFigure = true;
	end

	if nargin < 4
	    extraParams = struct;
	end

	% ------------------------------------------------------------------------------
	% Set extra parameters

	% Custom xTickLabels
	if ~isfield(extraParams,'xTickLabels')
		xTickLabels = [1:length(data{1})];
	else
	    xTickLabels = extraParams.xTickLabels;
	end

	% Custom XTickLabelRot
	if ~isfield(extraParams,'XTickLabelRot')
	    XTickLabelRot = 0;
	else
	    XTickLabelRot = extraParams.XTickLabelRot;
	end

	% Custom YTickLabelRot
	if ~isfield(extraParams,'YTickLabelRot')
	    YTickLabelRot = 0;
	else
	    YTickLabelRot = extraParams.YTickLabelRot;
	end

	% Custom xLimits
	if ~isfield(extraParams,'xLimits')
	    xLimits = [0 size(data{1},1)+1];
	else
	    xLimits = extraParams.xLimits;
	end

	% Custom yLimits
	if ~isfield(extraParams,'yLimits')
	    if length(data) == 1
	    	yLimits = [0 100];
	    elseif length(data) == 2
	    	yLimits = [-100 100];
		end
	else
	    yLimits = extraParams.yLimits;
	end

	% Custom Yaxisdir
	if ~isfield(extraParams,'Yaxisdir')
	    Yaxisdir = 'Normal';
	else
	    Yaxisdir = extraParams.Yaxisdir;
	end

	% Custom xLabel
	if ~isfield(extraParams,'xLabel')
	    xLabel = '';
	else
	    xLabel = extraParams.xLabel;
	end	

	% Custom yLabel
	if ~isfield(extraParams,'yLabel')
	    yLabel = '';
	else
	    yLabel = extraParams.yLabel;
	end	

	% Custom Title
	if ~isfield(extraParams,'Title')
	    Title = '';
	else
	    Title = extraParams.Title;
	end

	% Custom plotWidth
	if ~isfield(extraParams,'plotWidth')
	    plotWidth = 10.5; % width of an A4 page
	    % plotWidth = 21; % width of an A4 page
	else
	    plotWidth = extraParams.plotWidth;
	end

	% Custom plotHeight
	if ~isfield(extraParams,'plotHeight')
	    plotHeight = 10;
	else
	    plotHeight = extraParams.plotHeight;
	end

	% Custom FSize
	if ~isfield(extraParams,'FSize')
	    FSize = 10;
	else
	    FSize = extraParams.FSize;
	end

	% absolute y axis
	if ~isfield(extraParams,'makeABS')
		makeABS = false;
	else
		makeABS = extraParams.makeABS;
	end

	% add text
	if ~isfield(extraParams,'addText')
		addText = true;
	else
		addText = extraParams.addText;
	end

	if ~isfield(extraParams,'theColors')
		theColors = num2cell(repmat([0.3 0.3 0.3],numel(data{1}),1),2);
	else
		theColors = extraParams.theColors;
	end

	if ~isfield(extraParams,'theLines')
		[theLines{1:numel(data{1})}] = deal('-');
	else
		theLines = extraParams.theLines;
	end

	% ------------------------------------------------------------------------------

	% ------------------------------------------------------------------------------
	% Create figure and get axes
	% ------------------------------------------------------------------------------
	if makeFigure
		fHand = figure('color','w', 'units', 'centimeters', 'pos', [0 0 plotWidth plotHeight], 'name',['']);
		movegui(fHand,'center');
	end
	hold on
	box('on');
	
	% ------------------------------------------------------------------------------
	% Plot data
	% ------------------------------------------------------------------------------

	if length(data) == 1 || length(data) == 2
		% Plot bar chart 1 above 0 point on y-axis
		if size(data{1},2) == 1
			for i = 1:size(data{1},1)
				bar(i,data{1}(i),'FaceColor',theColors{i},'LineStyle',theLines{i},'LineWidth',2);
			end
		elseif size(data{1},2) == 3
			b1 = bar(data{1});
			b1(1).FaceColor = [0 0 0];
			b1(2).FaceColor = [0.5 0.5 0.5];
			b1(3).FaceColor = [1 1 1];			
		end

		if any(data_std{1})
			% eb1 = errorbar(data{1},data_std{1},'k.','MarkerFaceColor',[0.3 0.3 0.3],'MarkerEdgeColor',[0.3 0.3 0.3]);
			eb1 = errorbar(data{1},data_std{1},'k.','LineWidth',2);
		end

		if length(data) == 2
			% Plot bar chart 2 below 0 point on y-axis
			if size(data{2},2) == 1
				for i = 1:numel(data{2})
					bar(i,-data{2}(i),'FaceColor',theColors{i},'LineStyle',theLines{i},'LineWidth',2);
				end
			elseif size(data{2},2) > 1
				b2 = bar(-data{2});
				b2(1).FaceColor = [0 0 0];
				b2(2).FaceColor = [0.5 0.5 0.5];
				b2(3).FaceColor = [1 1 1];
			end

			if any(data_std{2})
				eb2 = errorbar(-data{2},-data_std{2},'k.','LineWidth',2);
			end
		end
	else
		error('TheBarChart only supports cells where length <=2')
	end

	ax = gca;

	% Set axis stuff
	ax.FontSize = FSize;
	XTicks = [1:size(data{1},1)];
	ax.XTick = XTicks;
	ax.XTickLabel = xTickLabels;
	ax.XTickLabelRotation = XTickLabelRot;
	ax.XLim = (xLimits);
	ax.YTickLabelRotation = YTickLabelRot;
	if yLimits(2) <= 0.5
		ax.YTick = ([yLimits(1):0.5:yLimits(2)]);
		% ax.YTick = ([yLimits(1):0.1:yLimits(2)]);
	elseif yLimits(2) <= 1
		ax.YTick = ([yLimits(1):0.5:yLimits(2)]);
	elseif yLimits(2) <= 20
		ax.YTick = ([yLimits(1):10:yLimits(2)]);
	elseif yLimits(2) == 25
		ax.YTick = ([yLimits(1):5:yLimits(2)]);
	elseif yLimits(2) <= 40
		ax.YTick = ([yLimits(1):20:yLimits(2)]);
	else
		% ax.YTick = ([yLimits(1):yLimits(2)]);
	end
	ax.YLim = (yLimits);
	ax.YDir = Yaxisdir;

	xlabel(xLabel)
	ylabel(yLabel)
	title(Title,'FontSize',FSize+2,'FontWeight','normal')

	% add legend
	if isfield(extraParams,'Legend')
		lh = legend(extraParams.Legend,'Location','northoutside');
	end

	% add text
	TextRotation = 0;

	% if yLimits(2) <= 0.5
	% 	strprec = '%0.3f';
	% elseif yLimits(2) <= 1 & yLimits(2) > 0.5
	% 	strprec = '%0.2f';
	% elseif yLimits(2) > 1
	% 	strprec = '%0.1f';
	% end
	strprec = '%0.2f';

	if addText == true
		if length(data) == 1 || length(data) == 2
			if size(data{1},2) == 1
				text(XTicks,repmat(yLimits(2) - yLimits(2)*.05,1,size(data{1},1)),num2str(data{1},strprec),... 
				'HorizontalAlignment','right',... 
				'VerticalAlignment','middle',...
				'Color','black',...
				'FontSize', FSize,...
				'Rotation',TextRotation)
			elseif size(data{1},2) == 2
				offset = [-0.125 0.125];
				for i = 1:2
					text(XTicks+offset(i),repmat(yLimits(2) - yLimits(2)*.05,1,size(data{1},1)),num2str(data{1}(:,i),strprec),... 
					'HorizontalAlignment','right',... 
					'VerticalAlignment','middle',...
					'Color','black',...
					'FontSize', FSize,...
					'Rotation',TextRotation)
				end	
			elseif size(data{1},2) == 3
				offset = [-0.225 0 0.225];
				for i = 1:3
					text(XTicks+offset(i),repmat(yLimits(2) - yLimits(2)*.05,1,size(data{1},1)),num2str(data{1}(:,i),strprec),... 
					'HorizontalAlignment','right',... 
					'VerticalAlignment','middle',...
					'Color','black',...
					'FontSize', FSize,...
					'Rotation',TextRotation)
				end
			end
		end

		if length(data) == 2
			if size(data{2},2) == 1
				text(XTicks,repmat(yLimits(1) - yLimits(1)*.05,1,size(data{2},1)),num2str(data{2},strprec),... 
				'HorizontalAlignment','left',... 
				'VerticalAlignment','middle',...
				'Color','black',...
				'FontSize', FSize,...
				'Rotation',TextRotation)
			elseif size(data{2},2) == 2
				offset = [-0.125 0.125];
				for i = 1:2
					text(XTicks+offset(i),repmat(yLimits(1) - yLimits(1)*.05,1,size(data{2},1)),num2str(data{2}(:,i),strprec),... 
					'HorizontalAlignment','left',... 
					'VerticalAlignment','middle',...
					'Color','black',...
					'FontSize', FSize,...
					'Rotation',TextRotation)
				end	
			elseif size(data{2},2) == 3
				offset = [-0.225 0 0.225];
				for i = 1:3
					text(XTicks+offset(i),repmat(yLimits(1) - yLimits(1)*.05,1,size(data{2},1)),num2str(data{2}(:,i),strprec),... 
					'HorizontalAlignment','left',... 
					'VerticalAlignment','middle',...
					'Color','black',...
					'FontSize', FSize,...
					'Rotation',TextRotation)
				end
			end
		end
	end

	% rotate
	view(90,90)

	% make y axis absolute
	if makeABS 
		ax.YTickLabel = abs(ax.YTick);
	end

