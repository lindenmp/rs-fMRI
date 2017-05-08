%% TheTMatrixPlot: 
function TheTMatrixPlot(TMatrix,ROIStruct,ColorLabelsWhere)

% ------------------------------------------------------------------------------
% Check inputs:
% ------------------------------------------------------------------------------
if nargin < 2
    fprintf(1,'We need colors and labels for each element of the matrix\n');
end
if nargin < 3 || isempty(ColorLabelsWhere)
    ColorLabelsWhere = 'left'; % left, bottom, both
end

% ------------------------------------------------------------------------------
% Flip stuff (pcolor is weird).
% ------------------------------------------------------------------------------
TMatrix = flipud(TMatrix); % first labels at the top
ROIStruct = ROIStruct(length(ROIStruct):-1:1); % first labels at the top

% ------------------------------------------------------------------------------
% assign each ROI its color
% ------------------------------------------------------------------------------
[C,ia,ic] = unique(ROIStruct);
numStructures = numel(unique(ROIStruct));
% theColors = colormap(lines(numStructures));
theColors = BF_getcmap('set5',numStructures);

ROIStructColorMap = [];
for i = 1:numStructures
	temp = repmat(theColors(i,:),sum(ic == i),1);
	ROIStructColorMap = [ROIStructColorMap; temp];
end

% ------------------------------------------------------------------------------
% Plot
% ------------------------------------------------------------------------------
MyColorMap = [flipud(BF_getcmap('blues',9,0));1,1,1;BF_getcmap('reds',9,0)];
pcolor([TMatrix, zeros(size(TMatrix,1),1); zeros(1,size(TMatrix,2)+1)]);
shading flat
colormap(MyColorMap);
% maxT = max(abs(TMatrix(:)));
% caxis([-maxT maxT]);
caxis([-5 5]);
colorbar('Ticks',[-5,-4,-3,-2,-1,0,1,2,3,4,5])
hold on

% ------------------------------------------------------------------------------
% Add structure labels
% ------------------------------------------------------------------------------
rectThickness = arrayfun(@(x)size(TMatrix,x)/50,1:2);

% Rows:
for i = 1:length(ROIStruct)
    if ismember(ColorLabelsWhere,{'both','left'})
        rectangle('Position',[1-rectThickness(2),i,rectThickness(2),1], ...
                    'FaceColor',ROIStructColorMap(i,:),'EdgeColor',ROIStructColorMap(i,:))
    end
end

% ------------------------------------------------------------------------------
% Adjust axes to show labels
% ------------------------------------------------------------------------------
scaleFactor = 0.2; % to see the little bit extra to capture the line thickness
switch ColorLabelsWhere
case 'left'
    set(gca,'XLim',[-rectThickness(2)*(scaleFactor),size(TMatrix,2)+rectThickness(2)*scaleFactor])
end
if ismember(ColorLabelsWhere,{'both','left'})
    % get_xlim = get(gca,'xlim');
    % get_xlim(1) = -rectThickness*(1+scaleFactor);
    set(gca,'XLim',[1-rectThickness(2)*(1+scaleFactor),size(TMatrix,2)+1]);
end

set(gca,'TickLength',[0,0])
