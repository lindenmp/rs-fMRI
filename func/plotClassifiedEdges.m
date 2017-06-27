function [out outPC outNorm] = plotClassifiedEdges(adj, ids, plotFig, labels)

% function [out outPC] = plotClassifiedEdges(adj, ids, plotFig, labels)

% This function is useful for understanding the results of NBS output.
% Given a classification of nodes into a subset of networks (e.g., DMN,
% FPN, CON, etc), it counts the number of edges present in a network both
% within and between networks and optionally plots the output as a matrix.
%
% The results are counted only at the level of binary topology
%
% -------
% INPUTS:
% -------
% adj       - binary N*N adjacency matrix, where N is number of nodes. This matrix
%           could be the output of an NBS analysis; i.e., the matrix stored
%           in nbs.NBS.con_may{1}.
%
% ids       - N*1 vector of network ids. Each network should be represented
%           as a unique number. Each row is a different node.
%
% plotFig   - set to 1 if you want to plot out; 2 for outPC; 3 for outNorm; 
%           0 otherwise. Default is 0.
%
% labels    - M*1 cell, where M is the number of different networks. Each 
%             cell contains a string representing the
%            network name. This is used to label the axes of the plot.
%
% -------
% OUTPUTS:
% -------
% out       - M*M adjacency matrix indicating the number of edges between
%           each network pair. Values along diagonal represent edges within
%           a network. 
%           NOTE: sum(sum(triu(out))) = sum(sum(adj))/2; in other 
%           words, the number of unique edges in adj equals the number of
%           edges in the upper triangle (including the diagonal) of out.
%
% outPC     - The out matrix normalized by the number of unqiue edge to
%           yield the proportion of edges for each category.
%
% outNorm   - this matrix is normalized separetly for each pair of regions
%           by the total number of edges between them. The values thus
%           represent the connection density of the subgraph of nodes
%           belonging to a given pair of modules. This normalization
%           accounts for differences in the size of modules, which can bias
%           the results of outPC.
%
%
% Alex Fornito, Monash University, Oct 2016
%
% Updates: 
% Feb 2017 added outNorm as output; options to plot different output
% matrices
%
%==========================================================================

% debugging
% adj = nbsAdj;
% ids = att.ci;
% plotFig = 0;
% labels = modLabels;

%==========================================================================
% Preliminaries
%==========================================================================

% set defaults
if nargin<3

    plotFig = 0;
    labels = {};
    
elseif nargin<4
    
    labels = {};

end

% if NBS output is used directly, convert matrix to double and symmetrize
if issparse(adj)==1
    
    adj = full(adj);
    adj = adj+adj';
    
end

% binarize
adj = double(logical(adj));

%==========================================================================
% Compute outputs
%==========================================================================

N = size(adj,1);                            % number of nodes
unq = unique(ids);                          % number of unique networks
out = zeros(length(unq));                   % initalize out

% sort edges
for i = 1:length(unq)
    
    for j = i:length(unq)
        
        inds_i = find(ids == i);
        inds_j = find(ids == j);

        edgeCount = sum(sum(adj(inds_i, inds_j)));

        if i==j
            edgeCount = edgeCount/2;                    % halve because diagnoal gets counted twice
        end
               
        out(i,j) = edgeCount;
        out(j,i) = edgeCount;
        
        nSubTot = length(inds_i) + length(inds_j);
        normFactor = (nSubTot^2 - nSubTot)/2;
        
        outNorm(i,j) = edgeCount/normFactor;
        outNorm(j,i) = edgeCount/normFactor;
        
    end
end

% get proportions
outPC = zeros(size(out,1));
outPC = out./(nnz(adj)/2);

%==========================================================================
% Plot output
%==========================================================================

if plotFig > 0 

    if plotFig == 1
        plotMat = out;
    elseif plotFig == 2
        plotMat = outPC;
    elseif plotFig == 3
        plotMat = outNorm;
    end
    plotMat = plotMat * 100;
    
    cmax = ceil(max(plotMat(:)));
    imagesc(plotMat, [0, cmax]); % Plot matrix, with limits between 0 and max val in nmatrix
        axis square
        axis tight
        colormap(BF_getcmap('reds',6,0))
        % caxis([0 900])
        caxis([0 cmax])
    % colormap(hot);
    hold on
    set(gca,'XTick',1:1:length(labels)); % Set x tick vals
    set(gca,'TickLength', [0 0])
    set(gca,'XTickLabel',labels,'FontSize',10,'FontWeight','normal','XTickLabelRotation',45); % set x tick labels
    set(gca,'YTick',1:1:length(labels)); % set y tick vals
    set(gca,'YTickLabel',labels,'FontSize',10,'FontWeight','normal'); % set y tick labels
    colorbar; % show colour bar

    % plot values in lower triangle of matrix
    % for i=1:size(plotMat,1)
    %    for j=i:size(plotMat,2)
    %         text(i,j,num2str(plotMat(i,j),'%0.0f'),'HorizontalAlignment','center',...
    %             'Color','k','FontSize',10,'FontWeight','normal');
    %    end
    % end

end



        