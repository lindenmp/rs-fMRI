%% GetTStats: 
function [T,GroupDiff,Denominator] = GetTStats(FCVec,Group,flipCon)
	if nargin < 3
		flipCon = false;
	end

	if flipCon
		% Split FC by group
		x = FCVec(Group == 2,:);
		y = FCVec(Group == 1,:);
	else
		% Split FC by group
		x = FCVec(Group == 1,:);
		y = FCVec(Group == 2,:);
	end

	% Get t stats
	[h,p,ci,stats] = ttest2(x,y,'Vartype','equal');
	T = stats.tstat;

	% Get t-test numerators
	xbar = mean(x,1);	
	ybar = mean(y,1);
	GroupDiff = xbar - ybar;
	GroupDiff = GroupDiff';

	% Get t-test denominators
	Sd = stats.sd;
	n1 = size(x,1);
	n2 = size(y,1);
	Denominator = Sd .* (sqrt(1/n1 + 1/n2));
	Denominator = Denominator';
