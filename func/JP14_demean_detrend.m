% ------------------------------------------------------------------------------
% Originally authored by Jonathan Power (see Power et al., 2014. NeuroImage).
% This code was downloaded from: http://www.jonathanpower.net/2014-ni-motion-2.html
% 
% Adapted by Linden Parkes for use in Parkes et al., 2017. bioRxiv
% This produces identical output to matlab's detrend function except also outputs the betas
% Note, Power assumes time on second dimension, whereas we work with time on first dimension.
% ------------------------------------------------------------------------------
function [tempbold,tempbetas] = JP14_demean_detrend(data,tmask);

    tempbold = data;
	[vox ts] = size(tempbold); 

	if nargin < 2
		% Note that in this tmask, 1 = volumes to keep, not volumes to censor!
	    tmask = ones(ts,1);
	end
	
	% make logical
	tmask = logical(tmask);

	linreg = [repmat(1,[ts 1]) linspace(0,1,ts)'];
	tempboldcell = num2cell(tempbold(:,tmask)',1);
	linregcell = repmat({linreg(tmask,:)},[1 vox]);
	tempbetas = cellfun(@mldivide,linregcell,tempboldcell,'uniformoutput',0);
	tempbetas = cell2mat(tempbetas);
	tempbetas = tempbetas';
	tempintvals = tempbetas * linreg';
	tempbold = tempbold - tempintvals;

end
