%% Computes the interclass correlations for a given functional connection
function ICC = GetICC(data,cse,typ)
% INPUT:
%   data - n data k (targets x raters) data matrix
%         in the context of fMRI, a target is a participant and a rater is their different fMRI scans
%         thus, n = participants and k = fmri scans
%         data should only contain information across subjects and scans for a SINGLE functional connection
%         Each target is assumed to be a randon sample from a population of targets.
%   cse - 1 2 or 3: 
%         1 = random effects model
%           1 if each target is measured by a different set of raters from a population of raters
%         2 = mixed effects model
%           2 if each target is measured by the same raters, but that these raters are sampled from a 
%           population of raters
%         3 = fixed effects model
%           3 if each target is measured by the same raters and 
%           these raters are the only raters of interest.
%   typ - 'single' or 'k': denotes whether the ICC is based on a single
%         measurement or on an average of k measurements, where 
%         k = the number of ratings/raters.
%    

if nargin < 2
    % for fMRI, typically we use random effects model since there is nothing inherently meaningful
    % about the scan labels (i.e., "first scan" or "second scan")
    cse = 1;
end
if nargin < 3
    typ = 'single';
end

[n,k] = size(data);

% mean per target
mpt = mean(data,2);
% mean per rater/rating
mpr = mean(data);
% get total mean
tm = mean(data(:));
% within target sum sqrs
tmp = (data - repmat(mpt,1,k)).^2;
WSS = sum(tmp(:));
% within target mean sqrs
WMS = WSS / (n*(k - 1));
% between rater sum sqrs
RSS = sum((mpr - tm).^2) * n;
% between rater mean sqrs
RMS = RSS / (k - 1);
% between target sum sqrs
BSS = sum((mpt - tm).^2) * k;
% between targets mean squares
BMS = BSS / (n - 1);
% residual sum of squares
ESS = WSS - RSS;
% residual mean sqrs
EMS = ESS / ((k - 1) * (n - 1));

switch cse
    case 1
        switch typ
            case 'single'
                ICC = (BMS - WMS) / (BMS + (k - 1) * WMS);
            case 'k'
                ICC = (BMS - WMS) / BMS;
            otherwise
               error('Wrong value for input typ') 
        end
    case 2
        switch typ
            case 'single'
                ICC = (BMS - EMS) / (BMS + (k - 1) * EMS + k * (RMS - EMS) / n);
            case 'k'
                ICC = (BMS - EMS) / (BMS + (RMS - EMS) / n);
            otherwise
               error('Wrong value for input typ') 
        end
    case 3
        switch typ
            case 'single'
                ICC = (BMS - EMS) / (BMS + (k - 1) * EMS);
            case 'k'
                ICC = (BMS - EMS) / BMS;
            otherwise
               error('Wrong value for input typ') 
        end
    otherwise
        error('Wrong value for input cse')
end