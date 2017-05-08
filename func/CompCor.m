function [vecs pc_var p_val] = CompCor(data,method,Nperms,keep)

% This function implements CompCor as desceribed in Behzadi et al. (2007)
% NeuroImage. Given a matrix of input nuisance voxel time series, it performs a 
% temporal PCA and then selects relevnt PCs for noise regression according to 
% on of several options, variously suggested in the literature.
%
% -------
% INPUTS:
% -------
%
% data = a N x M matrix, where N=number of time points and M = number of
% voxels. This could contain time series of voles showing the highest
% (e.g., top 2SD) temporal variability, or derived from an anatomical mask
% (e.g., 99% probability of white or csf tissue).
%
% method = string determining the method used to select principal components:
%   'retain' will keep the first X components, where X is set by the keep
%   argument
%
%   '50pc' will keep the firsy X components accounting for 50% of the
%   variance. This was recommended in Muschelli et al. (2014) NeuroImage.
%
%   'perm' will do monte carlo testing to identify components accoutning
%   for a significant proportion of variance
%
% Nperms = the number of permutations to run for the monte-carlo
% simulations. Default=1000. Only relevant if method = 'perm'
%
% keep = scalar value indicating how many components to keep. Only valid if 
% method = 'retain'. 
%
% --------
% OUTPUTS:
% --------
%
% vecs = a N x M matrix containing the time courses of each retained
% component.
% N= number of time points; M = number of components.
%
% pc_var = % variance explained by each pc
%
% pval = the pvals for each component. Only valid if method = 'perm'.
%
% Alex Fornito, Monash University, April 2014
% =========================================================================

p_val = [];

if nargin<2
    method = '50pc';
    Nperms=1000;
	keep = 5;
end

if nargin<3
    Nperms = 1000;
	keep=5;
end

if nargin<4
    keep=5;
end

pval = [];

    % Remove linear trend
    dt=detrend(data);
    z = zscore(dt);
    [coeff score latent] = pca_stats(z);
    
    pc_var = latent./sum(latent);
    
 switch method
     
     case {'retain'}
         
         vecs = score(:,1:keep);
         
     case {'50pc'}
    
         cs = cumsum(pc_var);
         inds = find(cs>.50);
         vecs = score(:,1:min(inds));
         
     case {'perm'}

        %Do Monte Carlo simulation to generate null distribution
        for i=1:Nperms
            
            Rnd=randn(size(data)); 

            Rdt=detrend(Rnd);
            Rz = zscore(Rdt);
            [Rcoeff Rscore Rlatent] = pca_stats(Rz);

            null(:,i) = latent./sum(latent);
  
        end

        %Compare the actual principal values to the null distribution 
        for i=1:length(Rpc_var); 
            p_val(i)=length(find(null(:,i)>pc_var(i)))/Nperms; 
        end

        %i is the number of significant principal components
        i=1;
        while p_val(i)<0.05
            i=i+1;
        end
        i=i-1;

        %vecs is the principal components stored as column vectors 
        vecs=score(:,[1:i]);

end

    
    