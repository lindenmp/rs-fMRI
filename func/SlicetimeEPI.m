function [] = SlicetimeEPI(epifile, nslices, TR, order, refslice, N)

% This function will perform standard SPM slice-timing correction.
%
% -------
% INPUTS:
% -------
% epifile     - a string containing the full path and filename of a 4D input 
%            EPI file. The images should have been co-registered to the T1. 
%            That is, they should be in T1 space; e.g.,
%            '/path/to/epi/4d.nii'
%
% nslices      - the number of slices in each volume
% 
% TR          - the repetition time of the acquisition in secs
%
% order       - the slice order of the acquisition. Some exmaples from SPM
%             help:
%
%             ascending (first slice=bottom): [1:1:nslices]
%
%             descending (first slice=top): [nslices:-1:1]
%
%             interleaved (middle-top):
%               for k = 1:nslices,
%                    round((nslices-k)/2 + (rem((nslices-k),2) * (nslices - 1)/2)) + 1;
%               end
%
%             interleaved (bottom -> up): [1:2:nslices 2:2:nslices]
%
%             interleaved (top -> down): [nslices:-2:1, nslices-1:-2:1]
%
% refslice    - reference slice usually taken in the middle of the
%             acquisition (e.g., if sequential, take round(nslices/2))
%
% N       - Lenth of time series (no of volumes in EPI 4D file).
%
% -------
% OUTPUTS:
% -------
%
% slice-time corrected images will be output with the prefix 'a'
%   
% =========================================================================
spm('defaults','fmri');
spm_jobman('initcfg')

%-----------------------------------------------------------------------
% Job configuration created by cfg_util
%-----------------------------------------------------------------------
%%


% define epi list
for i = 1:N
    
    matlabbatch{1}.spm.temporal.st.scans{1}{i,1} = [epifile,',',num2str(i)];
    
end

matlabbatch{1}.spm.temporal.st.nslices = nslices;
matlabbatch{1}.spm.temporal.st.tr = TR;
matlabbatch{1}.spm.temporal.st.ta = TR-(TR/nslices);
matlabbatch{1}.spm.temporal.st.so = order;
matlabbatch{1}.spm.temporal.st.refslice = refslice;
matlabbatch{1}.spm.temporal.st.prefix = 'a';

spm_jobman('run',matlabbatch)


