%% runQCFC: 
function [QCFCVec,NaNFilter,QCFC_PropSig,QCFC_PropSigUnc,QCFC_AbsMed,QCFC_DistDep,QCFC_DistDep_Pval,QCFC_AbsMed_exc] = RunQCFC(fdJenk_m,FC,ROIDistVec)

	[QCFC,P] = GetDistCorr(fdJenk_m,FC);

	% Flatten QCFC matrix
	QCFCVec = LP_FlatMat(QCFC);
	P = LP_FlatMat(P);
	
	% Filter out NaNs:
	NaNFilter = ~isnan(QCFCVec);
	if ~any(NaNFilter)
	    error('FATAL: No data left after filtering NaNs!');
	elseif any(NaNFilter)
		fprintf(1, '\tDetected %u bad ROIs: Removed %u NaN samples from data \n', sum(sum(isnan(QCFC)) == size(QCFC,1)),sum(~NaNFilter));
	    QCFCVec = QCFCVec(NaNFilter);
	    P = P(NaNFilter);
	end

	% correct p values using FDR
	P_corrected = mafdr(P,'BHFDR','true');
	QCFC_PropSig = round(sum(P_corrected<0.05) / numel(P_corrected) * 100,2);
	QCFC_PropSigUnc = round(sum(P<0.05) / numel(P) * 100,2);

	% Find absolute median
	QCFC_AbsMed = nanmedian(abs(QCFCVec));

	% Find nodewise correlation between distance and QC-FC
	[QCFC_DistDep,QCFC_DistDep_Pval] = corr(ROIDistVec(NaNFilter),QCFCVec,'type','Spearman');

	% Optionally get QCFC AbsMed as a function of increasingly stringent exclusion (one participant at a time)
	% [fd_sorted,idx] = sort(fdJenk_m,'descend');
	% FC_sorted = FC(:,:,idx);

	% numSubs = length(fd_sorted);
	% QCFC_AbsMed_exc = [];
	% for i = 1:numSubs-20
	% 	[QCFC_temp,~] = GetDistCorr(fd_sorted(i:end),FC_sorted(:,:,i:end));
	% 	% Flatten QCFC matrix
	% 	QCFCVec_temp = LP_FlatMat(QCFC_temp);

	% 	% Find absolute median
	% 	QCFC_AbsMed_exc(i) = nanmedian(abs(QCFCVec_temp));
	% end

end
