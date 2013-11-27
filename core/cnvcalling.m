function [ CNV, GN ] = cnvcalling(dataObj, params, hyperparams, options, chrNo)

CNV = cell(100, 1);
GN = cell(100, 1);
if isempty(dataObj)
	return;
end

fprintf('QuantiSNP. CNV Calling: ');
for chrNo = options.chrRange 

	fprintf('%2.0f ', chrNo);

	if isempty(dataObj{chrNo})
		continue;
	end

	nprobes = length(dataObj{chrNo}.pos);
	if nprobes < 10
		continue;
	end

	transitionMatrix = computeTransitionMatrix(dataObj{chrNo}.pos, hyperparams.nStates+1, options.L, hyperparams, options, chrNo);

	[obslike_sum, genotype_alloc_prob] = ...
		calcObslike(dataObj{chrNo}, params{chrNo}, hyperparams);

	[gamma, loglik(chrNo), fwdProb, backProb, fwdScale, backScale] = ...
		fwdback(hyperparams.prior, transitionMatrix, exp(obslike_sum));

	vpath = viterbi_path(hyperparams.prior, transitionMatrix, exp(obslike_sum));

	model_params.state_prior        = hyperparams.prior;
	model_params.transitionMatrix   = transitionMatrix;
	model_params.obsLike            = obslike_sum;
	model_params.fwdprob            = fwdProb;
	model_params.backprob           = backProb;
	model_params.fwdscale           = fwdScale;
	model_params.backscale          = backScale;
	model_params.gamma              = gamma;
	
	clear gamma fwdProb backProb fwdScale backScale obslike_sum transitionMatrix;

	if chrNo == options.chrX
		CNV{chrNo} = findDupDel(vpath, dataObj{chrNo}.rs, dataObj{chrNo}.pos, model_params, options.isMaleX, options.nStates+1, hyperparams);
	else
		CNV{chrNo} = findDupDel(vpath, dataObj{chrNo}.rs, dataObj{chrNo}.pos, model_params, 0, options.nStates+1, hyperparams);
	end

	if options.doGenotyping
		GN{chrNo} = genotypecalling(chrNo, vpath, genotype_alloc_prob, CNV{chrNo}, dataObj{chrNo}, hyperparams, options);
	end

	clear model_params;

end
fprintf('\n');
