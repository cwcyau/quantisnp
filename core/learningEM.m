function params = learningEM(dataObj, params, hyperparams, options)

for it = 1 : options.EMiters

	if options.doVerbose
		disp(['QuantiSNP. EM Iteration: ' num2str(it) ]);
	end
	
    %% calculate expectations
	if options.doVerbose
		fprintf('\tCalculating expectations: ');
	end
    loglik = zeros(1, max(options.chrRange));
    for chrNo = options.chrRange

		if options.doVerbose
			fprintf(' %g', chrNo);
		end

		expectationObj{chrNo} = [];

		if isempty(dataObj{chrNo})
			continue;
		end

        nprobes = length(dataObj{chrNo}.pos);
        if nprobes < 10
            continue;
        end

        transitionMatrix = computeTransitionMatrix(dataObj{chrNo}.pos, hyperparams.nStates+1, options.L, hyperparams, options, chrNo);

		[obslike_sum, genotype_alloc_prob, comp_alloc_prob, u] = ...
			calcObslike(dataObj{chrNo}, params, hyperparams);

        [gamma, loglik(chrNo)] = ...
            fwdback(hyperparams.prior, transitionMatrix, exp(obslike_sum));

        expectationObj{chrNo}.gamma 				= gamma;
        expectationObj{chrNo}.genotype_alloc_prob 	= genotype_alloc_prob;
        expectationObj{chrNo}.comp_alloc_prob 		= comp_alloc_prob;
        expectationObj{chrNo}.u 					= u;

		clear u gamma genotype_alloc_prob comp_alloc_prob

    end
	if options.doVerbose
		fprintf('\n');
	end

	if options.doVerbose
		dispstr('\tLog-likelihood: %s', sum(loglik));
	end

    %% update parameters
    if options.doVerbose
		dispstr('\tUpdating parameters.');
	end

    params = updateParams(expectationObj, dataObj, params, hyperparams, options);

	if options.doVerbose

		fprintf('\tOutlier Rate: %s\n', num2str(params.nu));
		
		fprintf('\tMean LRR Levels: ');
		for k = 1 : options.nStates
			fprintf('%1.3f ', params.m_r(k));
		end
		fprintf('\n');
		
		fprintf('\tMixture Component Parameters - \n');
		for j = 1 : 3

			fprintf('\t\t%g: Weights = ', j);
			for k = 1 : options.nComp
				fprintf('%s ', num2str( params.q{j}(k), '%1.3f' ));
			end
			fprintf('\n');
			
			fprintf('\t\t%g: Centres (LRR) = ', j);
			for k = 1 : options.nComp
				fprintf('%s ', num2str( params.d_r{j}(k), '%1.3f' ));
			end
			fprintf('\n');

			fprintf('\t\t%g: Centres (BAF) = ', j);
			for k = 1 : options.nComp
				fprintf('%s ', num2str( params.d_b{j}(k), '%1.3f' ));
			end
			fprintf('\n');

			fprintf('\t\t%g: Std. Dev. (LRR) = ', j);
			for k = 1 : options.nComp
				fprintf('%s ', num2str( chol(params.S{j}(2, 2, k)), '%1.3f' ));
			end		
			fprintf('\n');

			fprintf('\t\t%g: Std. Dev. (BAF) = ', j);
			for k = 1 : options.nComp
				fprintf('%s ', num2str( chol(params.S{j}(1, 1, k)), '%1.3f' ));
			end
			fprintf('\n');
			
		end
			
		fprintf('\tStd. Dev. (CNV) = %s\n', num2str(sqrt(params.S_r_cnv), '%1.3f') );
		fprintf('\tStd. Dev. (HomDel) = %s\n', num2str(sqrt(params.S_r_homdel), '%1.3f') );

		fprintf('\tGenotype Prob:\n');
		for i = 1 : options.nStates
			fprintf('\t\t%s: ', num2str(i) );
			for j = 1 : length(params.w{i})
				fprintf('%s ', num2str(params.w{i}(j), '%1.3f') );
			end
			fprintf('\n');
		end
		fprintf('\n');
		
		end
		
	end
	
end
