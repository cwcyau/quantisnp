function transitionMatrix = computeTransitionMatrix(pos, nStates, L, hyperparams, options, chr)

longChromosome = hyperparams.longChromosome;

if options.isMaleX == 1 & chr == options.chrX
    null_state = 2;
else
    null_state = 3;
end

% number of SNPs
nSnps = length(pos);

% compute transition probabilities using Haldane's Map function
d = diff(pos);
transitionProb = zeros(nSnps-1, nStates);
for i = 1 : nStates
	if i ~= null_state
		transitionProb(:, i) = 1-exp(-2*(d(:)/L)); %%% d (in cMs = 1 Mb)
	else
		transitionProb(:, i) = 1-exp(-2*(d(:)/longChromosome)); %%% d (in cMs = 1 Mb)
	end
end
transitionMatrix = zeros(nStates, nStates, nSnps-1);

% state change probability
for state = [1:nStates]
    transitionMatrix(state, :, :) = repmat( reshape(transitionProb(:, state)/(nStates-1), [1 1 nSnps-1]), [1 nStates 1]);
    transitionMatrix(state, state, :) = 1-transitionProb(:, state);
end
