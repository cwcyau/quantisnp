function CNV = findDupDel(stateSeq, rs, pos, model_params, isMaleX, nStates, hyperparams)

% get HMM probabilities
state_prior        = model_params.state_prior;
transitionMatrix   = model_params.transitionMatrix;
obsLike            = model_params.obsLike;
fwdprob            = model_params.fwdprob;
backprob           = model_params.backprob;
fwdscale           = model_params.fwdscale;
backscale          = model_params.backscale;
gamma              = model_params.gamma;

copy_number = hyperparams.copy_number;
cnloh = hyperparams.cnloh;

% find number of states and number of SNP markers
[nStates, nSnps] = size(obsLike);


% is this an X chromosome on a male? if so, then the null state is copy number 1
% (state 2), otherwise null state is copy number 2 (states 3).
if isMaleX == 1
    null_state = 2;
else
    null_state = 3;
end

%
% find each deletion or duplication
%

count = 1; % variable to count number of CNVs
start = 0; % variable to store start SNP of CNV 

% for each marker
CNV = [];
for i = 1 : nSnps	

    % is this SNP normal?
	isNormal = length( find(null_state == stateSeq(i)) );

    % if we detect a switch from the normal state then start a new CNV
    % region
	if isNormal == 0 & start == 0 & i ~= nSnps
		start = 1;        
		CNV{count}.cnloh = cnloh(stateSeq(i));
		CNV{count}.type = stateSeq(i);
        CNV{count}.copy = copy_number(stateSeq(i)); 
        CNV{count}.index(1) = i;
        CNV{count}.location(1) = pos(i);
		CNV{count}.rs(1) = rs(i);		
		continue;
	end 

    % if we are already in a CNV region
	if i > 1 & start == 1
        
        % if we switch back from an aberrant state to a normal state then record
        % CNV information
		if stateSeq(i) ~= stateSeq(i-1) & isNormal == 1 
			CNV{count}.index(2) = i-1;
        	CNV{count}.location(2) = pos(i-1);	
			CNV{count}.rs(2) = rs(i-1);

			change = CNV{count}.type;
        	t(1) = CNV{count}.index(1);
        	t(2) = CNV{count}.index(2);			
			s(1) = CNV{count}.type;
        	s(2) = stateSeq(t(2));

			bf = computeBayesFactor(change, s, t, null_state, model_params, nStates);
			CNV{count}.delta = bf;		

			start = 0;
			count = count + 1;
		end

        % if we move from one aberrant state to another, then store
        % information for the original CNV and start recording for the new
        % one
		if stateSeq(i) ~= stateSeq(i-1) & isNormal == 0  
            
            if isMaleX == 1 & stateSeq(i) == 3 & stateSeq(i-1) == 3 
                continue;    
            end
            
			CNV{count}.index(2) = i-1;
        	CNV{count}.location(2) = pos(i-1);
			CNV{count}.rs(2) = rs(i-1);

			change = CNV{count}.type;
        	t(1) = CNV{count}.index(1);
        	t(2) = CNV{count}.index(2); 			
			s(1) = CNV{count}.type;
        	s(2) = stateSeq(t(2));

			bf = computeBayesFactor(change, s, t, null_state, model_params, nStates);
			CNV{count}.delta = bf;

			count = count + 1;

			CNV{count}.cnloh = cnloh(stateSeq(i));
			CNV{count}.type = stateSeq(i);
        	CNV{count}.copy = copy_number(stateSeq(i)); 
        	CNV{count}.index(1) = i;
        	CNV{count}.location(1) = pos(i);
			CNV{count}.rs(1) = rs(i);

		end
	end

    % if we reach the end of the region and we are still in an aberrant
    % region then store CNV information
	if i == nSnps & start == 1
        
		CNV{count}.index(2) = i;
        CNV{count}.location(2) = pos(i);
		CNV{count}.rs(2) = rs(i);
		
		change = CNV{count}.type;
        t(1) = CNV{count}.index(1);
        t(2) = CNV{count}.index(2); 			
		s(1) = CNV{count}.type;
        s(2) = stateSeq(t(2));

		bf = computeBayesFactor(change, s, t, null_state, model_params, nStates);
		CNV{count}.delta = bf;

	end

end 



function bf = computeBayesFactor(change, s, t, null_state, model_params, nStates)
% 
% compute Bayes Factor
%

% get HMM probabilities
state_prior        = model_params.state_prior;
transitionMatrix   = model_params.transitionMatrix;
obsLike            = model_params.obsLike;
alpha              = model_params.fwdprob;
beta               = model_params.backprob;
alphaScale         = model_params.fwdscale;
betaScale          = model_params.backscale;
gamma              = model_params.gamma;

% get number of normal states
nNorm = length(null_state);  

% compute P(data|all paths that go through the normal state)
idx = [ t(1) max(t(1), t(2)-1) ];

% make sure we don't go off the end!
n = size(transitionMatrix, 3);
if ( ( idx(2) > n ) & ( idx(1) > n ) )
	idx(2) = n;
	idx(1) = n;
end

% compute P(data|all paths that go through the non-normal state)
for newState = 1 : nStates

	tMat = 0*transitionMatrix(:, :, idx(1):idx(2));
	tMat(newState, newState, :) = 1;
	transitionMatrix(:, :, idx(1):idx(2)) = tMat;

	% compute probability of data up to CNV
	fwdprob = sum(log(alphaScale(1:idx(1)))); 

	% compute backward probability of data up to CNV
	bkprob  = sum(log(betaScale(idx(2)+1:end))) + log(sum(beta(newState, idx(2)+1))) ;

 	% compute probability of data given that the CNV is aberrant
 	pr = tMat(:, :, 1)*alpha(:, idx(1));
 	[gammaS, segprob] = fwdback(pr, tMat(:, :, 2:end), exp(obsLike(:, idx(1)+1:idx(2)+1)));
 	bf(newState) = fwdprob + bkprob + segprob;
	
end

bf = bf - bf(null_state);

