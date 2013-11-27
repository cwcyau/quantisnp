% QuantiSNP v1.0
% --------------
% An Objective Bayesian Hidden Markov Model approach for detecting copy
% number variation from SNP genotyping data.
%
% Copyright (C) 2007  The University of Oxford.
%
% For all enquiries please consult the contacts page on the program website: 
%   http://www.well.ox.ac.uk/quantisnp 
%
function [path, loglik] = viterbi_path(prior, transmat, obslik)
%
% Viterbi Algorithm.
%
%
% Written on 19/01/2007
%
% Christopher Yau
% Department of Statistics  
% University of Oxford
%
% This is modified from Kevin Murphy's viterbi_path.m script in the HMM Toolbox.
% http://www.cs.ubc.ca/~murphyk/Software/HMM/hmm.html
%

% number of hidden states and datapoints
[Q, T] = size(obslik);

% prior distribution 
prior = prior(:);

% size of transition matrix
transmatsize = length(size(transmat));

% predefine matrices for storing traceback, path and log likelihood matrices
delta = zeros(Q, T);
psi   = zeros(Q, T);
path  = zeros(1, T);
scale = ones(1, T);

% calculate for first datapoint
delta(:, 1)      = prior .* obslik(:, 1);
[delta(:, 1), n] = normaliseC(delta(:, 1));
scale(1)         = n;
psi(:, 1)        = 0; 

% for each subsequent data point
for t = 2 : T
    
    % for each state
    for j = 1 : Q
        
        % choose homogeneous or heterogeneous transition matrix
        if transmatsize == 2
            [delta(j, t), psi(j, t)] = max(delta(:, t-1) .* transmat(:, j));
        else
            tmat(:, :) = transmat(:, j, t-1);
            [delta(j, t), psi(j, t)] = max(delta(:, t-1) .* tmat);
        end
        
        delta(j, t) = delta(j, t) * obslik(j, t);
    end
    
    [delta(:, t), n] = normaliseC(delta(:, t));
    scale(t) = n;   
    
end

% traceback
[p, path(T)] = max(delta(:, T));
for t = T-1 : -1 : 1
    path(t) = psi(path(t+1), t+1);
end

% compute log likelihood
loglik = sum(log(scale));

