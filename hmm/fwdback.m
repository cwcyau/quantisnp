function [gamma, loglik, eta,  alpha, beta, alphaScale, betaScale] = fwdback(prior, transmat, obslik)

[Q, T] = size(obslik);

[alpha, alphaScale, beta, betaScale, gamma, eta] = fwdhmmC(prior, transmat, obslik, Q, T);
if any(alphaScale == 0)
    loglik = -inf;
else
    loglik = sum(log(alphaScale));
end
