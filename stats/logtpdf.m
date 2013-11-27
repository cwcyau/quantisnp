function [y, u] = logtpdf(x, m, iS, v)

[n, p] = size(x);

if size(m) == size(x)
  z = x - m;
else
  z = x - repmat(m, [n 1]);
end

d = sum(z'.*(iS*z'), 1);

if p == 2
  S = inv(iS);
  logdetS = log(S(1, 1)*S(2, 2) - S(1, 2)*S(2, 1));
end
if p == 1
  S = 1/iS;
  logdetS = log(S);
end

k = -0.5*logdetS;
k = k - 0.5*p*log(pi*v) + 0.5*(v+p)*log(v);
k = k + gammaln(0.5*(v+p)) - gammaln(0.5*v);

y = k - 0.5*(v+p)*log(v + d);

u = (v + p)./(v + d);
