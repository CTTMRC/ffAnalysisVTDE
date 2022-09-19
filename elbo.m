function L = elbo(phi, Pi, beta, v, W, alpha0, beta0, v0, W0)
N = size(phi, 1);
D = size(W, 1);
K = size(Pi, 1);
L0 = gammaln(sum(alpha0)) - sum(gammaln(alpha0)) - gammaln(sum(Pi)) + sum(gammaln(Pi)) - (N * D * K * log(2 * pi)) / 2;
L1=zeros(4,K);
for k = 1: K
    L1(1,k) = (-(v0 * D * log(2)) / 2) + ((v(k, :) * D * log(2)) / 2);
    L1(2,k) = -multigammaln(v0 / 2, D) + multigammaln(v(k, :) / 2, D);
    L1(3,k) = (D / 2) * log(abs(beta0)) - (D / 2) * log(abs(beta(k, :)));
    L1(4,k) = (v0 / 2) * log(det(W0)) - (v(k, :) / 2) * log(det(W(:, :, k)));
end
L1=sum(sum(L1));
aux = dot(log(phi+eps)',phi');
L2 = sum(aux);
L = L0 + L1 -L2;
end

function MG = multigammaln(a, p)
MG1 = (1 / 4) * p * (p-1) * log(pi);
MG2=zeros(1,p);
for j = 1: p
    MG2(j)= gammaln(a + (1 - j) / 2);
end
MG = MG1 + sum(MG2);
end