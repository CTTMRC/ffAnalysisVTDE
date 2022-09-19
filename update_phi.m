function phi = update_phi(pI, m, v, W, beta, x, K, D)
phi=zeros(size(x,1),K);
for k = 1: K
    phi0 = psi(pI(k)+eps) - psi(sum(pI)+eps);
    phi1 = m(k, :) * v(k, :) / W(:, :, k)* x';
    phi2 = sum((1 / 2) * v(k, :)* pinv( W(:, :, k))* x' .* x', 1);
    phi3 = (D / 2) * (1/ beta(k, :));
    phi4 = (1 / 2) * v(k, :) * m(k, :) / (W(:, :, k)) *  m(k, :)';
    phi5 = (D / 2) * log(2);
    phi60=ones(1,D);
    for i = 1: D
        phi60(i) = psi((v(k, :) + 1 - i) / 2);
    end
    phi6 = (1 / 2) * sum(phi60);
    phi7 = (1 / 2) * log(abs(det(W(:, :, k))));
    phi(:, k) = phi0 + phi1' - phi2' - phi3 - phi4 + phi5 + phi6 - phi7;
end
phi = softmax(phi);
function sm = softmax(x)
ex = exp(x);
sumEx=sum(ex, 2)+eps;
sm = ex ./ sumEx;