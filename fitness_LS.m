function fit_p=fitness_LS(d_DL,X0,Y0,DL,w,tau)
N = size(Y0, 2);
D = size(w, 2);
X_s = X0(:, DL + 1  : end ) ;
for m = 1 : D
    X_s(m, :) = circshift(X0(m, DL + 1  : end ) ,  d_DL(m) );
end
% X_s=X_s(:,DL+1: end );
Q1 = N * (-1 * log(sqrt(2 * pi * tau)));
Q20 = sum((Y0 - w * X_s) .* (Y0 - w * X_s));
Q2 = (-N / (2 * tau)) * Q20;
fit_p = -(Q1 + Q2);
return;