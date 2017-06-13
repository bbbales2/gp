clear all; clc;

sigma = 0.25;
l = 0.20;

a = 1 / (4 * sigma^2);
b = 1 / (2 * l^2);
c = sqrt(a^2 + 2 * a * b);

epsilon = sqrt(b);
alpha = sqrt(2 * a);
beta = (1 + (2 * epsilon / alpha)^2)^0.25;
delta = sqrt(alpha^2 * (beta^2 - 1) / 2);

M = 10;
N = 51;
Q = zeros(N, M);
Q2 = zeros(N, M);
V = zeros(M, 1);
x = linspace(-0.5, 0.5, N);%4 * (rand(N, 1) - 0.5)';

tic
figure(1);
rho = alpha * exp(-alpha^2 * x.^2) / sqrt(pi);
H = zeros(N, M);
xp = alpha * beta * x;
H(:, 1) = sqrt(beta) .* exp(-delta^2 * x .* x)';%%sqrt(alpha * exp(-alpha^2 * x .* x) / sqrt(pi))' .* 
H(:, 2) = sqrt(1 / 2) * 2 * xp' .* H(:, 1);
for n = 3:M
   H(:, n) = sqrt(1 / (2 * (n - 1))) * 2 * xp' .* H(:, n - 1) - sqrt(1 / (4 * (n - 1) * (n - 2))) * 2 * (n - 2) * H(:, n - 2); 
end

Ht = zeros(N, M);
xp = alpha * beta * x;
f = sqrt(epsilon^2 / (alpha^2 + delta^2 + epsilon^2));
Ht(:, 1) = sqrt(sqrt(alpha^2 / (alpha^2 + delta^2 + epsilon^2))) * sqrt(beta) * exp(-delta^2 * x .* x)';
Ht(:, 2) = f * sqrt(1 / 2) * 2 * xp' .* Ht(:, 1);
for n = 3:M
   Ht(:, n) = f * sqrt(1 / (2 * (n - 1))) * 2 * xp' .* Ht(:, n - 1) - f^2 * sqrt(1 / (4 * (n - 1) * (n - 2))) * 2 * (n - 2) * Ht(:, n - 2); 
end

H3 = zeros(N, M);
xp = alpha * beta * x;
H3(:, 1) = exp(-x .* x / 2.0) / sqrt(sqrt(pi));
H3(:, 2) = 2 * xp' .* H3(:, 1) / sqrt(2);
for n = 3:M
   H3(:, n) = (2 * xp' .* H3(:, n - 1) - 2 * (n - 2) * H3(:, n - 2)) / sqrt((n - 1) * 2); 
end

for n = 1:M
    %%gam = sqrt(beta / (2^(n - 1) * gamma(n)));
    V(n) = sqrt(alpha^2 / (alpha^2 + delta^2 + epsilon^2)) * (epsilon^2 / (alpha^2 + delta^2 + epsilon^2))^(n - 1);

    %%Q2(:, n) = gam * H3(:, n) .* exp(-delta^2 * x .* x)';
    Q2(:, n) = H(:, n);
    %%gam * Ht(:, n) - H(:, n)
    Q(:, n) = H(:, n);
    plot(x, Q(:, n) * V(n));
    hold on;

end
hold off;
toc

T = Q * diag(sqrt(V));

Kt = T * T';
Ktt = Ht * Ht';

%Kt = Q * diag(V) * (Q');

K = zeros(N, N);
for n = 1:N
    for m = 1:N
        K(n, m) = exp(-epsilon^2 * (x(n) - x(m))^2);
    end
end

figure(2);
subplot(2, 2, 1);
imshow(Kt);

subplot(2, 2, 2);
imshow(K);

subplot(2, 2, 3);
imshow(log(abs(K - Ktt) + 1e-16), []);
colorbar();
title('Error in approx');

subplot(2, 2, 4);
imshow(log(abs(K - Kt) + 1e-16), []);
colorbar();
title('Error in approx');

Q3 = Q2 * diag(sqrt(V)) * Q2';
max(max(Q3 * Q3' / 2.0 - Kt));
max(max(Q2 * diag(V) * Q2' - Kt));