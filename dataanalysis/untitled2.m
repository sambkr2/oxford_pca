%% Find eigendecomp of linear operator for noisy matrix
mylambdas = [];

Y = velo;

% build data matrices
Y1 = Y(:, 1:end-1);
Y2 = Y(:, 2:end);

% compute DMD
[U, S, V] = svd(Y1, 'econ');
Atilde = U' * Y2 * V / S;

[W, D] = eig(Atilde);
lam = diag(D); % eigenvalues

mylambdas = [mylambdas; lam];

%% now forward/backward dmd

fblambdas = [];

% compute forward DMD
[U, S, V] = svd(Y1, 'econ');
f_Atilde = U' * Y2 * V / S;

% compute backward DMD
[U, S, V] = svd(Y2, 'econ');
b_Atilde = U' * Y1 * V / S;

% estimate forward/backward DMD
Atilde = (f_Atilde * inv(b_Atilde)) ^ 0.5;
[W, D] = eig(Atilde);
lam = diag(D);

fblambdas = [fblambdas; lam];

%% Plot eigenvalues
% neigh = sig * [-1 0.5];
gray = 0.7 * [1 1 1];

figure;
hold on;
plot(real(mylambdas), imag(mylambdas), '.', 'Color', 'b');
plot(real(fblambdas), imag(fblambdas), '+', 'Color', gray);
% plot(real(lambda), imag(lambda), 'r^', 'MarkerFaceColor', 'r');
rectangle('Position', [-1 -1 2 2], 'Curvature', 1, ...
    'LineStyle', '--');
axis square;
l1 = legend('DMD', 'fbDMD', 'true');
set(l1,'FontSize',14)
% xlim(real(lambda(1)) + neigh);
% ylim(imag(lambda(1)) + neigh);
% title(sprintf('Noise = %g', sig));