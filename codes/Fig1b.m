% ============================================================
%  Fig. 1b — g^(2)(0) as a function of (G, phi_a)
%            F = 0.1 is fixed.
% ============================================================
clear; clc; close all;
tic;

% ---------- Shared parameters ----------
kkA     = 1.0;
phi_a   = pi/12;                       % only used to fix Delta here
Delta   = kkA*sin(2*phi_a)/(2*cos(2*phi_a));
U       = 0.5;
F_fixed = 0.1;

M  = 1*101;
N1 = 4;

% =====================================================================
%  Operators and precomputed Liouvillian pieces
% =====================================================================
a  = tensor(destroy(N1), identity(N1));   % mode-1 annihilation operator
ad = a';
n  = ad*a;
n2 = ad*ad*a*a;                            % n(n-1)

comm   = @(A) -1i*(spre(A) - spost(A));
L1     = 0.5*kkA*(2*spre(a)*spost(ad) - spre(n) - spost(n));  % cavity decay
L_n    = comm(n);                                            % coefficient of Delta_a
L_Un2  = comm(n2);                                           % coefficient of U
L_G    = comm(1i*(ad*ad - a*a));                             % coefficient of G
L_ad   = comm(ad);                                           % coefficient of F*exp( i*phi_a)
L_a    = comm(a);                                            % coefficient of F*exp(-i*phi_a)

d_n  = full(diag(n));
d_n2 = full(diag(n2));

% ---------- Axes ----------
G_vec   = 0.101*(2*(1:M)/M - 1);   % x-axis: G
phi_vec = 3.14 *(2*(1:M)/M - 1);   % y-axis: phi_a

% ---------- Assemble the parameter-independent parts of L ----------
L_fixed = L1 + U*L_Un2 + Delta*L_n;                             % U, Delta fixed

% ---------- Output ----------
g02 = zeros(M, M);

% ---------- Two-layer sweep: outer = phi_a, inner (parallel) = G ----------
for m = 1:M
    ph    = phi_vec(m);
    L_row = L_fixed + F_fixed*(exp(1i*ph)*L_ad + exp(-1i*ph)*L_a);  % fixed row
    parfor m1 = 1:M
        G = G_vec(m1);
        L = L_row + G*L_G;                 % one sparse add per point
        rhoss = steady(L);
        rd    = full(diag(rhoss));
        g02(m, m1) = (d_n2.' * rd) / (d_n.' * rd)^2;
    end
    fprintf('Progress: %5.1f%%\n', m/M*100);
end

% ---------- Analytic boundary: G = 2F^2(cos2phi + sin2phi)/(kkA + 2*Delta) ----------
G_ana = 2*(F_fixed^2)*(cos(2*phi_vec) + sin(2*phi_vec))/(kkA + 2*Delta);

% ---------- Plot figure ----------
figure;
imagesc(G_vec, phi_vec, log10(abs(g02)));
axis xy;
xlabel('$G$',       'Interpreter', 'latex');
ylabel('$\phi_a$',  'Interpreter', 'latex');
colormap('jet');
cb = colorbar;
title(cb, '$\lg[g^{(2)}(0)]$', 'FontSize', 16, 'Interpreter', 'latex');
hold on;
plot(G_ana, phi_vec, 'w--', 'LineWidth', 2.5);
hold off;

fprintf('Fig. 1b done in %.1f s\n', toc);