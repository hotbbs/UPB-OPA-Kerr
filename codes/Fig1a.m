% ============================================================
%  Fig. 1a — g^(2)(0) as a function of (F_a, G)
%            phi_a = pi/12 is fixed.
% ============================================================
clear; clc; close all;
tic;

% ---------- Shared parameters ----------
kkA   = 1.0;
phi_a = pi/12;
Delta = kkA*sin(2*phi_a)/(2*cos(2*phi_a));
U     = 0.5;

M  = 1*101;        % grid points per axis
N1 = 4;            % photon-number cutoff

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
FAx = 0.151*(2*(1:M)/M - 1);   % x-axis: F_a
gx  = 0.101*(2*(1:M)/M - 1);   % y-axis: G

% ---------- Assemble the parameter-independent parts of L ----------
L_fixed  = L1 + U*L_Un2 + Delta*L_n;                            % U, Delta fixed
L_F_unit = exp(1i*phi_a)*L_ad + exp(-1i*phi_a)*L_a;             % phi_a fixed

% ---------- Output ----------
g02 = zeros(M, M);

% ---------- Two-layer sweep: outer = G, inner (parallel) = F_a ----------
for m = 1:M
    G     = gx(m);
    L_row = L_fixed + G*L_G;                % fixed for this row
    parfor m1 = 1:M
        F = FAx(m1);
        L = L_row + F*L_F_unit;             % one sparse add per point
        rhoss = steady(L);
        rd    = full(diag(rhoss));
        g02(m, m1) = (d_n2.' * rd) / (d_n.' * rd)^2;
    end
    fprintf('Progress: %5.1f%%\n', m/M*100);
end

% ---------- Analytic blocking boundary ----------
FA1 = 0.152*(2*(1:M)/M - 1);
gg1 = 2*(FA1.^2).*(cos(2*phi_a) + sin(2*phi_a))/(kkA + 2*Delta);

% ---------- Plot figure ----------
figure;
imagesc(FAx, gx, log10(abs(g02)));
axis xy;
xlabel('$F_a$', 'Interpreter', 'latex');
ylabel('$G$',   'Interpreter', 'latex');
colormap('jet');
cb = colorbar;
title(cb, '$\lg[g^{(2)}(0)]$', 'FontSize', 16, 'Interpreter', 'latex');
hold on;
plot(FA1, gg1, 'w--', 'LineWidth', 2.5);
hold off;

fprintf('Fig. 1a done in %.1f s\n', toc);