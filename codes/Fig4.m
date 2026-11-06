% ============================================================
%  g^(2)(0) over the (F_a, G) plane for four interaction strengths
%  U = {0.1, 1.0, 2.0, 5.0}, shown as a 2x2 tiled figure.
%
%  Parameters: kkA = 1, phi_a = pi/12, Delta_a = kkA*sin(2phi)/(2cos(2phi)).
% ============================================================
clear; clc; close all;
tic;

% ---------- System parameters ----------
kkA   = 1.0;
F_b   = 0.1;                    % (unused placeholder, kept for structure)
phi_a = pi/12;
Delta = kkA*sin(2*phi_a)/(2*cos(2*phi_a));

M    = 1*101;                   % 1*101 grid points per axis
N1   = 4;
U_All = [0.1, 1.0, 2.0, 5.0];

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

% ---------- Sweep axes (identical for every panel) ----------
F_vec = 0.152*(2*(1:M)/M - 1);      % x-axis: F_a
G_vec = 0.101*(2*(1:M)/M - 1);      % y-axis: G

% ---------- Precompute the Delta_a and phi_a block (U-independent) ----------
L_Delta  = L1 + Delta*L_n;
L_F_unit = exp(1i*phi_a)*L_ad + exp(-1i*phi_a)*L_a;   % unit-F drive

% ---------- Analytic boundary (U-independent) ----------
FA1 = 0.16*(2*(1:M)/M - 1);                              % as in the original
gg1 = 2*(FA1.^2).*exp(1i*2*phi_a)/(Delta + kkA);

% ---------- Figure ----------
figure;
TTT = tiledlayout(2, 2);
TTT.TileSpacing = 'compact';
TTT.Padding     = 'compact';

for N_u = 1:numel(U_All)
    U = U_All(N_u);

    % Fixed part for this panel: L1 + U*L_Un2 + Delta*L_n
    L_U = L_Delta + U*L_Un2;

    g02 = zeros(M, M);

    % -------- Two-layer sweep: outer = G, inner (parallel) = F_a --------
    for m = 1:M
        G     = G_vec(m);
        L_row = L_U + G*L_G;              % fixed for this row
        parfor m1 = 1:M
            L = L_row + F_vec(m1)*L_F_unit;   % one sparse add per point
            rhoss = steady(L);
            rd    = full(diag(rhoss));
            g02(m, m1) = (d_n2.' * rd) / (d_n.' * rd)^2;
        end
        if mod(m, 10) == 0
            fprintf('Panel %d (U=%.1f): %5.1f%%\n', N_u, U, m/M*100);
        end
    end

    % -------- Plot in the corresponding tile --------
    nexttile
    imagesc(F_vec, G_vec, log10(abs(g02)));
    axis xy;
    xlabel('$F_a$', 'Interpreter', 'latex');
    ylabel('$G$',    'Interpreter', 'latex');
    colormap('jet');
    cb = colorbar;
    title(cb, '$\lg[g^{(2)}(0)]$', 'FontSize', 16, 'Interpreter', 'latex');
    hold on;
    plot(FA1, real(gg1), 'w--', 'LineWidth', 1.5);
    hold off;
end

fprintf('Total elapsed: %.2f s\n', toc);