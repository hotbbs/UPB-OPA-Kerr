% ============================================================
%  Photon number N vs drive phase phi_a for four G values.
%  Parameters: kkA = 1, Delta_a = 0, F = 0.1, U = 0.5.

%  G is swept over 0.05, 0.10, 0.15, 0.20 and phi_a over [-pi, pi].
% ============================================================
clear; clc; close all;
tic;

% ---------- System parameters ----------
kkA   = 1.0;
Delta = 0;           % dropped out of L, but kept for structural uniformity
F     = 0.1;
U     = 0.5;

M  = 501;            % number of phi_a grid points
Mn = 4;              % number of G values
N1 = 4;              % photon-number cutoff

% =====================================================================
%  Operators and precomputed Liouvillian pieces
%  (identical in every script of this project — change here, change once)
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

% ---------- Sweep axes ----------
phi_vec = pi*(2*(1:M)/M - 1);       % phi_a in [-pi, pi] — matches the pi-based ticks
G_vec   = 0.05*(1:Mn);              % G in {0.05, 0.10, 0.15, 0.20}

% ---------- Assemble the parameter-independent part of L ----------
L_fixed = L1 + U*L_Un2 + Delta*L_n; % Delta = 0 here, but written for clarity

% ---------- Output arrays ----------
g02 = zeros(Mn, M);
N_p = zeros(Mn, M);

% =====================================================================
%  Two-layer sweep:  outer = G,  inner (parallel) = phi_a
% =====================================================================
for m = 1:Mn
    G     = G_vec(m);
    L_row = L_fixed + G*L_G;        % fixed across the phi_a scan
    parfor m1 = 1:M
        ph = phi_vec(m1);
        % drive term must be rebuilt because phi_a varies
        L  = L_row + F*(exp(1i*ph)*L_ad + exp(-1i*ph)*L_a);
        rhoss = steady(L);
        rd    = full(diag(rhoss));
        N_p(m, m1) = d_n.'  * rd;
        g02(m, m1) = (d_n2.' * rd) / N_p(m, m1)^2;
    end
    fprintf('Progress: %5.1f%%\n', m/Mn*100);
end

% =====================================================================
%  Plot:  log10 of photon number vs phi_a, one curve per G
% =====================================================================
figure;
colors  = {'r','g','b','k'};
leg_str = cell(1, Mn);
for m = 1:Mn
    plot(phi_vec, log10(abs(N_p(m,:))), colors{m}, 'LineWidth', 1.2);
    hold on;
    leg_str{m} = sprintf('G = %.2f', G_vec(m));
end
hold off;
grid on;
box on;
legend(leg_str, 'FontSize', 14, 'Location', 'best');

xlabel('$\phi_a$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\lg[N]$',  'Interpreter', 'latex', 'FontSize', 20);

xticks(-pi : pi/2 : pi);
xticklabels({ '$-\pi$', '$-\frac{\pi}{2}$', '$0$', '$\frac{\pi}{2}$', '$\pi$' });
ax = gca;
ax.TickLabelInterpreter = 'latex';
xlim([-pi, pi]);

fprintf('Total elapsed: %.2f s\n', toc);