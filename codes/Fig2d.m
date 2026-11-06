% ============================================================
%  Photon number N vs Delta_a for four U values.
%  Parameters: kkA = 1, G = 0.02, F = 0.1, phi_a = pi/8.
%
%  Here G, F, phi_a are fixed; U is swept over {0.5, 1.0, 1.5, 2.0}
%  and Delta_a over [-3.5, 3.5].
% ============================================================
clear; clc; close all;
tic;

% ---------- System parameters ----------
kkA   = 1.0;
G     = 0.02;
F     = 0.1;
phi_a = pi/8;

M  = 501;         % number of Delta_a grid points
Mn = 4;           % number of U values
N1 = 4;           % photon-number cutoff

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

% ---------- Sweep axes ----------
Delta_vec = 3.5*(2*(1:M)/M - 1);   % Delta_a in [-3.5, 3.5]
U_vec     = 0.5*(1:Mn);            % U in {0.5, 1.0, 1.5, 2.0}

% ---------- Assemble the U- and Delta-independent part of L ----------
L_fixed = L1 + G*L_G + F*(exp(1i*phi_a)*L_ad + exp(-1i*phi_a)*L_a);

% ---------- Output arrays ----------
g02 = zeros(Mn, M);
N_p = zeros(Mn, M);

% =====================================================================
%  Two-layer sweep:  outer = U,  inner (parallel) = Delta_a
% =====================================================================
for m = 1:Mn
    U     = U_vec(m);
    L_row = L_fixed + U*L_Un2;              % one sparse add per row
    parfor m1 = 1:M
        L = L_row + Delta_vec(m1)*L_n;      % one sparse add per point
        rhoss = steady(L);
        rd    = full(diag(rhoss));
        N_p(m, m1) = d_n.'  * rd;
        g02(m, m1) = (d_n2.' * rd) / N_p(m, m1)^2;
    end
    fprintf('Progress: %5.1f%%\n', m/Mn*100);
end

% =====================================================================
%  Plot:  log10 of photon number vs Delta_a, one curve per U
% =====================================================================
figure;
colors  = {'r','g','b','k'};
leg_str = cell(1, Mn);
for m = 1:Mn
    plot(Delta_vec, log10(abs(N_p(m,:))), colors{m}, 'LineWidth', 1.2);
    hold on;
    leg_str{m} = sprintf('U = %.1f', U_vec(m));
end
hold off;

legend(leg_str, 'FontSize', 14, 'Location', 'best');

xlabel('$\Delta_a$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\lg[N]$',   'Interpreter', 'latex', 'FontSize', 20);

fprintf('Total elapsed: %.2f s\n', toc);