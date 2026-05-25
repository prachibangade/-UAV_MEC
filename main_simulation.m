%% ========================================================================
%  MAIN_SIMULATION  -  Corrected two-stage UAV-MEC simulation
%
%  Implements the corrected paper:
%    - Two-stage offloading (offload ratio l, relay ratio phi)
%    - SP1 : inner BCD over {f},{tau,b},{l},{phi}   (Algorithm 1)
%    - SP2 : SCA with linearised log-rates           (Algorithm 2)
%    - Joint BCD-SCA                                 (Algorithm 3)
%
%  Produces:
%    fig2_convergence.png       - energy vs BCD iteration
%    fig3_trajectory.png        - optimised 2-D UAV trajectories
%    fig4_energy_vs_users.png   - energy vs number of users
%    fig5_energy_vs_tasksize.png- energy vs task data size
%    fig7_energy_breakdown.png  - E1..E5 component breakdown
%    fig8_relay_compare.png     - adaptive vs uniform relay (the novelty)
%
%  Requires CVX on the MATLAB path.
% ========================================================================
clc; clear; close all;

params = parameters();

fprintf('================================================================\n');
fprintf(' UAV-MEC corrected | U=%d M=%d N=%d T=%.0fs I=%.0fMbit\n', ...
        params.U, params.M, params.N, params.T, params.D(1)/1e6);
fprintf('================================================================\n\n');

%% ---- 1. Proposed scheme: BCD-SCA convergence --------------------------
fprintf('--- Proposed scheme (adaptive relay scheduling) ---\n');
[Q_prop, X_prop, hist_prop, Eb_prop] = bcd_sca_solve(params, 'adaptive', true);
E_prop = hist_prop(end);
fprintf('Proposed final energy: %.4e J\n\n', E_prop);

%% ---- 2. Baselines -----------------------------------------------------
fprintf('--- Baselines ---\n');
[E_local,  Eb_local ] = baseline_scheme('local',  params);
fprintf('  local    : %.4e J\n', E_local);
[E_full,   Eb_full  ] = baseline_scheme('full',   params);
fprintf('  full     : %.4e J\n', E_full);
[E_static, Eb_static] = baseline_scheme('static', params);
fprintf('  static   : %.4e J\n', E_static);
[E_equal,  Eb_equal ] = baseline_scheme('equal',  params);
fprintf('  equal    : %.4e J\n', E_equal);
fprintf('\n');

%% ---- 3. Key ablation: adaptive vs uniform relay scheduling ------------
fprintf('--- Ablation: uniform relay scheduling ---\n');
[E_uniform, Eb_uniform] = baseline_scheme('uniform', params);
fprintf('  uniform relay : %.4e J\n', E_uniform);
gain = 100 * (E_uniform - E_prop) / E_uniform;
fprintf('  >> adaptive relay scheduling saves %.2f %% of total energy\n', gain);
fprintf('  >> relay-energy: adaptive E4=%.3e vs uniform E4=%.3e\n\n', ...
        Eb_prop.E4, Eb_uniform.E4);

%% ---- 4. Sweep: energy vs number of users -----------------------------
fprintf('--- Sweep: energy vs U ---\n');
U_list = [4 6 8 10 12 14];
E_vs_U = zeros(5, numel(U_list));   % rows: proposed, local, full, static, uniform
for iu = 1:numel(U_list)
    pu = parameters(U_list(iu));
    [~, Xp, hp] = bcd_sca_solve(pu, 'adaptive', false);
    E_vs_U(1, iu) = hp(end);
    E_vs_U(2, iu) = baseline_scheme('local',   pu);
    E_vs_U(3, iu) = baseline_scheme('full',    pu);
    E_vs_U(4, iu) = baseline_scheme('static',  pu);
    E_vs_U(5, iu) = baseline_scheme('uniform', pu);
    fprintf('  U=%2d done\n', U_list(iu));
end
fprintf('\n');

%% ---- 5. Sweep: energy vs task data size ------------------------------
fprintf('--- Sweep: energy vs task size ---\n');
I_list = (20:15:80) * 1e6;          % bits
E_vs_I = zeros(3, numel(I_list));   % proposed, full, static
for ii = 1:numel(I_list)
    pI = parameters();
    pI.D = I_list(ii) * ones(pI.U, 1);
    [~, ~, hI] = bcd_sca_solve(pI, 'adaptive', false);
    E_vs_I(1, ii) = hI(end);
    E_vs_I(2, ii) = baseline_scheme('full',   pI);
    E_vs_I(3, ii) = baseline_scheme('static', pI);
    fprintf('  I=%2.0f Mbit done\n', I_list(ii)/1e6);
end
fprintf('\n');

%% ---- 6. Save results + plot ------------------------------------------
results.hist_prop  = hist_prop;
results.Q_prop     = Q_prop;
results.Eb_prop    = Eb_prop;
results.E_prop     = E_prop;
results.E_local    = E_local;
results.E_full     = E_full;
results.E_static   = E_static;
results.E_equal    = E_equal;
results.E_uniform  = E_uniform;
results.Eb_local   = Eb_local;
results.Eb_full    = Eb_full;
results.Eb_static  = Eb_static;
results.Eb_uniform = Eb_uniform;
results.U_list     = U_list;
results.E_vs_U     = E_vs_U;
results.I_list     = I_list;
results.E_vs_I     = E_vs_I;
results.params     = params;
save('results.mat', 'results');

plot_results(results);
fprintf('All figures saved. Done.\n');
