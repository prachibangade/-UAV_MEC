function params = parameters(U_override)
% PARAMETERS  Configuration for the corrected two-stage UAV-MEC framework.
%
%   Aligned with the corrected paper:
%     - Slotted timing model: T = N*dt is the deadline.
%     - Fixed MD transmit power P_u and fixed bandwidth B.
%     - Two-stage offloading: offload ratio l, relay ratio phi.
%     - b[n] (per-slot relay bits) is a true optimisation variable in SP1.
%
%   PARAMETER REGIME (important for the trajectory trade-off):
%     The horizon T is kept short (40 s) and the task load high (70 Mbit)
%     so that the UAV is time-pressured. Under a tight deadline the UAV
%     must fly toward users / the TBS to raise link rates and meet the
%     delivery constraints -- this is the coupling that makes trajectory
%     optimisation non-trivial (same mechanism as Ji 2021 [1], Zeng 2019
%     [16]). The propulsion constants k1, k2 are chosen for a light
%     rotary-wing UAV so that propulsion energy E5 is the SAME ORDER as
%     communication+relay energy E2+E4 (ratio ~0.8). This produces a
%     visible, defensible energy trade-off rather than propulsion
%     dominating everything.

if nargin < 1 || isempty(U_override)
    params.U = 10;
else
    params.U = U_override;
end

% ---------------- Geometry --------------------------------------------
rng(42);
params.area_side = 500;                              % m
params.user_pos  = params.area_side * (0.12 + 0.76 * rand(params.U, 2));

% TBS off-centre so UAV-to-TBS distance varies along the trajectory.
params.GS = [60, 250];                               % m

params.M = 2;
params.q_start = [ 120, 120;
                   120, 380];
params.q_end   = [ 420, 120;
                   420, 380];

% ---------------- Time discretisation ---------------------------------
params.N  = 40;                                      % time slots
params.T  = 60;                                      % s  (horizon = deadline)
params.dt = params.T / params.N;
params.H  = 100;                                     % m altitude

% ---------------- Wireless --------------------------------------------
params.B     = 2e6;                                  % Hz
params.noise = 1e-9;                                 % W
params.beta0 = 1e-3;                                 % reference channel gain @1 m

params.P_u   = 0.5;                                  % MD transmit power (FIXED)
params.P_uav = 2.0;                                  % UAV relay transmit power (FIXED)

% ---------------- Computation -----------------------------------------
params.f_max_user = 1.0e9;                           % 1 GHz
params.f_max_uav  = 5e9;                             % 5 GHz
params.kappa_user = 1e-28;
params.kappa_uav  = 1e-28;
params.cycles_per_bit = 1000;                        % O_u

% Task data size per MD (bits).
params.D = 50e6 * ones(params.U, 1);

% ---------------- UAV propulsion (light rotary-wing regime) -----------
params.k1    = 0.0020;                               % blade profile (kappa1)
params.k2    = 3.0;                                  % induced power  (kappa2)
params.v_max = 25;                                   % m/s
params.v_min = 1;                                    % m/s numerical floor

% ---------------- Collision avoidance ---------------------------------
params.d_min    = 20;                                % m
params.d_min_sq = params.d_min^2;

% ---------------- Optimisation controls -------------------------------
params.max_iter = 20;                                % outer BCD iterations
params.sca_iter = 8;                                 % inner SCA iterations (SP2)
params.sp1_iter = 6;                                 % inner BCD iterations (SP1)
params.conv_tol = 1e-4;

params.l_min   = 0.01;
params.l_max   = 0.99;
params.phi_min = 0.01;
params.phi_max = 0.99;

% ---------------- Relay scheduling mode -------------------------------
params.relay_mode = 'adaptive';

end
