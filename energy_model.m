function [E_total, E_break] = energy_model(X, Q, params)
% ENERGY_MODEL  Total system energy for the two-stage UAV-MEC system.
%
% Implements the five energy components of the corrected paper (Eq. 33):
%   E1 : MD local computation       gamma_u * f_u^2 * (1-l) * I * O
%   E2 : MD->UAV offload transmit   P_u * sum_n tau_off[n]
%   E3 : UAV onboard computation    gamma_m * f_uav^2 * l*(1-phi) * I * O
%   E4 : UAV->TBS relay transmit    P_uav * sum_n tau_TBS[n]
%   E5 : UAV propulsion             sum_n dt*(k1 v^3 + k2/v)
%
% X must contain (produced by SP1_resource_allocation):
%   X.l       [U x 1]   offload ratio
%   X.phi     [U x 1]   relay ratio (fraction of offloaded bits sent to TBS)
%   X.assign  [U x 1]   user -> UAV assignment
%   X.f_loc   [U x 1]   MD CPU frequency
%   X.f_uav   [M x 1]   UAV CPU frequency (aggregate per UAV)
%   X.tau_off [U x N]   per-slot offload time
%   X.tau_tbs [M x N]   per-slot relay time
%   X.b       [U x N]   per-slot relay bits
%
% Outputs:
%   E_total  scalar
%   E_break  struct with fields E1..E5

M = params.M; U = params.U; N = params.N; dt = params.dt;
ku = params.kappa_user;
km = params.kappa_uav;
O  = params.cycles_per_bit;
D  = params.D(:);

% ---- E1: MD local computation -----------------------------------------
E1 = 0;
for u = 1:U
    D_loc = (1 - X.l(u)) * D(u);
    f_u   = max(X.f_loc(u), 1);
    E1 = E1 + ku * f_u^2 * O * D_loc;
end

% ---- E2: MD -> UAV offload transmission -------------------------------
%   E2 = P_u * sum_n tau_off[n]   (power fixed)
E2 = params.P_u * sum(X.tau_off(:));

% ---- E3: UAV onboard computation --------------------------------------
%   bits computed onboard UAV m = sum over its users of l*(1-phi)*I
E3 = 0;
for m = 1:M
    users_m = find(X.assign == m);
    if isempty(users_m), continue; end
    D_uav = sum( D(users_m) .* X.l(users_m) .* (1 - X.phi(users_m)) );
    f_uav = max(X.f_uav(m), 1);
    E3 = E3 + km * f_uav^2 * O * D_uav;
end

% ---- E4: UAV -> TBS relay transmission --------------------------------
%   E4 = P_uav * sum_n tau_TBS[n]
E4 = params.P_uav * sum(X.tau_tbs(:));

% ---- E5: UAV propulsion -----------------------------------------------
E5 = 0;
for m = 1:M
    for n = 2:N
        dq = squeeze(Q(m, :, n)) - squeeze(Q(m, :, n - 1));
        v  = max(norm(dq) / dt, params.v_min);
        E5 = E5 + dt * (params.k1 * v^3 + params.k2 / v);
    end
end

E_total = E1 + E2 + E3 + E4 + E5;

E_break.E1 = E1;
E_break.E2 = E2;
E_break.E3 = E3;
E_break.E4 = E4;
E_break.E5 = E5;

end
