function X = SP1_resource_allocation(Q, params, X_init)
% SP1_RESOURCE_ALLOCATION  Resource allocation for fixed trajectory (Algorithm 1).
%
% Solves SP1 by an INNER BCD over four sub-blocks of X:
%   Block 1: {f_loc, f_uav}            -> convex (quadratic energy, c/f delay)
%   Block 2: {tau_off, tau_tbs, b}     -> linear program
%   Block 3: {l}                       -> convex (with phi, b, tau fixed)
%   Block 4: {phi}                     -> convex (with l, b, tau fixed)
%
% Each block is a convex CVX program; the bilinear term l*phi in
% constraint (P1c) becomes affine once the partner variable is fixed.
%
% Inputs:
%   Q       - fixed UAV trajectory [M x 2 x N]
%   params  - parameter struct
%   X_init  - (optional) warm-start X struct
%
% Output:
%   X - resource allocation struct (see energy_model.m for fields)

M = params.M; U = params.U; N = params.N; dt = params.dt;
O  = params.cycles_per_bit;
D  = params.D(:);
ku = params.kappa_user;
km = params.kappa_uav;

% ---- user -> UAV assignment (nearest mean position) -------------------
Q_mean = squeeze(mean(Q, 3));            % [M x 2]
assign = zeros(U, 1);
for u = 1:U
    d2 = sum((Q_mean - params.user_pos(u, :)).^2, 2);
    [~, assign(u)] = min(d2);
end

% ---- constant rates for this fixed trajectory -------------------------
[~, Ruser] = channel_model(Q, params, params.P_u * ones(U, 1));  % [M x U x N]
Rrelay = relay_rate(Q, params);                                  % [M x N]

% per-(u,n) uplink rate for the assigned UAV
Ru = zeros(U, N);
for u = 1:U
    Ru(u, :) = squeeze(Ruser(assign(u), u, :))';
end
Ru = max(Ru, 1e2);

% ---- initial feasible X -----------------------------------------------
if nargin >= 3 && ~isempty(X_init)
    X = X_init;
    X.assign = assign;        % assignment may change with Q
else
    X = init_feasible_X(assign, Ru, Rrelay, params);
end
X.Rrelay_cache = Rrelay;      % used by Blocks 3 and 4

% ====================== INNER BCD ======================================
for it = 1:params.sp1_iter

    % ---- Block 1: CPU frequencies {f_loc, f_uav} ----------------------
    X = block_frequencies(X, params, O, D, ku, km);

    % ---- Block 2: {tau_off, tau_tbs, b} as a linear program -----------
    X = block_times_bits(X, params, D, Ru, Rrelay);

    % ---- Block 3: offload ratio {l} -----------------------------------
    X = block_offload_ratio(X, params, O, D, ku, km, Ru);

    % ---- Block 4: relay ratio {phi} -----------------------------------
    X = block_relay_ratio(X, params, O, D, km);

end

X.assign = assign;
X.Ru     = Ru;
X.Rrelay = Rrelay;

end


% ======================================================================
%  Initial feasible point
% ======================================================================
function X = init_feasible_X(assign, Ru, Rrelay, params)
M = params.M; U = params.U; N = params.N; dt = params.dt;
D = params.D(:);

X.assign = assign;
X.l   = 0.5 * ones(U, 1);
X.phi = 0.5 * ones(U, 1);

% offload time: split required bits evenly across slots
X.tau_off = zeros(U, N);
for u = 1:U
    need = X.l(u) * D(u);                  % bits to offload
    % even split, capped by per-slot budget
    t_each = need / sum(Ru(u, :));
    X.tau_off(u, :) = min(t_each, 0.45 * dt);
end

% relay time and bits: even split
X.tau_tbs = zeros(M, N);
X.b       = zeros(U, N);
for u = 1:U
    m = assign(u);
    relay_bits = X.l(u) * X.phi(u) * D(u);
    X.b(u, :)  = relay_bits / N;
end
for m = 1:M
    users_m = find(assign == m);
    for n = 1:N
        bsum = sum(X.b(users_m, n));
        X.tau_tbs(m, n) = min(bsum / Rrelay(m, n), 0.45 * dt);
    end
end

% frequencies: enough to meet deadline
X.f_loc = zeros(U, 1);
for u = 1:U
    C = params.cycles_per_bit * (1 - X.l(u)) * D(u);
    X.f_loc(u) = min(max(C / params.T, 1), params.f_max_user);
end
X.f_uav = zeros(M, 1);
for m = 1:M
    users_m = find(assign == m);
    C = params.cycles_per_bit * sum(D(users_m) .* X.l(users_m) .* (1 - X.phi(users_m)));
    X.f_uav(m) = min(max(C / params.T, 1), params.f_max_uav);
end

end


% ======================================================================
%  Block 1 : CPU frequencies  (convex QP, closed form)
% ======================================================================
function X = block_frequencies(X, params, O, D, ku, km)
% Energy gamma*f^2*C is increasing in f; the delay constraint C/f <= T
% forces f >= C/T. The energy-minimal choice is therefore f = C/T exactly
% (clipped to [1, f_max]). This is the convex-subproblem optimum.
U = params.U; M = params.M;

for u = 1:U
    C = O * (1 - X.l(u)) * D(u);
    X.f_loc(u) = min(max(C / params.T, 1), params.f_max_user);
end
for m = 1:M
    users_m = find(X.assign == m);
    C = O * sum(D(users_m) .* X.l(users_m) .* (1 - X.phi(users_m)));
    X.f_uav(m) = min(max(C / params.T, 1), params.f_max_uav);
end
end


% ======================================================================
%  Block 2 : {tau_off, tau_tbs, b}  -- LINEAR PROGRAM (CVX)
% ======================================================================
function X = block_times_bits(X, params, D, Ru, Rrelay)
M = params.M; U = params.U; N = params.N; dt = params.dt;
l = X.l; phi = X.phi; assign = X.assign;

cvx_clear
cvx_begin quiet
    variable tau_off(U, N) nonnegative
    variable tau_tbs(M, N) nonnegative
    variable b(U, N) nonnegative

    % objective: minimise transmission energy E2 + E4 (power fixed)
    minimize( params.P_u * sum(sum(tau_off)) + params.P_uav * sum(sum(tau_tbs)) )

    subject to
        % (P1b) offloading delivers l*I bits
        for u = 1:U
            sum( tau_off(u, :) .* Ru(u, :) ) >= l(u) * D(u);
        end
        % (P1c) all relay bits forwarded:  sum_n b = l*phi*I
        for u = 1:U
            sum( b(u, :) ) == l(u) * phi(u) * D(u);
        end
        % (P1d) per-slot relay capacity
        for m = 1:M
            users_m = find(assign == m);
            for n = 1:N
                sum( b(users_m, n) ) <= tau_tbs(m, n) * Rrelay(m, n);
            end
        end
       
        % (P1g) per-slot time budget
        for m = 1:M
            users_m = find(assign == m);
            for n = 1:N
                sum( tau_off(users_m, n) ) + tau_tbs(m, n) <= dt;
            end
        end
cvx_end

if strcmp(cvx_status,'Solved') || contains(cvx_status,'Inaccurate')
    X.tau_off = max(tau_off, 0);
    X.tau_tbs = max(tau_tbs, 0);
    X.b       = max(b, 0);
else
    warning('SP1 Block-2 LP status: %s -- keeping previous {tau,b}.', cvx_status);
end
end


% ======================================================================
%  Block 3 : offload ratio {l}  (convex CVX program)
% ======================================================================
function X = block_offload_ratio(X, params, O, D, ku, km, Ru)
% With phi fixed, the bilinear product l*phi in (P1c) becomes AFFINE in l.
% We jointly update {l, b} so that (P1c) stays satisfied exactly:
%   minimise   E1(l) + E3(l) + E4-proxy
%   E1 = ku*f_loc^2*O*(1-l)*D     (affine, decreasing in l)
%   E3 = km*f_uav^2*O*l*(1-phi)*D (affine, increasing in l)
% subject to (P1b) offload feasibility and (P1c) with b re-scaled.
% Frequencies are held at the Block-1 values (recomputed next sweep).
U = params.U; N = params.N;
phi = X.phi;  f_loc = X.f_loc;  f_uav_of = zeros(U,1);
for u = 1:U, f_uav_of(u) = X.f_uav(X.assign(u)); end

cvx_clear
cvx_begin quiet
    variable l(U)
    variable b(U, N) nonnegative

    e1 = 0; e3 = 0;
    for u = 1:U
        e1 = e1 + ku * f_loc(u)^2 * O * (1 - l(u)) * D(u);
        e3 = e3 + km * f_uav_of(u)^2 * O * l(u) * (1 - phi(u)) * D(u);
    end
    minimize( e1 + e3 )

    subject to
        l >= params.l_min;  l <= params.l_max;
        for u = 1:U
            % (P1b) offload feasibility (tau fixed)
            sum( X.tau_off(u,:) .* Ru(u,:) ) >= l(u) * D(u);
            % (P1c) total relay bits = l*phi*I  (affine: phi fixed)
            sum( b(u,:) ) == l(u) * phi(u) * D(u);
            % per-slot relay capacity (tau_tbs fixed)
            % (kept as a soft check; b shape free here)
        end
        for m = 1:params.M
            users_m = find(X.assign == m);
            for n = 1:N
                sum( b(users_m,n) ) <= X.tau_tbs(m,n) * X.Rrelay_cache(m,n);
            end
        end
cvx_end

if strcmp(cvx_status,'Solved') || contains(cvx_status,'Inaccurate')
    X.l = min(max(l, params.l_min), params.l_max);
    X.b = max(b, 0);
else
    warning('SP1 Block-3 status: %s -- keeping previous l.', cvx_status);
end
end


% ======================================================================
%  Block 4 : relay ratio {phi}  (convex CVX program)
% ======================================================================
function X = block_relay_ratio(X, params, O, D, km)
% With l fixed, the bilinear product l*phi becomes AFFINE in phi.
% Update {phi, b} jointly:
%   minimise E3(phi)   (E3 decreases as phi increases -> more sent to TBS)
% subject to (P1c) and per-slot relay capacity.
U = params.U; N = params.N;
l = X.l;  f_uav_of = zeros(U,1);
for u = 1:U, f_uav_of(u) = X.f_uav(X.assign(u)); end

cvx_clear
cvx_begin quiet
    variable phi(U)
    variable b(U, N) nonnegative

    e3 = 0;
    for u = 1:U
        e3 = e3 + km * f_uav_of(u)^2 * O * l(u) * (1 - phi(u)) * D(u);
    end
    % E4 proxy: relaying more bits costs relay time; discourage extreme phi
    minimize( e3 )

    subject to
        phi >= params.phi_min;  phi <= params.phi_max;
        for u = 1:U
            sum( b(u,:) ) == l(u) * phi(u) * D(u);     % (P1c) affine in phi
        end
        for m = 1:params.M
            users_m = find(X.assign == m);
            for n = 1:N
                sum( b(users_m,n) ) <= X.tau_tbs(m,n) * X.Rrelay_cache(m,n);
            end
        end
cvx_end

if strcmp(cvx_status,'Solved') || contains(cvx_status,'Inaccurate')
    X.phi = min(max(phi, params.phi_min), params.phi_max);
    X.b   = max(b, 0);
else
    warning('SP1 Block-4 status: %s -- keeping previous phi.', cvx_status);
end
end
