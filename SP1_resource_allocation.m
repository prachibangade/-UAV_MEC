function X = SP1_resource_allocation(Q, params, X_init)
M = params.M; U = params.U; N = params.N; dt = params.dt;
O  = params.cycles_per_bit;
D  = params.D(:);
ku = params.kappa_user;
km = params.kappa_uav;

% user -> UAV assignment
Q_mean = squeeze(mean(Q, 3));
assign = zeros(U, 1);
for u = 1:U
    d2 = sum((Q_mean - params.user_pos(u, :)).^2, 2);
    [~, assign(u)] = min(d2);
end

% rates for fixed trajectory
[~, Ruser] = channel_model(Q, params, params.P_u * ones(U, 1));
Rrelay = relay_rate(Q, params);

Ru = zeros(U, N);
for u = 1:U
    Ru(u, :) = squeeze(Ruser(assign(u), u, :))';
end
Ru = max(Ru, 1e2);

if nargin >= 3 && ~isempty(X_init)
    X = X_init;
    X.assign = assign;
else
    X = init_feasible_X(assign, Ru, Rrelay, params);
end
X.Rrelay_cache = Rrelay;

for it = 1:params.sp1_iter
    X = block_frequencies(X, params, O, D, ku, km);
    X = block_times_bits(X, params, D, Ru, Rrelay);
    X = block_offload_ratio(X, params, D, Ru);
    X = block_relay_ratio(X, params, O, D, km);
end

X.assign = assign;
X.Ru     = Ru;
X.Rrelay = Rrelay;
end


function X = init_feasible_X(assign, Ru, Rrelay, params)
M = params.M; U = params.U; N = params.N; dt = params.dt;
D = params.D(:);

X.assign = assign;
X.l   = 0.5 * ones(U, 1);
X.phi = 0.5 * ones(U, 1);

X.tau_off = zeros(U, N);
for u = 1:U
    need = X.l(u) * D(u);
    t_each = need / max(sum(Ru(u, :)), 1);
    X.tau_off(u, :) = min(t_each, 0.45 * dt);
end

X.tau_tbs = zeros(M, N);
X.b       = zeros(U, N);
for u = 1:U
    relay_bits = X.l(u) * X.phi(u) * D(u);
    X.b(u, :)  = relay_bits / N;
end
for m = 1:M
    users_m = find(assign == m);
    for n = 1:N
        bsum = sum(X.b(users_m, n));
        X.tau_tbs(m, n) = min(bsum / max(Rrelay(m, n), 1), 0.45 * dt);
    end
end

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


function X = block_frequencies(X, params, O, D, ku, km)
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


function X = block_times_bits(X, params, D, Ru, Rrelay)
M = params.M; U = params.U; N = params.N; dt = params.dt;
l = X.l; phi = X.phi; assign = X.assign;

cvx_clear
cvx_begin quiet
    variable tau_off(U, N) nonnegative
    variable tau_tbs(M, N) nonnegative
    variable b(U, N) nonnegative

    minimize( params.P_u * sum(sum(tau_off)) + params.P_uav * sum(sum(tau_tbs)) )

    subject to
        for u = 1:U
            sum( tau_off(u, :) .* Ru(u, :) ) >= l(u) * D(u);
        end
        for u = 1:U
            sum( b(u, :) ) == l(u) * phi(u) * D(u);
        end
        for m = 1:M
            users_m = find(assign == m);
            for n = 1:N
                sum( b(users_m, n) ) <= tau_tbs(m, n) * Rrelay(m, n);
            end
        end
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


function X = block_offload_ratio(X, params, D, Ru)
% Closed-form update for l using (P1c) and (P1b)
U = params.U;
for u = 1:U
    phi_u = max(X.phi(u), params.phi_min);
    b_sum = sum(X.b(u,:));
    l_from_c = b_sum / max(phi_u * D(u), 1);
    l_cap = sum(X.tau_off(u,:) .* Ru(u,:)) / max(D(u), 1);
    X.l(u) = min(max(min(l_from_c, l_cap), params.l_min), params.l_max);
end
end


function X = block_relay_ratio(X, params, O, D, km)
U = params.U; N = params.N;
l = X.l;
f_uav_of = zeros(U, 1);
for u = 1:U
    f_uav_of(u) = X.f_uav(X.assign(u));
end

cvx_clear
cvx_begin quiet
    variable phi(U)
    variable b(U, N) nonnegative

    e3 = 0;
    for u = 1:U
        e3 = e3 + km * f_uav_of(u)^2 * O * l(u) * (1 - phi(u)) * D(u);
    end
    minimize( e3 )

    subject to
        phi >= params.phi_min;
        phi <= params.phi_max;
        for u = 1:U
            sum( b(u,:) ) == l(u) * phi(u) * D(u);
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