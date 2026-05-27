function [max_viol, detail] = sp2_feasibility(X, Q, params)
% SP2_FEASIBILITY  Maximum constraint violation of a candidate solution
% evaluated against the TRUE (nonlinear) rates and the real geometry.
%
%   SP2 solves a convex surrogate in which the rates are frozen at the
%   reference trajectory. After SP2 returns a candidate (Q, tau, b) we
%   must check the candidate against the genuine constraints before
%   accepting it. This function returns the largest violation (<=0 means
%   feasible up to tolerance) across all SP2 constraints.
%
% Inputs:
%   X      - resource struct: needs l, phi, assign, tau_off, tau_tbs, b
%   Q      - candidate trajectory [M x 2 x N]
%   params - parameter struct
%
% Outputs:
%   max_viol - scalar, max_i ( g_i(x) ) where g_i <= 0 is the i-th
%              constraint; positive => infeasible
%   detail   - struct of the worst violation per constraint family

M = params.M; U = params.U; N = params.N; dt = params.dt;
D = params.D(:);
assign = X.assign;

% true rates at the candidate trajectory
[~, Ruser_full] = channel_model(Q, params, params.P_u * ones(U,1));
Rrelay = relay_rate(Q, params);
Ru = zeros(U, N);
for u = 1:U
    Ru(u, :) = squeeze(Ruser_full(assign(u), u, :))';
end

v = zeros(1, 0);

% ---- (33d) offload delivery:  l*D - sum(tau_off.*Ru) <= 0 -------------
g_off = -inf;
for u = 1:U
    g = X.l(u)*D(u) - sum(X.tau_off(u,:) .* Ru(u,:));
    g_off = max(g_off, g);
end

% ---- (33e) relay total bits:  | sum(b) - l*phi*D | <= tol -------------
g_relaytot = -inf;
for u = 1:U
    g = abs(sum(X.b(u,:)) - X.l(u)*X.phi(u)*D(u));
    g_relaytot = max(g_relaytot, g);
end

% ---- (33f) relay capacity:  sum_u b - tau_tbs*Rrelay <= 0 -------------
g_cap = -inf;
for m = 1:M
    um = find(assign == m);
    for n = 1:N
        g = sum(X.b(um,n)) - X.tau_tbs(m,n)*Rrelay(m,n);
        g_cap = max(g_cap, g);
    end
end

% ---- (33g) relay causality: cumsum(b) - cumsum_prev(tau_off.*Ru) <=0 --
g_caus = -inf;
for u = 1:U
    delivered = X.tau_off(u,:) .* Ru(u,:);
    cum_recv  = [0, cumsum(delivered(1:N-1))];   % bits available by slot n
    cum_relay = cumsum(X.b(u,:));
    g_caus = max(g_caus, max(cum_relay - cum_recv));
end

% ---- (33j) per-slot time budget: sum(tau_off)+tau_tbs - dt <= 0 -------
g_budget = -inf;
for m = 1:M
    um = find(assign == m);
    for n = 1:N
        g = sum(X.tau_off(um,n)) + X.tau_tbs(m,n) - dt;
        g_budget = max(g_budget, g);
    end
end

% ---- (33k) speed limit: ||q[n+1]-q[n]|| - v_max*dt <= 0 ---------------
g_speed = -inf;
for m = 1:M
    for n = 1:N-1
        step = norm(squeeze(Q(m,:,n+1)) - squeeze(Q(m,:,n)));
        g_speed = max(g_speed, step - params.v_max*dt);
    end
end

% ---- (33m) collision: d_min^2 - ||q_m-q_j||^2 <= 0 --------------------
g_coll = -inf;
for m = 1:M
    for j = m+1:M
        for n = 1:N
            d2 = sum((squeeze(Q(m,:,n)) - squeeze(Q(j,:,n))).^2);
            g_coll = max(g_coll, params.d_min_sq - d2);
        end
    end
end
if M < 2, g_coll = -inf; end

detail.offload    = g_off;
detail.relay_tot  = g_relaytot;
detail.relay_cap  = g_cap;
detail.causality  = g_caus;
detail.budget     = g_budget;
detail.speed      = g_speed;
detail.collision  = g_coll;

max_viol = max([g_off, g_relaytot, g_cap, g_caus, g_budget, g_speed, g_coll]);

end
