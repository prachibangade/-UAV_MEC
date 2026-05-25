function X = SP1_uniform_relay(Q, params)
% SP1_UNIFORM_RELAY  Baseline resource allocation with UNIFORM relay
% scheduling. Identical to SP1_resource_allocation except that the
% per-slot relay bits b[n] are forced to be equal across all slots
% (b[u,n] = l*phi*I / N) instead of being optimised.
%
% This isolates the contribution of the proposed adaptive relay
% scheduling: trajectory and offloading are optimised the same way,
% only b[n] differs.

% First get a normal SP1 solution (for l, phi, f, assignment).
X = SP1_resource_allocation(Q, params);

M = params.M; U = params.U; N = params.N; dt = params.dt;
D = params.D(:);
Rrelay = X.Rrelay;
assign = X.assign;

% Force uniform b[n].
X.b = zeros(U, N);
for u = 1:U
    relay_bits = X.l(u) * X.phi(u) * D(u);
    X.b(u, :)  = relay_bits / N;            % equal split
end

% Recompute relay time tau_tbs to carry the uniform b under the per-slot
% capacity:  tau_tbs[m,n] = sum_u b[u,n] / Rrelay[m,n].
X.tau_tbs = zeros(M, N);
for m = 1:M
    users_m = find(assign == m);
    for n = 1:N
        bsum = sum(X.b(users_m, n));
        X.tau_tbs(m, n) = bsum / Rrelay(m, n);
    end
end

% Note: uniform scheduling ignores time-varying relay rate, so it spends
% more relay time (hence energy) in low-rate slots. That is exactly the
% inefficiency the proposed adaptive scheme removes.

end
