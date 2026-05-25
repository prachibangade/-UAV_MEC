function [E, E_break, Q, X] = baseline_scheme(name, params)
% BASELINE_SCHEME  Compute total energy for a named baseline.
%
%   name (string):
%     'local'      - all tasks computed locally (l = l_min ~ 0)
%     'full'       - all offloaded (l = l_max ~ 1), trajectory via SCA
%     'static'     - UAVs hover at midpoint, resource allocation optimised
%     'equal'      - fixed l = 0.5, initial straight-line trajectory
%     'uniform'    - full BCD-SCA but UNIFORM relay scheduling (key ablation)
%
% Outputs: total energy E, breakdown struct, trajectory Q, allocation X.

switch lower(name)

    case 'local'
        % minimal offloading -> almost all local computation
        Q = trajectory_init(params);
        p = params;  p.l_max = p.l_min + 1e-3;
        X = SP1_resource_allocation(Q, p);
        [E, E_break] = energy_model(X, Q, params);

    case 'full'
        % force maximum offloading, optimise trajectory
        p = params;  p.l_min = p.l_max - 1e-3;
        Q = trajectory_init(p);
        X = SP1_resource_allocation(Q, p);
        Q = SP2_trajectory_SCA(X, Q, p);
        X = SP1_resource_allocation(Q, p);
        [E, E_break] = energy_model(X, Q, params);

    case 'static'
        % UAVs hover at midpoint of start/end (no trajectory optimisation)
        Q = static_trajectory(params);
        X = SP1_resource_allocation(Q, params);
        [E, E_break] = energy_model(X, Q, params);

    case 'equal'
        % fixed l = 0.5, initial trajectory, single resource solve
        Q = trajectory_init(params);
        p = params;  p.l_min = 0.5; p.l_max = 0.5;
        X = SP1_resource_allocation(Q, p);
        [E, E_break] = energy_model(X, Q, params);

    case 'uniform'
        % proposed framework but uniform relay scheduling
        [Q, X, ~, E_break] = bcd_sca_solve(params, 'uniform', false);
        [E, ~] = energy_model(X, Q, params);

    otherwise
        error('baseline_scheme: unknown baseline "%s"', name);
end

end


function Q = static_trajectory(params)
% UAVs hover at the midpoint of their start/end positions.
M = params.M; N = params.N;
Q = zeros(M, 2, N);
for m = 1:M
    mid = 0.5 * (params.q_start(m,:) + params.q_end(m,:));
    for n = 1:N
        Q(m, :, n) = mid;
    end
end
end
