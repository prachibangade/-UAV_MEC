function [E, E_break, Q, X] = baseline_scheme(name, params)
% BASELINE_SCHEME  Compute total energy for a named baseline.

switch lower(name)

    case 'local'
        % Local only: no offloading, UAVs hover at midpoint (zero flight energy)
        Q = zeros(params.M, 2, params.N);
        for m = 1:params.M
            mid = 0.5 * (params.q_start(m,:) + params.q_end(m,:));
            for n = 1:params.N
                Q(m,:,n) = mid;
            end
        end
        p = params; p.l_max = p.l_min + 1e-3;
        X = SP1_resource_allocation(Q, p);
        [E, E_break] = energy_model(X, Q, params);

    case 'full'
        % Force maximum offloading, optimise trajectory
        p = params; p.l_min = p.l_max - 1e-3;
        Q = trajectory_init(p);
        X = SP1_resource_allocation(Q, p);
        Q = SP2_trajectory_SCA(X, Q, p);
        X = SP1_resource_allocation(Q, p);
        [E, E_break] = energy_model(X, Q, params);

    case 'static'
        % UAVs hover at midpoint, resource allocation optimised
        Q = static_trajectory(params);
        X = SP1_resource_allocation(Q, params);
        [E, E_break] = energy_model(X, Q, params);

    case 'equal'
        % Fixed l = 0.5, initial trajectory
        Q = trajectory_init(params);
        p = params; p.l_min = 0.5; p.l_max = 0.5;
        X = SP1_resource_allocation(Q, p);
        [E, E_break] = energy_model(X, Q, params);

    case 'uniform'
        % Full BCD-SCA but uniform relay scheduling
        [Q, X, ~, E_break] = bcd_sca_solve(params, 'uniform', false);
        [E, ~] = energy_model(X, Q, params);

    otherwise
        error('baseline_scheme: unknown baseline "%s"', name);
end
end


function Q = static_trajectory(params)
M = params.M; N = params.N;
Q = zeros(M, 2, N);
for m = 1:M
    mid = 0.5 * (params.q_start(m,:) + params.q_end(m,:));
    for n = 1:N
        Q(m, :, n) = mid;
    end
end
end