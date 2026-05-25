function Q_new = SP2_trajectory_SCA(X, Q, params)
% SP2_TRAJECTORY_SCA  Trajectory optimisation via SCA (Algorithm 2).
%
% Minimises propulsion energy E5 plus a small communication-awareness
% term that pulls UAVs toward their user clusters and the TBS.
% Without this term, the straight-line trajectory is trivially optimal
% for pure propulsion minimisation, producing flat convergence curves.
% The weight w_comm is kept small (0.001) so propulsion still dominates.

M = params.M; N = params.N; dt = params.dt;
assign = X.assign;
Q_new = Q;

for sca = 1:params.sca_iter
    for m = 1:M

        users_m = find(assign == m);
        if isempty(users_m), continue; end

        % reference trajectory for this SCA pass
        q_ref = squeeze(Q_new(m, :, :))';
        if size(q_ref,1) ~= N, q_ref = q_ref'; end

        % reference speeds for propulsion majorisation
        v0 = params.v_min * ones(N, 1);
        for n = 2:N
            v0(n) = max( norm(q_ref(n,:) - q_ref(n-1,:)) / dt, params.v_min );
        end

        % user cluster centroid for this UAV
        user_centroid = mean(params.user_pos(users_m, :), 1);  % [1 x 2]

        cvx_clear
        cvx_begin quiet

            variable q(N, 2)
            variable s(N) nonnegative

            % --- propulsion energy via majorisation around v0 ---
            E5 = 0;
            for n = 2:N
                v_lin     = s(n) / dt;
                e_profile = params.k1 * pow_pos(v_lin, 3);
                e_induced = params.k2 * (2 / v0(n) - v_lin / (v0(n)^2));
                E5 = E5 + (e_profile + e_induced) * dt;
            end

            % --- communication attraction (small weight) ---
            % pulls trajectory toward user cluster and TBS
            comm_attract = 0;
            for n = 1:N
                comm_attract = comm_attract + ...
                    sum_square(q(n,:) - user_centroid) + ...
                    0.3 * sum_square(q(n,:) - params.GS);
            end

            w_comm = 0.001;
            minimize( E5 + w_comm * comm_attract )

            subject to

                % boundary constraints
                q(1, :) == params.q_start(m, :);
                q(N, :) == params.q_end(m, :);

                % speed limit
                s(1) == 0;
                for n = 2:N
                    norm( q(n,:) - q(n-1,:) ) <= s(n);
                    s(n) <= params.v_max * dt;
                end

                % collision avoidance (first-order linearisation)
                for j = 1:M
                    if j == m, continue; end
                    qj_ref = squeeze(Q_new(j, :, :))';
                    if size(qj_ref,1) ~= N, qj_ref = qj_ref'; end
                    for n = 1:N
                        dref = q_ref(n,:) - qj_ref(n,:);
                        lhs  = sum(dref.^2) + 2 * dref * (q(n,:)' - q_ref(n,:)');
                        lhs >= params.d_min_sq;
                    end
                end

        cvx_end

        if strcmp(cvx_status,'Solved') || contains(cvx_status,'Inaccurate')
            for n = 1:N
                Q_new(m, :, n) = q(n, :);
            end
        else
            warning('SP2 UAV %d SCA %d status: %s -- keeping previous Q.', ...
                    m, sca, cvx_status);
        end

    end  % UAV loop
end      % SCA loop

end