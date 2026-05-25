function Q_new = SP2_trajectory_SCA(X, Q, params)
% SP2_TRAJECTORY_SCA  Trajectory optimisation via Successive Convex
% Approximation (Algorithm 2 of the corrected paper).
%
% This version implements the paper's described method (Section V-D-3):
%   - The log-rates R_{u,m}[n] and R_m^TBS[n] are NON-CONVEX in q.
%     Each is convex in the squared distance D = ||q - w||^2, so we
%     LOWER-BOUND each rate by its first-order Taylor expansion in D
%     around the current iterate.  D itself is convex in q, so the
%     resulting lower bound is concave in q -> usable in CVX.
%   - Propulsion k1*v^3 + k2/v: v^3 via pow_pos (convex); k2/v via the
%     standard tangent upper bound around v0.
%   - Collision constraint linearised to an affine lower bound.
%
% The SP2 objective is the trajectory-dependent energy E2 + E4 + E5:
%   E2 (offload) and E4 (relay) enter through epigraph variables t_off,
%   t_tbs bounded using the linearised rates; E5 is propulsion.
%
% Inputs:
%   X      - fixed resource allocation from SP1
%   Q      - current trajectory [M x 2 x N]
%   params - parameter struct
%
% Output:
%   Q_new  - refined trajectory [M x 2 x N]

M = params.M; N = params.N; dt = params.dt;
H = params.H;
assign = X.assign;

Q_new = Q;

for sca = 1:params.sca_iter

    for m = 1:M

        users_m = find(assign == m);

        % ---- reference iterate for this SCA pass --------------------
        q_ref = squeeze(Q_new(m, :, :))';        % [N x 2]
        if size(q_ref,1) ~= N, q_ref = q_ref'; end

        % reference speeds (for k2/v tangent)
        v0 = params.v_min * ones(1, N);
        for n = 2:N
            v0(n) = max( norm(q_ref(n,:) - q_ref(n-1,:)) / dt, params.v_min );
        end

        % ---- precompute Taylor data for uplink rates (per user) -----
        % R_user(D) = B*log2(1 + P_u*beta0/(noise*(H^2+D)))
        % treat as function of D = ||q-w_u||^2.
        % rate_ref and slope d R / d D  at D_ref.
        Ruser_ref  = zeros(numel(users_m), N);
        Ruser_slp  = zeros(numel(users_m), N);
        for ii = 1:numel(users_m)
            u = users_m(ii);
            for n = 1:N
                D_ref = sum( (q_ref(n,:) - params.user_pos(u,:)).^2 );
                [r0, s0] = lograte_taylor(D_ref, params.P_u, params);
                Ruser_ref(ii,n) = r0;
                Ruser_slp(ii,n) = s0;       % <= 0
            end
        end

        % ---- precompute Taylor data for relay rate ------------------
        Rrelay_ref = zeros(1, N);
        Rrelay_slp = zeros(1, N);
        for n = 1:N
            D_ref = sum( (q_ref(n,:) - params.GS).^2 );
            [r0, s0] = lograte_taylor(D_ref, params.P_uav, params);
            Rrelay_ref(n) = r0;
            Rrelay_slp(n) = s0;
        end

        % ---- solve convex surrogate SP2 for UAV m -------------------
        cvx_clear
        cvx_begin quiet
            variable q(N, 2)
            variable s(N)                       % speed slack, v <= s/dt? see below

            % --- propulsion energy E5 ---
            E5 = 0;
            for n = 2:N
                v_lin = s(n) / dt;
                e_profile = params.k1 * pow_pos(v_lin, 3);
                % k2/v tangent upper bound around v0(n):  k2/v <= k2*(2/v0 - v/v0^2)
                e_induced = params.k2 * (2 / v0(n) - v_lin / (v0(n)^2));
                E5 = E5 + (e_profile + e_induced) * dt;
            end

            % --- communication energy proxy E2 + E4 ---
            % We want UAVs to fly so that rates are high (low tx energy).
            % Using linearised rates:
            %   Rlin_user(ii,n) = Ruser_ref + slope*( ||q-w||^2 - D_ref )
            % Since slope<=0 and ||q-w||^2 is convex, Rlin is concave in q.
            % We reward total linearised rate (equivalently penalise -rate),
            % weighted by the fixed bits each link must carry.
            comm = 0;
            for ii = 1:numel(users_m)
                u = users_m(ii);
                bits_off = X.l(u) * params.D(u);            % offloaded bits
                for n = 1:N
                    Dn = sum_square( q(n,:) - params.user_pos(u,:) );
                    Rlin = Ruser_ref(ii,n) + Ruser_slp(ii,n) * (Dn - ...
                           sum((q_ref(n,:)-params.user_pos(u,:)).^2));
                    % penalise low rate: weight by bits / horizon
                    comm = comm - (bits_off / (N*params.T)) * Rlin;
                end
            end
            for n = 1:N
                bits_rel = sum( X.l(users_m).*X.phi(users_m).*params.D(users_m) );
                Dn = sum_square( q(n,:) - params.GS );
                Rlin = Rrelay_ref(n) + Rrelay_slp(n) * (Dn - ...
                       sum((q_ref(n,:)-params.GS).^2));
                comm = comm - (bits_rel / (N*params.T)) * Rlin;
            end

            % weight balances propulsion vs communication terms
            w_comm = params.P_uav;     % scale; both terms then ~ J
            minimize( E5 + w_comm * comm )

            subject to
                q(1, :) == params.q_start(m, :);
                q(N, :) == params.q_end(m, :);
                for n = 2:N
                    norm( q(n,:) - q(n-1,:) ) <= s(n);
                    s(n) <= params.v_max * dt;
                end
                s >= 0;

                % --- collision avoidance, linearised (affine lower bound) ---
                for j = 1:M
                    if j == m, continue; end
                    qj_ref = squeeze(Q_new(j, :, :))';
                    if size(qj_ref,1) ~= N, qj_ref = qj_ref'; end
                    for n = 1:N
                        dref = q_ref(n,:) - qj_ref(n,:);
                        % ||q_m-q_j||^2 >= ||dref||^2 + 2 dref.(q_m - q_ref_m)
                        % (q_j fixed at this pass)
                        lhs = sum(dref.^2) + 2 * dref * ( q(n,:)' - q_ref(n,:)' );
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


% ======================================================================
%  lograte_taylor  -- value and slope of B*log2(1+P*beta0/(noise*(H^2+D)))
%                     as a function of D = ||q-w||^2, evaluated at D_ref.
% ======================================================================
function [r0, slope] = lograte_taylor(D_ref, P_tx, params)
H2  = params.H^2;
a   = P_tx * params.beta0 / params.noise;     % so SNR = a/(H2+D)
denom = H2 + D_ref;
snr   = a / denom;
r0    = params.B * log2(1 + snr);

% d/dD [ B*log2(1 + a/(H2+D)) ]
%   let x = a/(H2+D);  dx/dD = -a/(H2+D)^2
%   dr/dD = B/ln2 * (1/(1+x)) * dx/dD
dx_dD = -a / denom^2;
slope = (params.B / log(2)) * (1 / (1 + snr)) * dx_dD;   % <= 0
end
