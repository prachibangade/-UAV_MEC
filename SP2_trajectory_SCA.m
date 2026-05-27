function [Q_out, X_out, diag] = SP2_trajectory_SCA(X, Q, params)
% SP2_TRAJECTORY_SCA  Corrected trajectory subproblem (Algorithm 2).
%
%   SP2 minimises the TRAJECTORY-DEPENDENT restriction of the true total
%   system energy:
%
%       min   E2 + E4 + E5
%        q,tau_off,tau_tbs,b
%
%   where, with {l, phi, f} fixed from SP1, E1 and E3 are constants and
%   are dropped (dropping a constant does not change the minimiser, which
%   is exactly what keeps SP2 consistent with SP1: both blocks minimise
%   restrictions of the SAME E_tot = E1+E2+E3+E4+E5).
%
%       E2 = P_u   * sum_{u,n} tau_off[u,n]      (MD -> UAV offload)
%       E4 = P_uav * sum_{m,n} tau_tbs[m,n]      (UAV -> TBS relay)
%       E5 = sum_{m,n} dt*( k1 v^3 + k2 / v )    (propulsion)
%
%   OPTIMISED here : q, tau_off, tau_tbs, b, and propulsion speed slacks.
%   FIXED here      : l, phi, f_loc, f_uav  (held at SP1's output).
%
%   Non-convexities and how they are handled:
%     - log-rate in the offload/relay constraints -> Option A: rates are
%       evaluated at the SCA reference trajectory (rate_surrogate.m) and
%       treated as constants within one CVX solve, refreshed each pass.
%       The Taylor expansion of a convex-decreasing rate is a global
%       lower bound, so this is conservative.
%     - propulsion k1 v^3 + k2/v : written EXACTLY as pow_pos(v,3) and
%       inv_pos(v). Both are convex for v>0; NO majorisation is used or
%       needed. v >= v_min guards inv_pos.
%     - collision ||q_m-q_j||^2 >= d_min^2 : reverse-convex, linearised
%       by a first-order lower bound around the reference (valid: tangent
%       of a convex function lies below it, so the surrogate implies the
%       true constraint).
%
%   Stabilisation:
%     - per-slot trust region ||q[n]-q_ref[n]|| <= Delta, adapted by a
%       predicted-vs-actual ratio test.
%     - monotonicity safeguard: a candidate is ACCEPTED only if it is
%       truly feasible (checked against the real nonlinear rates) AND it
%       does not increase the true total energy. Otherwise it is rejected
%       and the trust region is shrunk.
%
%   NO heuristic attraction term. NO silent fallback. Solver failures are
%   surfaced and handled explicitly.
%
% Inputs:
%   X      - resource struct from SP1. Must contain:
%              X.l, X.phi          [U x 1]  offload / relay ratios (FIXED)
%              X.assign            [U x 1]  user -> UAV assignment
%              X.f_loc             [U x 1]  MD CPU frequency  (FIXED)
%              X.f_uav             [M x 1]  UAV CPU frequency (FIXED)
%              X.tau_off, X.tau_tbs, X.b    warm-start values (optional)
%   Q      - incoming trajectory [M x 2 x N] (BCD reference)
%   params - parameter struct
%
% Outputs:
%   Q_out  - refined trajectory  [M x 2 x N]
%   X_out  - X with tau_off, tau_tbs, b updated by SP2 (l,phi,f unchanged)
%   diag   - per-SCA-iteration diagnostics struct array

M = params.M; U = params.U; N = params.N; dt = params.dt;
D = params.D(:);
assign = X.assign;

% ----- trust-region parameters ----------------------------------------
Delta      = 1.0 * params.v_max * dt;   % initial radius (m)
Delta_min  = 0.05 * params.v_max * dt;  % stop when radius this small
Delta_max  = 5.0  * params.v_max * dt;
eta_lo     = 0.25;                      % ratio test thresholds
eta_hi     = 0.75;
gamma_dec  = 0.5;                       % shrink factor
gamma_inc  = 2.0;                       % expand factor
feas_tol   = 1e-3;                      % true-feasibility tolerance
sca_max    = params.sca_iter;

% ----- accepted incumbent ---------------------------------------------
Q_acc = Q;
X_acc = X;
% energy of the incumbent (true). Uses X's current tau/b -- if SP1 just
% ran, these are SP1's values, which are the correct baseline.
[~, Eb_acc] = energy_model(X_acc, Q_acc, params);
E_acc = Eb_acc.E2 + Eb_acc.E4 + Eb_acc.E5;

diag = struct('sca',{},'cvx_status',{},'surrogate',{}, ...
              'E_true',{},'pred_dec',{},'act_dec',{}, ...
              'ratio',{},'max_viol',{},'Delta',{}, ...
              'accepted',{},'viol_detail',{});
for j = 1:sca_max

    % ===== 1. reference rates at the current incumbent ================
    S = rate_surrogate(Q_acc, params, assign);   % S.Ruser, S.Rrelay
    Rref_u   = S.Ruser;     % [U x N]
    Rref_tbs = S.Rrelay;    % [M x N]

    % reference geometry for the collision linearisation
    Qref = Q_acc;

    % ===== 2. solve the convex SP2 surrogate ==========================
    cvx_clear
    cvx_begin quiet
        variable q(M, 2, N)
        variable s(M, N)        nonnegative      % per-slot displacement
        variable v(M, N)        nonnegative      % per-slot speed
        variable tau_off(U, N)  nonnegative
        variable tau_tbs(M, N)  nonnegative
        variable b(U, N)        nonnegative

        % ---- objective: E2 + E4 + E5 ---------------------------------
        E2 = params.P_u  * sum(sum(tau_off));
        E4 = params.P_uav * sum(sum(tau_tbs));
        E5 = 0;
        for m = 1:M
            for n = 2:N
                % exact convex propulsion: no linearisation
                E5 = E5 + dt * ( params.k1 * pow_pos(v(m,n), 3) ...
                               + params.k2 * inv_pos(v(m,n)) );
            end
        end
        minimize( E2 + E4 + E5 )

        subject to

            % ---- boundary conditions (33l) --------------------------
            for m = 1:M
                reshape(q(m,:,1),1,2) == params.q_start(m,:);
                reshape(q(m,:,N),1,2) == params.q_end(m,:);
            end

            % ---- speed / displacement / kinematics (33k) ------------
            for m = 1:M
                v(m,1) == params.v_min;          % no slot-1 transition
                s(m,1) == 0;
                for n = 2:N
                    step = reshape(q(m,:,n) - q(m,:,n-1), 1, 2);
                    norm( step ) <= s(m,n);
                    s(m,n) <= v(m,n) * dt;       % v is the speed of slot n
                    v(m,n) >= params.v_min;      % guards inv_pos
                    v(m,n) <= params.v_max;
                end
            end

            % ---- offload delivery (33d), linear in tau_off ----------
            % rate frozen at reference -> affine constraint
            for u = 1:U
                sum( tau_off(u,:) .* Rref_u(u,:) ) >= X.l(u) * D(u);
            end

            % ---- relay total bits (33e), l & phi fixed -> linear ----
            for u = 1:U
                sum( b(u,:) ) == X.l(u) * X.phi(u) * D(u);
            end

            % ---- relay per-slot capacity (33f), linear in tau_tbs ---
            for m = 1:M
                um = find(assign == m);
                for n = 1:N
                    sum( b(um,n) ) <= tau_tbs(m,n) * Rref_tbs(m,n);
                end
            end

            % ---- relay causality (33g) ------------------------------
            % cumulative relayed bits up to slot n cannot exceed
            % cumulative offloaded bits received by slot n-1
            for u = 1:U
                delivered = tau_off(u,:) .* Rref_u(u,:);   % affine row
                for n = 1:N
                    if n == 1
                        sum( b(u,1) ) <= 0;                % nothing yet
                    else
                        sum( b(u,1:n) ) <= sum( delivered(1:n-1) );
                    end
                end
            end

            % ---- per-slot time budget (33j) -------------------------
            for m = 1:M
                um = find(assign == m);
                for n = 1:N
                    sum( tau_off(um,n) ) + tau_tbs(m,n) <= dt;
                end
            end

            % ---- collision avoidance (33m), first-order lower bound -
            % f(q) = ||q_m-q_jj||^2 is convex; its tangent at the
            % reference is a global lower bound, so enforcing
            %   f(qref) + grad.(q-qref) >= d_min^2
            % implies the true non-convex constraint f(q) >= d_min^2.
            for m = 1:M
                for jj = m+1:M
                    for n = 1:N
                        % reference difference, forced to a 1x2 ROW
                        drow = reshape(Qref(m,:,n) - Qref(jj,:,n), 1, 2);
                        % current difference, a CVX expression (1x2 row)
                        dqrow = reshape(q(m,:,n) - q(jj,:,n), 1, 2);
                        % drow*(dqrow-drow).'  ->  scalar inner product
                        sum(drow.^2) ...
                            + 2 * drow * (dqrow - drow).' ...
                            >= params.d_min_sq;
                    end
                end
            end

            % ---- trust region --------------------------------------
            for m = 1:M
                for n = 1:N
                    dtr = reshape(q(m,:,n), 1, 2) - ...
                          reshape(Qref(m,:,n), 1, 2);
                    norm( dtr ) <= Delta;
                end
            end

    cvx_end

    status    = cvx_status;
    surrogate = cvx_optval;

    % ===== 3. handle solver outcome ===================================
    solved = strcmpi(status,'Solved') || ...
             ~isempty(strfind(lower(status),'inaccurate'));

    if ~solved
        % NO silent fallback. Surface it, shrink the trust region,
        % and retry on the next SCA pass from the same incumbent.
        rec = make_diag(j, status, surrogate, E_acc, NaN, NaN, NaN, ...
                        NaN, Delta, false);
        rec.viol_detail = struct();
        diag(end+1) = rec; %#ok<AGROW>
        print_sp2_sca_diag(rec, NaN);
        warning('SP2:solverFailed', ...
            ['SP2 SCA %d: CVX status = %s. Shrinking trust region ', ...
             'Delta %.2f -> %.2f and retrying.'], ...
             j, status, Delta, gamma_dec*Delta);
        Delta = max(gamma_dec * Delta, Delta_min);
        if Delta <= Delta_min
            warning('SP2:trustRegionCollapsed', ...
                ['SP2: trust region collapsed without a solve at SCA ', ...
                 '%d. Returning last accepted incumbent.'], j);
            break;
        end
        continue;
    end

    % ===== 4. assemble candidate and check TRUE feasibility ===========
    X_cand          = X_acc;          % inherits l, phi, f, assign
    X_cand.tau_off  = max(tau_off, 0);
    X_cand.tau_tbs  = max(tau_tbs, 0);
    X_cand.b        = max(b, 0);
    Q_cand          = q;

    [max_viol, vdetail] = sp2_feasibility(X_cand, Q_cand, params);

    % ===== 5. monotonicity safeguard ==================================
    [~, Eb_cand] = energy_model(X_cand, Q_cand, params);
    E_true = Eb_cand.E2 + Eb_cand.E4 + Eb_cand.E5;

    pred_dec = E_acc - surrogate;
    act_dec  = E_acc - E_true;                      % actual drop
    if abs(pred_dec) < 1e-12
        ratio = 1.0;            % nothing predicted; treat as neutral
    else
        ratio = act_dec / pred_dec;
    end

    feasible_true = (max_viol <= feas_tol);
    improves      = (E_true <= E_acc - 1e-9);
    accept        = feasible_true && improves;
    traj_move     = norm(Q_cand(:) - Q_acc(:));

    % ===== 6. trust-region update + accept/reject =====================
    if accept
        Q_acc = Q_cand;
        X_acc = X_cand;
        E_acc = E_true;
        if ratio > eta_hi
            Delta = min(gamma_inc * Delta, Delta_max);   % expand
        elseif ratio < eta_lo
            Delta = max(gamma_dec * Delta, Delta_min);   % shrink
        end
        % eta_lo <= ratio <= eta_hi -> keep Delta
    else
        % reject: incumbent unchanged, shrink trust region
        Delta = max(gamma_dec * Delta, Delta_min);
    end

    rec = make_diag(j, status, surrogate, E_true, pred_dec, act_dec, ...
                    ratio, max_viol, Delta, accept);
    rec.viol_detail = vdetail;
    diag(end+1) = rec; %#ok<AGROW>
    print_sp2_sca_diag(rec, traj_move);

    % ===== 7. stopping tests ==========================================
    if accept && abs(act_dec) < params.conv_tol * max(abs(E_acc),1)
        break;   % converged: accepted step but negligible improvement
    end
    if Delta <= Delta_min
        break;   % trust region exhausted
    end

end

Q_out = Q_acc;
X_out = X_acc;

end


% ======================================================================
function rec = make_diag(j, status, surrogate, E_true, pred_dec, ...
                         act_dec, ratio, max_viol, Delta, accepted)
% MAKE_DIAG  Pack one SCA-iteration diagnostic record.
rec.sca        = j;
rec.cvx_status = status;
rec.surrogate  = surrogate;
rec.E_true     = E_true;
rec.pred_dec   = pred_dec;
rec.act_dec    = act_dec;
rec.ratio      = ratio;
rec.max_viol   = max_viol;
rec.Delta      = Delta;
rec.accepted   = accepted;
end


% ======================================================================
function print_sp2_sca_diag(rec, traj_move)
% PRINT_SP2_SCA_DIAG  Console diagnostics for one SP2 SCA iteration.
fprintf('\n--- SP2 SCA iter %d ---\n', rec.sca);
fprintf('  cvx_status:          %s\n', rec.cvx_status);
fprintf('  surrogate objective: %.6e\n', rec.surrogate);
if isnan(rec.E_true)
    fprintf('  E_true:              NaN\n');
else
    fprintf('  E_true:              %.6e\n', rec.E_true);
end
if isnan(rec.pred_dec)
    fprintf('  predicted decrease:  NaN\n');
else
    fprintf('  predicted decrease:  %.6e\n', rec.pred_dec);
end
if isnan(rec.act_dec)
    fprintf('  actual decrease:     NaN\n');
else
    fprintf('  actual decrease:     %.6e\n', rec.act_dec);
end
if isnan(rec.max_viol)
    fprintf('  max feas violation:  NaN\n');
else
    fprintf('  max feas violation:  %.6e\n', rec.max_viol);
end
fprintf('  Delta:               %.4f\n', rec.Delta);
if rec.accepted
    fprintf('  accepted:            accepted\n');
else
    fprintf('  accepted:            rejected\n');
end
if isnan(traj_move)
    fprintf('  ||Q_new-Q_old||:     NaN\n');
else
    fprintf('  ||Q_new-Q_old||:     %.6e\n', traj_move);
end
end