function [Q, X, hist, E_break] = bcd_sca_solve(params, relay_mode, verbose)
% BCD_SCA_SOLVE  Joint BCD-SCA solver for problem P1 (Algorithm 3).
%
%   Alternates:
%     SP1 : resource allocation + relay scheduling  (inner BCD)
%     SP2 : trajectory optimisation                 (SCA)
%   until the relative decrease of total energy falls below conv_tol.
%
% Inputs:
%   params      - parameter struct
%   relay_mode  - 'adaptive' (proposed) or 'uniform' (baseline)
%   verbose     - (optional) true to print per-iteration log
%
% Outputs:
%   Q       - final trajectory
%   X       - final resource allocation
%   hist    - energy per outer iteration  (honest, NOT post-processed)
%   E_break - final energy breakdown struct

if nargin < 2 || isempty(relay_mode), relay_mode = 'adaptive'; end
if nargin < 3, verbose = true; end

Q = trajectory_init(params);
X = [];
hist = zeros(params.max_iter, 1);
E_prev = inf;

for k = 1:params.max_iter

    % ---- SP1: resource allocation + relay scheduling ------------------
    if strcmpi(relay_mode, 'uniform')
        X = SP1_uniform_relay(Q, params);
    else
        X = SP1_resource_allocation(Q, params, X);
    end

    % ---- SP2: trajectory optimisation via SCA -------------------------
    [Q, X, sp2diag] = SP2_trajectory_SCA(X, Q, params);

    % ---- evaluate true total energy -----------------------------------
    [E_tot, E_break] = energy_model(X, Q, params);
    hist(k) = E_tot;

    if verbose
        fprintf(['  iter %2d | E=%.4e J | E1=%.3e E2=%.3e E3=%.3e ', ...
                 'E4=%.3e E5=%.3e\n'], k, E_tot, ...
                 E_break.E1, E_break.E2, E_break.E3, E_break.E4, E_break.E5);
    end

    % ---- convergence test ---------------------------------------------
    if k > 1
        rel = abs(E_prev - E_tot) / max(abs(E_prev), 1);
        if rel < params.conv_tol
            hist = hist(1:k);
            if verbose
                fprintf('  converged at iteration %d (rel change %.2e)\n', k, rel);
            end
            break;
        end
    end
    E_prev = E_tot;

end

hist = hist(1:k);

end
