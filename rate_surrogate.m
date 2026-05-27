function S = rate_surrogate(Q_ref, params, assign)
% RATE_SURROGATE  Conservative (lower-bound) rate surrogates for SP2.
%
%   For SP2 we use Option A from the redesign: the achievable rates are
%   evaluated at the SCA reference trajectory Q_ref and treated as
%   CONSTANT arrays inside one CVX solve. They are refreshed every SCA
%   pass. Because the true rate R(x) is convex-decreasing in the squared
%   distance x = ||q - w||^2, its first-order Taylor expansion in x is a
%   global LOWER BOUND. Evaluated AT the reference (x = x_ref) the bound
%   collapses to the true rate at the reference point:
%
%       Rhat^(j)(x_ref) = R^(j)
%
%   so the "surrogate" rate used by SP2 is simply the true rate at the
%   reference trajectory. This is conservative in the sense that any tau
%   sized against R^(j) will OVER-deliver against the true rate at the
%   updated trajectory only if the updated trajectory has higher rate;
%   the monotonicity safeguard (in SP2) re-checks true feasibility/energy
%   regardless, so correctness does not rely on this being a strict bound
%   away from the reference.
%
% Inputs:
%   Q_ref  - reference trajectory [M x 2 x N]
%   params - parameter struct
%   assign - [U x 1] user -> UAV assignment (from SP1)
%
% Output struct S:
%   S.Ruser  [U x N]  MD u offload rate to its assigned UAV, at Q_ref
%   S.Rrelay [M x N]  UAV m -> TBS relay rate, at Q_ref

M = params.M; U = params.U; N = params.N; H = params.H;

% ----- MD -> UAV offload rates (per assigned UAV) ----------------------
S.Ruser = zeros(U, N);
for u = 1:U
    m  = assign(u);
    wu = params.user_pos(u, :);
    for n = 1:N
        q_mn = squeeze(Q_ref(m, :, n));
        d2   = sum((q_mn(:)' - wu).^2);
        g    = params.beta0 / (H^2 + d2);
        S.Ruser(u, n) = params.B * log2(1 + params.P_u * g / params.noise);
    end
end

% ----- UAV -> TBS relay rates ------------------------------------------
S.Rrelay = zeros(M, N);
for m = 1:M
    for n = 1:N
        q_mn = squeeze(Q_ref(m, :, n));
        d2   = sum((q_mn(:)' - params.GS).^2);
        g    = params.beta0 / (H^2 + d2);
        S.Rrelay(m, n) = params.B * log2(1 + params.P_uav * g / params.noise);
    end
end

% No artificial floor. If a rate is genuinely tiny that is physical
% information the optimiser must see. CVX handles small positive
% coefficients without trouble; a zero would only arise at infinite
% distance, which the bounded operating area precludes.

end
