function [G, R] = channel_model(Q, params, P_tx)
% CHANNEL_MODEL  Air-to-ground channel gains and achievable uplink rates.
%
%   MD-to-UAV channel (Stage 1):
%       h_{u,m}[n] = beta0 / ( H^2 + || q_m[n] - w_u ||^2 )
%       R_{u,m}[n] = B * log2( 1 + P_u * h_{u,m}[n] / noise )
%
% Inputs:
%   Q      - UAV trajectory  [M x 2 x N]
%   params - parameter struct
%   P_tx   - (optional) per-user transmit power [U x 1]; defaults to params.P_u
%
% Outputs:
%   G - channel gains  [M x U x N]
%   R - data rates     [M x U x N]   (bits/s)

M = params.M;
U = params.U;
N = params.N;
H = params.H;

if nargin < 3 || isempty(P_tx)
    P_tx = params.P_u * ones(U, 1);
end
P_tx = P_tx(:);

G = zeros(M, U, N);
R = zeros(M, U, N);

for n = 1:N
    for m = 1:M
        q_mn = Q(m, :, n);
        for u = 1:U
            w_u = params.user_pos(u, :);
            d2  = sum((q_mn - w_u).^2);
            g   = params.beta0 / (H^2 + d2);
            G(m, u, n) = g;
            R(m, u, n) = params.B * log2(1 + P_tx(u) * g / params.noise);
        end
    end
end

end
