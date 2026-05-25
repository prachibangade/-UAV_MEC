function R_relay = relay_rate(Q, params)
% RELAY_RATE  UAV-to-TBS relay rate per slot, per UAV.
%
%   h_{m,0}[n] = beta0 / ( H^2 + || q_m[n] - w_TBS ||^2 )
%   R_m^TBS[n] = B * log2( 1 + P_uav * h_{m,0}[n] / noise )
%
% Input:
%   Q      - UAV trajectory [M x 2 x N]
%   params - parameter struct
%
% Output:
%   R_relay - relay rate [M x N]  (bits/s)

M = params.M;
N = params.N;
H = params.H;

R_relay = zeros(M, N);
for m = 1:M
    for n = 1:N
        q_mn = squeeze(Q(m, :, n));
        d2 = sum((q_mn(:)' - params.GS).^2);
        g  = params.beta0 / (H^2 + d2);
        R_relay(m, n) = params.B * log2(1 + params.P_uav * g / params.noise);
    end
end
R_relay = max(R_relay, 1e3);   % numerical floor

end
