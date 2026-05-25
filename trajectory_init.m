function Q = trajectory_init(params)
% TRAJECTORY_INIT  Initial feasible UAV trajectory.
%
%   Straight line from q_start to q_end for each UAV, sampled at N slots.
%   This is a feasible starting point for the BCD-SCA loop (it satisfies
%   the boundary constraints and, for the chosen geometry, the speed limit).

M = params.M;
N = params.N;
Q = zeros(M, 2, N);

for m = 1:M
    qs = params.q_start(m, :);
    qe = params.q_end(m, :);
    for n = 1:N
        t = (n - 1) / max(N - 1, 1);
        Q(m, :, n) = (1 - t) * qs + t * qe;
    end
end

end
