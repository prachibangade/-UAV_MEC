function plot_results(results)
p = results.params;

C.proposed = [0     0.447 0.741];
C.local    = [0.850 0.325 0.098];
C.full     = [0.466 0.674 0.188];
C.static   = [0.494 0.184 0.556];
C.uniform  = [0.929 0.694 0.125];
lw = 1.8;  ms = 7;  fs = 12;

% FIG 2: convergence
f = figure('Position',[100 100 560 420]);
plot(1:numel(results.hist_prop), results.hist_prop, '-o', ...
     'Color', C.proposed, 'LineWidth', lw, 'MarkerSize', ms, ...
     'MarkerFaceColor', C.proposed);
xlabel('BCD Iteration Index','FontSize',fs);
ylabel('Total Energy Consumption (J)','FontSize',fs);
title('Convergence of Proposed BCD-SCA Algorithm','FontSize',fs+1);
grid on; box on; set(gca,'FontSize',fs-1);
exportgraphics(f,'fig2_convergence.png','Resolution',300);
fprintf('  saved fig2_convergence.png\n');

% FIG 3: trajectories
f = figure('Position',[150 100 560 420]); hold on;
scatter(p.user_pos(:,1), p.user_pos(:,2), 70, [0.3 0.3 0.3], ...
        'filled', 'DisplayName','Users');
plot(p.GS(1), p.GS(2), 'p', 'MarkerSize',18, 'MarkerFaceColor',C.uniform, ...
     'MarkerEdgeColor','k', 'DisplayName','TBS');
uav_c = {C.proposed, C.local};
Q = results.Q_prop;
for m = 1:p.M
    xm = squeeze(Q(m,1,:));  ym = squeeze(Q(m,2,:));
    plot(xm, ym, '-', 'Color', uav_c{min(m,2)}, 'LineWidth', 2.2, ...
         'DisplayName', sprintf('UAV %d path', m));
    plot(xm(1),  ym(1),  'o', 'Color',uav_c{min(m,2)}, 'MarkerSize',10, ...
         'MarkerFaceColor',uav_c{min(m,2)}, 'HandleVisibility','off');
    plot(xm(end),ym(end),'s', 'Color',uav_c{min(m,2)}, 'MarkerSize',10, ...
         'MarkerFaceColor',uav_c{min(m,2)}, 'HandleVisibility','off');
end
xlabel('x (m)','FontSize',fs); ylabel('y (m)','FontSize',fs);
title('Optimized UAV Trajectories','FontSize',fs+1);
legend('Location','best','FontSize',fs-2);
grid on; box on; axis equal; set(gca,'FontSize',fs-1);
exportgraphics(f,'fig3_trajectory.png','Resolution',300);
fprintf('  saved fig3_trajectory.png\n');

% FIG 4: energy vs number of users
f = figure('Position',[200 100 560 420]); hold on;
labels = {'Proposed','Local only','Full offload','Static UAV','Uniform relay'};
cols   = {C.proposed, C.local, C.full, C.static, C.uniform};
mk     = {'o','s','d','^','v'};
for j = 1:size(results.E_vs_U,1)
    plot(results.U_list, results.E_vs_U(j,:), '-', 'Color',cols{j}, ...
         'Marker',mk{j}, 'LineWidth',lw, 'MarkerSize',ms, ...
         'MarkerFaceColor',cols{j}, 'DisplayName',labels{j});
end
xlabel('Number of Users U','FontSize',fs);
ylabel('Total Energy Consumption (J)','FontSize',fs);
title('Total Energy vs Number of Users','FontSize',fs+1);
legend('Location','northwest','FontSize',fs-2);
grid on; box on; set(gca,'FontSize',fs-1);
exportgraphics(f,'fig4_energy_vs_users.png','Resolution',300);
fprintf('  saved fig4_energy_vs_users.png\n');

% FIG 5: energy vs task size
f = figure('Position',[250 100 560 420]); hold on;
labels = {'Proposed','Full offload','Static UAV'};
cols   = {C.proposed, C.full, C.static};
mk     = {'o','d','^'};
for j = 1:size(results.E_vs_I,1)
    plot(results.I_list/1e6, results.E_vs_I(j,:), '-', 'Color',cols{j}, ...
         'Marker',mk{j}, 'LineWidth',lw, 'MarkerSize',ms, ...
         'MarkerFaceColor',cols{j}, 'DisplayName',labels{j});
end
xlabel('Task Data Size I_u (Mbit)','FontSize',fs);
ylabel('Total Energy Consumption (J)','FontSize',fs);
title('Total Energy vs Task Data Size','FontSize',fs+1);
legend('Location','northwest','FontSize',fs-2);
grid on; box on; set(gca,'FontSize',fs-1);
exportgraphics(f,'fig5_energy_vs_tasksize.png','Resolution',300);
fprintf('  saved fig5_energy_vs_tasksize.png\n');

% FIG 7: energy breakdown
f = figure('Position',[300 100 620 420]);
schemes = {'Proposed','Local','Full','Static','Uniform'};
EB = [ struct2vec(results.Eb_prop);
       struct2vec(results.Eb_local);
       struct2vec(results.Eb_full);
       struct2vec(results.Eb_static);
       struct2vec(results.Eb_uniform) ];
bar(EB, 'stacked');
set(gca,'XTickLabel',schemes,'FontSize',fs-1);
ylabel('Energy (J)','FontSize',fs);
title('Energy Breakdown by Component','FontSize',fs+1);
legend({'E1 local comp','E2 offload','E3 UAV comp','E4 relay','E5 propulsion'}, ...
       'Location','eastoutside','FontSize',fs-3);
grid on; box on;
exportgraphics(f,'fig7_energy_breakdown.png','Resolution',300);
fprintf('  saved fig7_energy_breakdown.png\n');

% FIG 8: relay energy E4 comparison only -- FIX: use categorical x-axis
f = figure('Position',[350 100 560 420]);
E4_vals = [results.Eb_prop.E4, results.Eb_uniform.E4];
b = bar(E4_vals);
b(1).FaceColor = C.proposed;
b(1).CData = [C.proposed; C.uniform];
b(1).FaceColor = 'flat';
set(gca,'XTickLabel',{'Adaptive (proposed)','Uniform (baseline)'},'FontSize',fs-1);
ylabel('Relay Energy E_4 (J)','FontSize',fs);
title('Adaptive vs Uniform Relay Scheduling','FontSize',fs+1);
% add percentage label
pct = 100*(E4_vals(2)-E4_vals(1))/E4_vals(2);
text(1.5, max(E4_vals)*0.6, sprintf('%.1f%% saving', pct), ...
     'HorizontalAlignment','center','FontSize',fs,'FontWeight','bold','Color','k');
grid on; box on;
exportgraphics(f,'fig8_relay_compare.png','Resolution',300);
fprintf('  saved fig8_relay_compare.png\n');

end


function v = struct2vec(Eb)
v = [Eb.E1, Eb.E2, Eb.E3, Eb.E4, Eb.E5];
end