figure1 = figure( 'Color', [ 1 1 1 ] );
axes1 = axes( 'Parent', figure1, 'FontName', 'serif', 'FontSize', 14,'linewidth',1);
hold( axes1, 'on' )
box( axes1, 'on')
plot(riT1rkeSie(:,1),riT1rkeSie(:,2),'linewidth',1.2)
plot(riT1(:,1),riT1(:,2),'linewidth',1.2)
% plot(riT7(:,1),riT7(:,2),'linewidth',1.2)
% xlim([-270 cadFirst+cadStep*(length(files)-1)])
xlim([-270 -35])
ylim([0.5 1])
xlabel('Crank angle degrees aTDCf [CAD]')
ylabel('Relevance index [-]')
legend('RKE2L init','RKE2L def','location','northwest')
set(axes1, 'xtick', -360:45:360)