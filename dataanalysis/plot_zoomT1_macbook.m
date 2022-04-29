figure1 = figure( 'Color', [ 1 1 1 ] );
axes1 = axes( 'Parent', figure1, 'FontName', 'serif', 'FontSize', 14,'linewidth',1);
hold( axes1, 'on' )
box( axes1, 'on')
plot( -360:1:360, EngineData.PIV.stats.PcylMean, 'Parent', axes1,'color','black','linestyle','-.','linewidth',2, 'DisplayName', 'Exp mean')
plot(T1press.cad,T1press.p,'parent',axes1,'color','red','linewidth',2, 'DisplayName', 'CFD')
% plot(pmotoredRKE.rke_cad,pmotoredRKE.rke_p,'parent',axes1,'color','red','linewidth',2, 'DisplayName', 'RKE')
% plot(pmotoredRKE2L.rke2l_cad,pmotoredRKE2L.rke2l_p,'parent',axes1,'color','blue','linewidth',2, 'DisplayName', 'Angelberger')
% plot(-360:1:360, quantile_PIV_PCyl(:,1),'k-','linewidth',1.2,'parent',axes1)
% plot(-360:1:360, quantile_PIV_PCyl(:,2),'k-','linewidth',1.2,'parent',axes1)
temp_PIV_PCyl = EngineData.PIV.cycles.filt.Pcyl;
std_PIV_PCyl = std(temp_PIV_PCyl,0,2);
plot_error_bar( -360:1:360, EngineData.PIV.stats.PcylMean, std_PIV_PCyl, 2, 'shaded', [ 0.3 0.3 0.3 ], 'Exp +/- 2 std.' )

axis( [-360 180 0 10] )
set( axes1, 'xtick', -360:90:360 )
xlabel( 'Crank Angle (CAD aTDCf)' )
ylabel( 'In-Cylinder Pressure (bar)' )
l = legend('show');
set(l,'location','northeast')
% title( 'T1 C33 DVA Piston B' )

zoomin_color = [0.30 0.75 0.93];
zoomin_xlimMin = -12.5;
zoomin_xlimMax = 12.5;
zoomin_ylimMin = 8.8;
zoomin_ylimMax = 9.4;
line( [ zoomin_xlimMin, zoomin_xlimMax ], [ zoomin_ylimMin, zoomin_ylimMin ], 'Parent', axes1, 'Color', zoomin_color, 'LineStyle', '-', 'LineWidth', 1.5, 'HandleVisibility', 'off' )
line( [ zoomin_xlimMin, zoomin_xlimMax ], [ zoomin_ylimMax, zoomin_ylimMax ], 'Parent', axes1, 'Color', zoomin_color, 'LineStyle', '-', 'LineWidth', 1.5, 'HandleVisibility', 'off' )
line( [ zoomin_xlimMin, zoomin_xlimMin ], [ zoomin_ylimMin, zoomin_ylimMax ], 'Parent', axes1, 'Color', zoomin_color, 'LineStyle', '-', 'LineWidth', 1.5, 'HandleVisibility', 'off' )
line( [ zoomin_xlimMax, zoomin_xlimMax ], [ zoomin_ylimMin, zoomin_ylimMax ], 'Parent', axes1, 'Color', zoomin_color, 'LineStyle', '-', 'LineWidth', 1.5, 'HandleVisibility', 'off' )

axes2 = axes( 'Parent', figure1, 'FontName', 'serif', 'FontSize', 12, 'Position', [ 0.24 0.45 0.3 0.4 ],'linewidth',1);
hold( axes2, 'on' )
box( axes2, 'on' )
plot( -360:1:360, EngineData.PIV.stats.PcylMean, 'Parent', axes2,'LineStyle', '-.','color','black','linewidth',2)
plot(T1press.cad,T1press.p,'parent',axes2,'color','red','linewidth',2, 'DisplayName', 'CFD')
% plot(pmotoredRKE.rke_cad,pmotoredRKE.rke_p,'parent',axes2,'color','red','linewidth',2, 'DisplayName', 'RKE')
% plot(pmotoredRKE2L.rke2l_cad,pmotoredRKE2L.rke2l_p,'parent',axes2,'color','blue','linewidth',2, 'DisplayName', 'Angelberger')
% plot(-360:1:360, quantile_PIV_PCyl(:,1),'k-','linewidth',1.2,'parent',axes2)
% plot(-360:1:360, quantile_PIV_PCyl(:,2),'k-','linewidth',1.2,'parent',axes2)
plot_error_bar( -360:1:360, EngineData.PIV.stats.PcylMean, std_PIV_PCyl, 2, 'shaded', [ 0.3 0.3 0.3 ], [] )

axis( [zoomin_xlimMin zoomin_xlimMax zoomin_ylimMin zoomin_ylimMax] )
xlabel( {'Crank Angle','(CAD aTDCf)'} )
ylabel( 'Cylinder {\it{P}} (bar)' )
% set( axes2, 'xtick', -360:10:360, 'ytick', 0:0.5:15 )
axes2.XGrid = 'on';
axes2.YGrid = 'on';
% Create lines
annotation(figure1,'line',[0.54 0.63], [0.45 0.825], 'Color', zoomin_color, 'LineStyle', '--');
annotation(figure1,'line',[0.54 0.63], [0.85 0.875], 'Color', zoomin_color, 'LineStyle', '--');