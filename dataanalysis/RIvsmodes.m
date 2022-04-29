figure1 = figure( 'Color', [ 1 1 1 ] );
axes1 = axes( 'Parent', figure1, 'FontName', 'serif', 'FontSize', 16,'linewidth',1);
hold( axes1, 'on' )
box( axes1, 'on')
bar( 1:1:299, n_ke,'r')
xlabel( 'Number of POD modes','FontSize',22)
ylabel({'Fraction of kinetic'; 'energy captured [-]'},'FontSize',22 )