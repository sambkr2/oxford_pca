%% Set variables
ke = PODResult{1,1}.KE_mode;
sum_ke = sum(ke);
n_ke = ke./sum_ke;

sv = diag(PODResult{1,1}.svdS);
sum_sv = sum(sv);
n_sv = sv./sum_sv;
cumn_sv = cumsum(sv)./sum_sv;
cum_sv = cumsum(n_sv);

%% Plot
figure1 = figure( 'Color', [ 1 1 1 ] );
axes1 = axes( 'Parent', figure1, 'FontName', 'serif', 'FontSize', 16,'linewidth',1);
hold( axes1, 'on' )
box( axes1, 'on')
bar( 1:1:299, n_ke,'r')
xlabel( 'Number of POD modes','FontSize',22)
ylabel({'Fraction of kinetic'; 'energy captured [-]'},'FontSize',22 )

figure2 = figure( 'Color', [ 1 1 1 ] );
axes2 = axes( 'Parent', figure2, 'FontName', 'serif', 'FontSize', 16,'linewidth',1);
hold( axes2, 'on' )
box( axes2, 'on')
plot(1:1:299, cum_sv, 'linewidth',1.5)
xlabel( 'Number of POD modes','FontSize',22)
ylabel( {'Cumulative fraction of'; 'kinetic energy captured [-]'},'FontSize',22 )
ylim([0 1])

% figure2 = figure( 'Color', [ 1 1 1 ] );
% axes2 = axes( 'Parent', figure2, 'FontName', 'serif', 'FontSize', 22,'linewidth',1);
% % hold( axes2, 'on' )
% box( axes2, 'on')
% bar(1:1:299, n_sv,'b')
% xlabel( 'Singular value' )
% ylabel( 'Fraction of variance captured' )

% axis( [-360 180 0 10] )
% set( axes1, 'xtick', -360:90:360 )
% l = legend('show');
% set(l,'location','northeast')
% title( 'T1 C33 DVA Piston B' )
