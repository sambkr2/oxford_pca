%% Load data
load('../dataformatting/x20190517_CrossTumble_CR12p5_T1_C33_DVA_Motored_Processed_all.mat')

%% Parameters setting
AnalysisResult.CrankAngle = [ -285 ];                                   % Change this line to allow more crank angles (avaiable from -295 to -60 CAD aTDCf)
AnalysisResult.CycleNo = 1:300;                                   % Do not change this line

[ ~, AnalysisResult.CrankAngleIndex ] = ismember( AnalysisResult.CrankAngle, PODData.CrankAngle );

t_cad = num2str(AnalysisResult.CrankAngle);
ccm_cad = strrep(t_cad,'-','m');

%% POD Analysis
AnalysisResult.PODResult = cell( length( AnalysisResult.CrankAngleIndex ), 1 );
for ca_No = 1 : length( AnalysisResult.CrankAngleIndex )
    % Perform POD
    CurrentCrankAngle = AnalysisResult.CrankAngle( ca_No );
    fprintf( 'CA = %.0f CAD aTDCf \n', CurrentCrankAngle )

    temp_velo_data = complex( PODData.U{ AnalysisResult.CrankAngleIndex( ca_No ) }, PODData.V{ AnalysisResult.CrankAngleIndex( ca_No ) } );
    temp_PODResult = Perform_POD_SB( temp_velo_data, 'Centered', 'Direct' );

    temp_PODResult.X_PODGrid = PODData.X{ AnalysisResult.CrankAngleIndex( ca_No ) };
    temp_PODResult.Y_PODGrid = PODData.Y{ AnalysisResult.CrankAngleIndex( ca_No ) };
    temp_PODResult.IndexInOriginalGrid = PODData.IndexInOriginal{ AnalysisResult.CrankAngleIndex( ca_No ) };
    temp_PODResult.nRowsInOriginalGrid = size( MaskedData.X, 1 );
    temp_PODResult.nColsInOriginalGrid = size( MaskedData.X, 2 );
    temp_PODResult.X_OriginalGrid = MaskedData.X;
    temp_PODResult.Y_OriginalGrid = MaskedData.Y;

    AnalysisResult.PODResult{ ca_No } = temp_PODResult;
end

PODResult = AnalysisResult.PODResult;

%% Gavish Donoho
svs = diag(AnalysisResult.PODResult{1,1}.svdS);
beta = PODResult{1,1}.nColsInOriginalGrid/PODResult{1,1}.nRowsInOriginalGrid;
tau = optimal_SVHT_coef(beta,0) * median(svs); % find cut-off tau
modes = svs(svs>tau);
GDmode = length(modes); % Gavish Donoho threshold mode

%% Fourier
a = PODResult{1,1}.svdS * PODResult{1,1}.svdV';
maxa = max(abs(a),[],'all');
norma = a/maxa;
Vt = PODResult{1,1}.svdV';

dt = 0.08;
t = (0:N-1)*dt;     % time vector
Fs = 1/dt;      % sampling frequency
afft = fft(a,[],2);
Vfft = fft(Vt,[],2);
N = 300;  % number of samples
freqs = Fs*(0:(N/2))/N;
dimt = ones(1,length(t));
dimf = ones(1,length(freqs)); % initialise POD number dim

% Plot time series of POD modes
figure1 = figure( 'Color', [ 1 1 1 ] );
axes('Parent', figure1, 'FontName', 'Times', 'FontSize', 20,'linewidth',1.5)
view([1,-1,3])
grid on
hold on
for ii = 1:5
    plot3(t,dimt*ii,norma(ii,:), 'LineWidth',1.5)
end
yl = ylabel('POD modes');
yl.Position = yl.Position + [-1 0.8 0];
xl = xlabel('Time (s)');
xl.Position = xl.Position + [-5 0.2 0];
zlabel('Normalised amplitude')

xpower = abs(afft(:,1:N/2+1))*2/N;
maxpower = max(abs(xpower),[],'all');
normpower = xpower/maxpower;

Vxpower = abs(Vfft(:,1:N/2+1))*2/N;
Vmaxpower = max(abs(Vxpower),[],'all');
Vnormpower = Vxpower/Vmaxpower;

% Plot frequency series of POD modes
figure2 = figure( 'Color', [ 1 1 1 ] );
axes('Parent', figure2,'FontName', 'Times', 'FontSize', 20,'linewidth',1.5)
view([1,-1,3])
grid on
hold on
for jj = 1:5
    plot3(freqs,dimf*jj,normpower(jj,:),'LineWidth',1.5)
end
xlim([0 6.1])
xticks([0 1 2 3 4 5 6])
yticks([1 2 3 4 5])
yl2 = ylabel('POD modes');
yl2.Position = yl2.Position + [-0.3 0.8 0];
xl2 = xlabel('Frequency (Hz)');
xl2.Position = xl2.Position + [-1 0.2 0];
zlabel('Normalised amplitude')

% Plot frequency series of POD modes
figure3 = figure( 'Color', [ 1 1 1 ] );
axes('Parent', figure3,'FontName', 'Times', 'FontSize', 20,'linewidth',1.5)
view([1,-1,3])
grid on
hold on
for jj = 1:5
    plot3(freqs,dimf*jj,Vnormpower(jj,:),'LineWidth',1.5)
end
xlim([0 6.1])
xticks([0 1 2 3 4 5 6])
yticks([1 2 3 4 5])
yl2 = ylabel('POD modes');
yl2.Position = yl2.Position + [-0.3 0.8 0];
xl2 = xlabel('Frequency (Hz)');
xl2.Position = xl2.Position + [-1 0.2 0];
zlabel('Normalised amplitude')

%% Save Fourier
name1 = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/timeseries_podmodes.png'];
name2 = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/fseries_podmodes.png'];

exportgraphics(figure1,name1,'resolution',600)
exportgraphics(figure2,name2,'resolution',600)

%% POD approx paramters
% nModes = [ 1 2 5 10 50 299 ];
nModes = 0;
CycleNo = 36;
cycstr = num2str(CycleNo);
cyc = ['cyc',cycstr];

%% 
figureprop.axes_lim = [ -25 25 -20 2 ];
figureprop.xlabel = '{\it y} (mm)';
figureprop.ylabel = '{\it z} (mm)';
% figureprop.ylabel = [];
figureprop.velocity_normalisation = 5;
figureprop.sparse_vector = 2;
figureprop.Clim = [ 0 50 ];

%%
% Check before overwriting images due to automatic save

%% Plot POD approx
for mm = 1 : length( nModes )
    [ PODApprox ] = Calc_PODApprox( PODResult{1,1}, nModes(mm), CycleNo );
%     PODApprox.X = MaskedData.X;
%     PODApprox.Y = MaskedData.Y;
    PODApprox.X = InterpolatedData.X;
    PODApprox.Y = InterpolatedData.Y;
    PODVel.(['POD',num2str(nModes(mm))]).u = PODApprox.U;
    PODVel.(['POD',num2str(nModes(mm))]).v = PODApprox.V;
    
    figure_output = ColourQuiver( PODApprox.X, PODApprox.Y, PODApprox.U, PODApprox.V, figureprop );
    set(gca,'fontsize',26)
%     annotation('line', [0.5 0.55],[0.645 0.695],'linewidth',7,'color','red')
%     annotation('line', [0.5 0.55],[0.695 0.645],'linewidth',7,'color','red')

%     sparse_vector = 3;
% 
%     quivx = PODApprox.X(1:sparse_vector:end);
%     quivy = PODApprox.Y(1:sparse_vector:end);
%     quivu = PODApprox.U(1:sparse_vector:end);
%     quivv = PODApprox.V(1:sparse_vector:end);
% 
%     angleMap_quiver = angleMap( 1:sparse_vector:end, 1:sparse_vector:end );
%     speedMap_quiver = speedMap( 1:sparse_vector:end, 1:sparse_vector:end );
% 
%     figure1 = figure('color',[1 1 1]);
%     q = quiver(quivx,quivy,quivu,quivv);
%     q.LineWidth = 1.2;
%     set(gca,'fontname','times','fontsize',26,'linewidth',1.2)
%     set(gcf, 'position', [357 379 762 497])
%     set(gcf,'OuterPosition',[357 379 762 497])
%     set(gca,'Position',[0.1225 0.1100 0.7302 0.8150])
%     axis equal
%     xlabel('{\it y} (mm)');
%     ylabel('{\it z} (mm)');
%     yticks([-100:10:100])
%     ylim([-20 2])
%     xlim([-25 25])
%     pivname = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/CTP/CTP_',cyc,'_',ccm_cad,'_POD',num2str( nModes(mm) ),'.png'];
%     exportgraphics(gcf,pivname,'resolution',600)

%     title( [ 'POD Approx., Order = ', num2str( nModes(mm) )] );
%     title(['PIV ensemble mean ',num2str(CurrentCrankAngle),' CAD'])
%     t = {['PIV cycle A, ',num2str(CurrentCrankAngle),' CAD'],['POD Approx., Order = ', num2str( nModes(mm) )]};
%     title(t)
%     export_fig( [ 'TP Cycle ', num2str( cycle_No ), ' POD Approx at -270 CAD aTDCf' ], '-pdf', '-nocrop', '-append' )
%     close all
end

%%
pivname = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/CTP/CTP_',cyc,'_',ccm_cad,'_POD',num2str( nModes(mm) ),'_annotated.png'];
exportgraphics(gcf,pivname,'resolution',600)

%% multi-cycle
nModes = [ 299 ];
CycleNo = [1:1:50];
% CycleNo = [35 36];

%% Plot POD approx cycles
for mm = 1 : length( CycleNo )
    [ PODApprox ] = Calc_PODApprox( PODResult{1,1}, nModes, CycleNo(mm) );
    cyc = num2str(CycleNo(mm));
%     PODApprox.X = MaskedData.X;
%     PODApprox.Y = MaskedData.Y;
    PODApprox.X = InterpolatedData.X;
    PODApprox.Y = InterpolatedData.Y;
%     PODVel.(['POD',num2str(nModes(mm))]).u = PODApprox.U;
%     PODVel.(['POD',num2str(nModes(mm))]).v = PODApprox.V;
    
    figure_output = ColourQuiver( PODApprox.X, PODApprox.Y, PODApprox.U, PODApprox.V, figureprop );
    ylim([-20 2])
%     set(gcf, 'position', [440 377 560 290])
%     pivname = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/CTP/CTP_cyc',cyc,'_',ccm_cad,'_POD',num2str( nModes ),'.png'];
%     exportgraphics(gcf,pivname,'resolution',600)

%     title( [ 'POD Approx., Order = ', num2str( nModes(mm) )] );
%     title(['PIV ensemble mean ',num2str(CurrentCrankAngle),' CAD'])
%     t = {['PIV cycle A, ',num2str(CurrentCrankAngle),' CAD'],['POD Approx., Order = ', num2str( nModes(mm) )]};
%     title(t)
%     export_fig( [ 'TP Cycle ', num2str( cycle_No ), ' POD Approx at -270 CAD aTDCf' ], '-pdf', '-nocrop', '-append' )
    %     close all
    set(gcf, 'position', [357 379 762 497])
    set(gcf,'OuterPosition',[357 379 762 497])
    set(gca,'Position',[0.1225 0.1100 0.7302 0.8150])
end

%%
pivname = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/CTP/CTP_',cyc,'_',ccm_cad,'_POD',num2str( nModes(mm) ),'.png'];
exportgraphics(gcf,pivname,'resolution',600)

%% Load CFD
load('rrT1rng/ccm_T1_mot_CTP.mat')

%% Parameters setting
AnalysisResult.CrankAngle = [ -285 ];                                   % Change this line to allow more crank angles (avaiable from -295 to -60 CAD aTDCf)
AnalysisResult.CycleNo = 1:300;                                   % Do not change this line

[ ~, AnalysisResult.CrankAngleIndex ] = ismember( AnalysisResult.CrankAngle, PODData.CrankAngle );

t_cad = num2str(AnalysisResult.CrankAngle);
ccm_cad = strrep(t_cad,'-','m');

%% Plot CFD 
temp_x = CFDData.Data.y_PIVGrid;
temp_y = CFDData.Data.z_PIVGrid;
temp_CCM_u = ccmdata.(ccm_cad).v;
temp_CCM_v = ccmdata.(ccm_cad).w;
temp_CCM_SpeedMap = abs( complex( temp_CCM_u, temp_CCM_v ) );

% Add PIV mask to CFD 
[ ~, PIV_CAindex ] = ismember( AnalysisResult.CrankAngle, InterpolatedData.CrankAngle );

temp_PIVem_u = mean( InterpolatedData.U( :,:,PIV_CAindex,: ), 4, 'omitnan' );
temp_PIVem_v = mean( InterpolatedData.V( :,:,PIV_CAindex,: ), 4, 'omitnan' );
temp_PIVem_SpeedMap = abs( complex( temp_PIVem_u, temp_PIVem_v ) );

temp_PIV_mask = ~isnan( temp_PIVem_SpeedMap );
temp_PIV_mask = double( temp_PIV_mask );
temp_PIV_mask( temp_PIV_mask==0 ) = NaN;
% temp_CCM_u = temp_CCM_u .* temp_PIV_mask;
% temp_CCM_v = temp_CCM_v .* temp_PIV_mask;
% temp_CCM_SpeedMap = temp_CCM_SpeedMap .* temp_PIV_mask;

ColourQuiver(temp_x, temp_y, temp_CCM_u ,temp_CCM_v, figureprop)
% ylim([-20 2])
% title(['CFD ',num2str(CurrentCrankAngle),' CAD'])
set(gcf, 'position', [440 377 560 290])

%%
name = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/CTP/CTP_CFD_',ccm_cad,'.png'];
exportgraphics(gcf,name,'resolution',600)
