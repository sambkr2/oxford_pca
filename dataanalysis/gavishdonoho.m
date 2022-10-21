%% Load tumble
load('../dataformatting/x20180706_Tumble_CR12p5_T1_C33_DVA_Motored_Processed_all_masked.mat');

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
% InterpolatedData = myData.InterpolatedData;

%% Gavish Donoho
svs = diag(AnalysisResult.PODResult{1,1}.svdS);
beta = PODResult{1,1}.nColsInOriginalGrid/PODResult{1,1}.nRowsInOriginalGrid;
tau = optimal_SVHT_coef(beta,0) * median(svs); % find cut-off tau
modes = svs(svs>tau);
GDmode = length(modes); % Gavish Donoho threshold mode

%% POD approx paramters
% nModes = [ 0 1 2 10 GDmode 299 ];
% nModes = [ 1 2 5 10 50 299 ];
nModes = [ 0 ];
CycleNo = 53;
cycstr = num2str(CycleNo);
cyc = ['cyc',cycstr];

%% 
figureprop.axes_lim = [ -25 25 -30 10 ];
figureprop.xlabel = '{\it x} (mm)';
figureprop.ylabel = '{\it z} (mm)';
figureprop.velocity_normalisation = 5;
figureprop.sparse_vector = 2;
figureprop.Clim = [ 0 50 ];

%% Plot POD approx

for mm = 1 : length( nModes )
    [ PODApprox ] = Calc_PODApprox( PODResult{1,1}, nModes(mm), CycleNo );
%     PODApprox.X = MaskedData.X;
%     PODApprox.Y = MaskedData.Y;
    PODApprox.X = InterpolatedData.X;
    PODApprox.Y = InterpolatedData.Y;
    PODVel.(['POD',num2str(nModes(mm))]).u = PODApprox.U;
    PODVel.(['POD',num2str(nModes(mm))]).v = PODApprox.V;
%     figure_output = ColourQuiver_SB( PODApprox.X, PODApprox.Y, PODApprox.U, PODApprox.V, figureprop );
    figure_output = ColourQuiver( PODApprox.X, PODApprox.Y, PODApprox.U, PODApprox.V, figureprop );    
    ylim([-30 10])
    set(gca,'fontsize',26)
%     title(['PIV ensemble mean ',num2str(CurrentCrankAngle),' CAD'])
%     title( [ 'POD Approx., Order = ', num2str( nModes(mm) )] );
%     annotation('line', [0.49 0.54],[0.3 0.35],'linewidth',7,'color','red')
%     annotation('line', [0.49 0.54],[0.35 0.3],'linewidth',7,'color','red')
end

%%
pivname = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/TP/TP_',cyc,'_',ccm_cad,'_POD',num2str( nModes(mm) ),'_annotated.png'];
exportgraphics(gcf,pivname,'resolution',600)

%% POD approx paramters
% nModes = [ 0 1 2 5 8 20 299 GDmode ];
nModes = [ 299 ];
CycleNo = [1:300];

%% 
figureprop.axes_lim = [ -25 25 -30 10 ];
figureprop.xlabel = '{\it x} (mm)';
figureprop.ylabel = '{\it z} (mm)';

figureprop.sparse_vector = 2;
figureprop.Clim = [ 0 50 ];

%% PODVel
for mm = 1 : length( CycleNo )
    [ PODApprox ] = Calc_PODApprox( PODResult{1,1}, nModes, CycleNo(mm) );
%     PODApprox.X = MaskedData.X;
%     PODApprox.Y = MaskedData.Y;
    PODApprox.X = InterpolatedData.X;
    PODApprox.Y = InterpolatedData.Y;
    PODVel.(['cyc',num2str(CycleNo(mm))]).u = PODApprox.U;
    PODVel.(['cyc',num2str(CycleNo(mm))]).v = PODApprox.V;
%     figure_output = ColourQuiver_SB( PODApprox.X, PODApprox.Y, PODApprox.U, PODApprox.V, figureprop );
%     ylim([-30 10])
%     title( [ num2str( CurrentCrankAngle ),' CAD, ', 'POD Approx., Order = ', num2str( nModes(mm) )] );
%     export_fig( [ 'TP Cycle ', num2str( cycle_No ), ' POD Approx at -270 CAD aTDCf' ], '-pdf', '-nocrop', '-append' )
%     close all
end

%% Calculate ensemble mean of GDmode
PODCyc = nan(PODData.nRowsInOriginal * PODData.nColsInOriginal, 300);

for ii = 1:300
    [ PODgd ] = Calc_PODApprox( PODResult{1,1}, 2, ii );
    PODgd.X = InterpolatedData.X;
    PODgd.Y = InterpolatedData.Y;

    % Turn velocity values into column vectors
    tempu = reshape(PODgd.U,[],1); 
    tempv = reshape(PODgd.V,[],1);
    
    % Store as complex number with each column as a new cycle
%     PODcyc.(['cyc',num2str(ii)]) = complex(tempu,tempv);
    PODCyc(:,ii) = complex(tempu,tempv);
    
    % Average along the rows to get the ensemble mean
    colGDem = mean(PODCyc,2);
    
    % Reshape into matrix
    GDem = reshape(colGDem, PODData.nRowsInOriginal, PODData.nColsInOriginal);
end

%% Plot GDem
% PODApprox.X = InterpolatedData.X;
% PODApprox.Y = InterpolatedData.Y;
gdU = real(GDem);
gdV = imag(GDem);
figure_output = ColourQuiver_SB( PODgd.X, PODgd.Y, gdU, gdV, figureprop );
ylim([-30 10])
line1 = ['Gavish Donoho ensemble mean'];
line2 = [t_cad,' CAD, ','POD Order = ', num2str( GDmode )];
title({line1,line2});
set(gcf, 'position', [440 352 560 446])

%% Plot CFD 
temp_x = CFDData.Data.x_PIVGrid;
temp_y = CFDData.Data.z_PIVGrid;
temp_CCM_u = ccmdata.(ccm_cad).u;
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
temp_CCM_u = temp_CCM_u .* temp_PIV_mask;
temp_CCM_v = temp_CCM_v .* temp_PIV_mask;
temp_CCM_SpeedMap = temp_CCM_SpeedMap .* temp_PIV_mask;

ColourQuiver(temp_x, temp_y, temp_CCM_u ,temp_CCM_v, figureprop)
ylim([-30 10])
% title(['CFD ',num2str(CurrentCrankAngle),' CAD'])
% set(gcf, 'position', [440 377 560 290])
name = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/TP/TP_CFD_',ccm_cad,'.png'];

%%
exportgraphics(gcf,name,'resolution',600)
