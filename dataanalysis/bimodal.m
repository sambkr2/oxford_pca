%% Tumble
load('../dataformatting/x20180706_Tumble_CR12p5_T1_C33_DVA_Motored_Processed_all_masked.mat');
name = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/hist_tumble.png'];

% Index at (0,-22)
xind = 37;
zind = 37;

%% Cross-tumble
load('../dataformatting/x20190517_CrossTumble_CR12p5_T1_C33_DVA_Motored_Processed_all.mat')
name = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/hist_crosstumble.png'];

% Index at (0,-1)
xind = 39;
zind = 8;

%% POD
AnalysisResult.CrankAngle = [ -285 ];                                   % Change this line to allow more crank angles (avaiable from -295 to -60 CAD aTDCf)
AnalysisResult.CycleNo = 1:300;                                   % Do not change this line

[ ~, AnalysisResult.CrankAngleIndex ] = ismember( AnalysisResult.CrankAngle, PODData.CrankAngle );

t_cad = num2str(AnalysisResult.CrankAngle);
ccm_cad = strrep(t_cad,'-','m');

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

%% Analysis
nModes = [ 299 ];
CycleNo = [1:1:300];

% Horizontal u components
horiz = zeros(length(CycleNo),2);
for mm = 1 : length( CycleNo )
    [ PODApprox ] = Calc_PODApprox( PODResult{1,1}, nModes, CycleNo(mm) );
    PODApprox.X = InterpolatedData.X;
    PODApprox.Y = InterpolatedData.Y;
    horiz(mm,1) = mm;
    horiz(mm,2) = PODApprox.U(xind,zind);
end

maxu = max(abs(horiz(:,2)));
normu = horiz(:,2)./maxu;

edges = [-1:0.1:1];

figure1 = figure( 'Color', [ 1 1 1 ] );
box on
h = histfit(normu,20,'normal');
h(1).FaceColor = [0, 0.4470, 0.7410];
h(1).LineWidth = 1;
h(2).LineWidth = 3;
% h = histogram(normu, edges);
% hold on

set(gca, 'FontName','times', 'FontSize', 20,'linewidth',1.5)
grid on
xlabel('Speed (n.u.)')
ylabel('No. of Cycles')
xlim([-1.5 1.5])

%% Export fig
exportgraphics(gcf,name,'resolution',600)
