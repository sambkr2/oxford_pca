%% Pre-requisites
% Run main_SB first to get PODApprox

%% Load cross tumble
load('../dataformatting/x20190517_CrossTumble_CR12p5_T1_C33_DVA_Motored_Processed_all.mat')
load('/Users/sambaker/Documents/Oxford-Uni/Technical/PCA/dataanalysis/rrT1rng/ccm_T1_mot_CTP.mat')

%% Set Parameters
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

%% POD approx paramters
nModes = [ 0:299 ];
% CycleNo = [ 1:300 ];
% nModes = [ 0 1 2 5 10 50 299 ];
CycleNo = [35 36];
cycstr = num2str(CycleNo);

%% POD approx
for m = 1:length(CycleNo)
    for mm = 1 : length( nModes )
        [ PODApprox ] = Calc_PODApprox( PODResult{1,1}, nModes(mm), CycleNo(m) );
        PODVel.(['Cycle',num2str(CycleNo(m))]).(['POD',num2str(nModes(mm))]).u = PODApprox.U;
        PODVel.(['Cycle',num2str(CycleNo(m))]).(['POD',num2str(nModes(mm))]).v = PODApprox.V;
    end
    disp(['Cycle ' num2str(CycleNo(m)) ' complete'])
end

% %%
% figure_output = ColourQuiver_SB( PODApprox.X, PODApprox.Y, temp_piv_u, temp_piv_v, figureprop );
% ylim([-30 10])
% 
% %%
% figure_output = ColourQuiver_SB( PODApprox.X, PODApprox.Y, temp_ccm_u, temp_ccm_v, figureprop );
% ylim([-30 10])

%% Set ccm crank angle
cad = AnalysisResult.CrankAngle;
cadstr = num2str(cad);
ccm_cad = strrep(cadstr,'-','m');

%% One cycle
ri = zeros(length(nModes), 2);

for i = 1:length(nModes)
    
    temp_mode = num2str(nModes(i));
    temp_index = ['POD',temp_mode];
    cyc_ind = ['Cycle',cycstr];
   
    temp_piv_u = PODVel.(cyc_ind).(temp_index).u;
    temp_piv_v = PODVel.(cyc_ind).(temp_index).v;
    temp_PIVem_SpeedMap = abs( complex( temp_piv_u, temp_piv_v ) );
    
    temp_ccm_u = ccmdata.(ccm_cad).v;
    temp_ccm_v = ccmdata.(ccm_cad).w;
    
    temp_PIV_mask = ~isnan( temp_PIVem_SpeedMap );
    temp_PIV_mask = double( temp_PIV_mask );
    temp_PIV_mask( temp_PIV_mask==0 ) = NaN;
    temp_ccm_u = temp_ccm_u .* temp_PIV_mask;
    temp_ccm_v = temp_ccm_v .* temp_PIV_mask;
    
    temp_RI_piv_ccm = Find_Relevance_Index( complex( temp_piv_u, temp_piv_v ), complex( temp_ccm_u, temp_ccm_v ), 'normal_num' );
    ri(i,2) = temp_RI_piv_ccm;
    ri(i,1) = nModes(i);
    
end

%% All cycles

% take CFD velocities
temp_ccm_u = ccmdata.(ccm_cad).v;
temp_ccm_v = ccmdata.(ccm_cad).w;

% for cycle = 1:300
for cycle = 35:36
    
    a = num2str(cycle);
    b = ['cycle' a];
    ri2.(b) = zeros(length(nModes), 2);
    
    temp_cycle = ['Cycle',num2str(cycle)];
    
    for i = 1:length(nModes)
        
        % which mode are we using
        temp_mode = num2str(nModes(i));
        temp_index = ['POD',temp_mode];
        
        % take PIV velocities at this mode
        temp_piv_u = PODVel.(temp_cycle).(temp_index).u;
        temp_piv_v = PODVel.(temp_cycle).(temp_index).v;
        temp_PIVem_SpeedMap = abs( complex( temp_piv_u, temp_piv_v ) );
    
        % mask CFD data
        temp_PIV_mask = ~isnan( temp_PIVem_SpeedMap );
        temp_PIV_mask = double( temp_PIV_mask );
        temp_PIV_mask( temp_PIV_mask==0 ) = NaN;
        temp_ccm_u = temp_ccm_u .* temp_PIV_mask;
        temp_ccm_v = temp_ccm_v .* temp_PIV_mask;
        
        % find RI between masked CFD data and PIV at this mode
        temp_RI_piv_ccm = Find_Relevance_Index( complex( temp_piv_u, temp_piv_v ), complex( temp_ccm_u, temp_ccm_v ), 'normal_num' );
        ri2.(b)(i,2) = temp_RI_piv_ccm;
        ri2.(b)(i,1) = nModes(i);
    end
    
    disp(['Cycle ' a ' complete'])
    
end

%% Plot all cycles
figure3 = figure( 'Color', [ 1 1 1 ] );
axes3 = axes( 'Parent', figure3, 'FontSize', 24,'linewidth',1);
hold( axes3, 'on' )
box( axes3, 'on')
for plotcyc = 1:300
% for plotcyc = 35:36
    
    c = num2str(plotcyc);
    d = ['cycle' c];
    plot(ri2.(d)(:,1),ri2.(d)(:,2),'linewidth',1.5)
    
end
xlim([0 10])
% ylim([0.4 0.9])
xlabel('Number of POD modes')
ylabel('Relevance index')
% t = {['Cross-tumble CFD against PIV POD modes'], ['all cycles, -285 CAD zoom']};
% title(t,'FontSize',20)
fontname(figure3,"times")
% legend('Cycle A', 'Cycle B')
% grid on
% set(gca, 'xtick', 0:50:300)
    
%% Save
% save('rrT1rng/ri_pod_ccm_all_m270.mat','ri')
pivname = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/CTP/CTP_allcycles_RI_zoom.png'];
exportgraphics(gcf,pivname,'resolution',600)

%% Plot
figure3 = figure( 'Color', [ 1 1 1 ] );
axes3 = axes( 'Parent', figure3, 'FontName', 'times', 'FontSize', 24,'linewidth',1);
hold( axes3, 'on' )
box( axes3, 'on')
plot(ri2.cycle1(:,1),ri2.cycle1(:,2),'linewidth',1.5)

xlabel('Number of POD modes')
ylabel('Relevance index')
% legend('RKE2L init','RKE2L def','location','northwest')
% set(axes1, 'xtick', -360:45:360)
cycstr = num2str(CycleNo);
% title(['Cycle ',cycstr,', ',cadstr,' CAD'])

%%
% figure1 = figure( 'Color', [ 1 1 1 ] );
% axes1 = axes( 'Parent', figure1, 'FontName', 'serif', 'FontSize', 14,'linewidth',1);
% hold( axes1, 'on' )
% box( axes1, 'on')
% for k = 1:length(CycleNo)
%     hold on
%     c = num2str(k);
%     d = ['cycle' c];
%     plot(ri.(d)(:,1),ri.(d)(:,2))
% end