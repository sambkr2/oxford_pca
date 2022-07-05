%% Pre-requisites
% Run main_CTP_analysis first to get PODVel

%% Set crank angle
cad = -285;
cadstr = num2str(cad);
ccm_cad = strrep(cadstr,'-','m');

%% Load ccm data
% load('ccm_T1_mot.mat')
load('rrT1rng/ccm_T1_mot_CTP.mat')

%% One cycle
ri = zeros(length(nModes), 2);

for i = 1:length(nModes)
    
    temp_mode = num2str(nModes(i));
    temp_index = ['POD',temp_mode];
   
    temp_piv_u = PODVel.(temp_index).u;
    temp_piv_v = PODVel.(temp_index).v;
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

for cycle = 1:1:length(CycleNo)
    
    a = num2str(cycle);
    b = ['cycle' a];
    ri.(b) = zeros(length(nModes), 2);
    
    for i = 1:length(nModes)
        
        % which mode are we using
        temp_mode = num2str(nModes(i));
        temp_index = ['POD',temp_mode];
        
        % take PIV velocities at this mode
        temp_piv_u = PODVel.(temp_index).u;
        temp_piv_v = PODVel.(temp_index).v;
        temp_PIVem_SpeedMap = abs( complex( temp_piv_u, temp_piv_v ) );
        
        % take CFD velocities
        temp_ccm_u = ccmdata.(ccm_cad).u;
        temp_ccm_v = ccmdata.(ccm_cad).w;
    
        % mask CFD data
        temp_PIV_mask = ~isnan( temp_PIVem_SpeedMap );
        temp_PIV_mask = double( temp_PIV_mask );
        temp_PIV_mask( temp_PIV_mask==0 ) = NaN;
        temp_ccm_u = temp_ccm_u .* temp_PIV_mask;
        temp_ccm_v = temp_ccm_v .* temp_PIV_mask;
        
        % find RI between masked CFD data and PIV at this mode
        temp_RI_piv_ccm = Find_Relevance_Index( complex( temp_piv_u, temp_piv_v ), complex( temp_ccm_u, temp_ccm_v ), 'normal_num' );
        ri.(b)(i,2) = temp_RI_piv_ccm;
        ri.(b)(i,1) = nModes(i);
    end
    disp(['Cycle ' a ' complete'])
end

%% Save
save('rrT1rng/ri_pod_ccm_all_m270.mat','ri')

%% Plot
figure1 = figure( 'Color', [ 1 1 1 ] );
axes1 = axes( 'Parent', figure1, 'FontName', 'serif', 'FontSize', 14,'linewidth',1);
hold( axes1, 'on' )
box( axes1, 'on')
plot(ri(:,1),ri(:,2),'linewidth',1.5)

xlabel('Number of POD modes','FontSize',22)
ylabel('Relevance index [-]','FontSize',22)
% legend('RKE2L init','RKE2L def','location','northwest')
% set(axes1, 'xtick', -360:45:360)
cycstr = num2str(CycleNo);
title(['RANS vs POD modes, ',cadstr,' CAD'])

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