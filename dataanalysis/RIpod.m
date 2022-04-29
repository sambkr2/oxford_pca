%% Pre-requisites
% Run main_SB first to get PODApproxx

%% Set crank angle
cad = -90;
cadstr = num2str(cad);
ccm_cad = strrep(cadstr,'-','m');

%% Load ccm data
load('ccm_T1_mot.mat')

%% 
relevance_index = zeros(length(nModes), 2);

for i = 1:length(nModes)
    
    temp_mode = num2str(nModes(i));
    temp_index = ['POD',temp_mode];
   
    temp_piv_u = PODVel.(temp_index).u;
    temp_piv_v = PODVel.(temp_index).v;
    temp_PIVem_SpeedMap = abs( complex( temp_piv_u, temp_piv_v ) );
    
    temp_ccm_u = ccmdata.(ccm_cad).u;
    temp_ccm_v = ccmdata.(ccm_cad).w;
    
    temp_PIV_mask = ~isnan( temp_PIVem_SpeedMap );
    temp_PIV_mask = double( temp_PIV_mask );
    temp_PIV_mask( temp_PIV_mask==0 ) = NaN;
    temp_ccm_u = temp_ccm_u .* temp_PIV_mask;
    temp_ccm_v = temp_ccm_v .* temp_PIV_mask;
    
    temp_RI_piv_ccm = Find_Relevance_Index( complex( temp_piv_u, temp_piv_v ), complex( temp_ccm_u, temp_ccm_v ), 'normal_num' );
    relevance_index(i,2) = temp_RI_piv_ccm;
    relevance_index(i,1) = nModes(i);
    
end

%% Save
save('T1rke/ri_pod_ccm_cyc96_m90.mat','relevance_index')

%% Plot
figure1 = figure( 'Color', [ 1 1 1 ] );
axes1 = axes( 'Parent', figure1, 'FontName', 'serif', 'FontSize', 14,'linewidth',1);
hold( axes1, 'on' )
box( axes1, 'on')
plot(relevance_index(:,1),relevance_index(:,2),'linewidth',1.5)

xlabel('Number of POD modes','FontSize',22)
ylabel('Relevance index [-]','FontSize',22)
% legend('RKE2L init','RKE2L def','location','northwest')
% set(axes1, 'xtick', -360:45:360)