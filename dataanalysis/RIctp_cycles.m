%%
% Run main_CTP_analysis first to get PODVel

% load ccm
load('rrT1rng/ccm_T1_mot_CTP.mat')

%%
% CCM data
ri_ccm = zeros(300, 2);
CrankAngle_plot = -285;
[ ~, CrankAngleNo_plot ] = ismember( CrankAngle_plot, MaskedData.CrankAngle );
str_cad = num2str(CrankAngle_plot);
ccm_cad = strrep(str_cad,'-','m');
ccm_u = ccmdata.(ccm_cad).v;
ccm_v = ccmdata.(ccm_cad).w;

% PIV ensemble mean
piv_u_em = nanmean( MaskedData.U( :,:,CrankAngleNo_plot,: ), 4 );
piv_v_em = nanmean( MaskedData.V( :,:,CrankAngleNo_plot,: ), 4 );

for k = 1:300
    
    temp_piv_u = MaskedData.U( :, :, CrankAngleNo_plot, k );
    temp_piv_v = MaskedData.V( :, :, CrankAngleNo_plot, k );
    temp_PIVem_SpeedMap = abs( complex( temp_piv_u, temp_piv_v ) );

%     temp_mode = num2str(k-1);
%     temp_index = ['POD',temp_mode];
% 
%     % take PIV velocities at this mode
%     temp_piv_u = PODVel.(temp_index).u;
%     temp_piv_v = PODVel.(temp_index).v;
%     temp_PIVem_SpeedMap = abs( complex( temp_piv_u, temp_piv_v ) );

    temp_PIV_mask = ~isnan( temp_PIVem_SpeedMap );
    temp_PIV_mask = double( temp_PIV_mask );
    temp_PIV_mask( temp_PIV_mask==0 ) = NaN;
    temp_ccm_u = ccm_u .* temp_PIV_mask;
    temp_ccm_v = ccm_v .* temp_PIV_mask;

    temp_RI_piv_ccm = Find_Relevance_Index( complex( temp_piv_u, temp_piv_v ), complex( temp_ccm_u, temp_ccm_v ), 'normal_num' );
    ri_ccm(k,1) = k;
    ri_ccm(k,2) = temp_RI_piv_ccm;
    
    temp_RI_piv_em = Find_Relevance_Index( complex( temp_piv_u, temp_piv_v ), complex( piv_u_em, piv_v_em ), 'normal_num' );
    ri_em(k,1) = k;
    ri_em(k,2) = temp_RI_piv_em;
    
    str = string(sprintf('%d',k));
    app = append('Cycle ', str, ' complete')
    
end

%% Plot
figure
axes('fontname','serif','fontsize',20,'linewidth',1)
hold on
scatter(ri_ccm(:,1),ri_ccm(:,2),'filled')
scatter(ri_em(:,1),ri_em(:,2),'filled')
% yline(0.7)
box on
l = legend('ccm vs ind cycle','em vs ind cycle');
% set(l,'location','southwest')
title([str_cad,' CAD'])
xlabel('Cycle Number')
ylabel('Relevance Index')
ylim([0 1])

%% 
% row counter
r = 0;
for i = 1:300
    if ri_ccm(i,2) > ri_em(i,2)
      r = r+1;
      ri_increase(r,1) = i;
      increase = ri_ccm(i,2) - ri_em(i,2);
      ri_increase(r,2) = increase;
    end
end

%% Plot
figure
axes('fontname','serif','fontsize',20,'linewidth',1)
hold on
scatter(ri_increase(:,1),ri_increase(:,2),'filled')
% scatter(ri_em(:,1),ri_em(:,2),'filled')
% yline(0.7)
box on
% l = legend('ccm vs ind cycle','em vs ind cycle');
% set(l,'location','southwest')
title([str_cad,' CAD'])
xlabel('Cycle Number')
ylabel('CCM and PIV RI diff')
ylim([0 0.25])