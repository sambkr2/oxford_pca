%% RI calcs
% RI.piv_cd = Find_Relevance_Index( complex( temp_PIVem_u, temp_PIVem_v ), complex( temp_CFD_u, temp_CFD_v ), 'normal_num' );
relevance_index = zeros(length(files), 1);

for i = cadFirst:cadStep:cadFirst+cadStep*(length(files)-1)
    
    [ ~, PIV_CAindex ] = ismember( i, PIVData.Data.CrankAngle );
    temp_piv_u = nanmean( PIVData.Data.u( :,:,PIV_CAindex,: ), 4 );
    temp_piv_v = nanmean( PIVData.Data.w( :,:,PIV_CAindex,: ), 4 );
    temp_PIVem_SpeedMap = abs( complex( temp_piv_u, temp_piv_v ) );

    temp_ccm_num = num2str(i);
    temp_ccm_cad = strrep(temp_ccm_num,'-','m');
    temp_ccm_u = ccmdata.(temp_ccm_cad).u;
    temp_ccm_v = ccmdata.(temp_ccm_cad).w;

    temp_PIV_mask = ~isnan( temp_PIVem_SpeedMap );
    temp_PIV_mask = double( temp_PIV_mask );
    temp_PIV_mask( temp_PIV_mask==0 ) = NaN;
    temp_ccm_u = temp_ccm_u .* temp_PIV_mask;
    temp_ccm_v = temp_ccm_v .* temp_PIV_mask;

    temp_RI_piv_ccm = Find_Relevance_Index( complex( temp_piv_u, temp_piv_v ), complex( temp_ccm_u, temp_ccm_v ), 'normal_num' );
    relevance_index(PIV_CAindex) = temp_RI_piv_ccm;
    
%     str = string(sprintf('%d',i));
%     app = append('Crank angle ', str, ' complete')
    
end

riT1rkeSie = [(cadFirst:cadStep:cadFirst+cadStep*(length(files)-1))', relevance_index];
% riT1 = [(-270:5:-30)',relevance_index(13:end)];
save('riT1rkeSie.mat', 'riT1rkeSie')