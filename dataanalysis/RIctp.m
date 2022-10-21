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

%%
riT1rkeSie = [(cadFirst:cadStep:cadFirst+cadStep*(length(files)-1))', relevance_index];
% riT1 = [(-270:5:-30)',relevance_index(13:end)];
save('riT1rkeSie.mat', 'riT1rkeSie')

%%
first = MaskedData.CrankAngle(1);
last = MaskedData.CrankAngle(1,end);
len = length([first:5:last]);
ri_ccm = zeros(len,2);

%% d
% CCM data

for i = 1:len

    CrankAngle_plot = MaskedData.CrankAngle(1,i);
%     [ ~, CrankAngleNo_plot ] = ismember( CrankAngle_plot, MaskedData.CrankAngle );
    str_cad = num2str(CrankAngle_plot);
    ccm_cad = strrep(str_cad,'-','m');
    ccm_u = ccmdata.(ccm_cad).v;
    ccm_v = ccmdata.(ccm_cad).w;
    
    % PIV ensemble mean
    piv_u_em = mean( MaskedData.U( :,:,i,: ), 4, 'omitnan' );
    piv_v_em = mean( MaskedData.V( :,:,i,: ), 4, 'omitnan' );
        
    % temp_piv_u = MaskedData.U( :, :, CrankAngleNo_plot, k );
    % temp_piv_v = MaskedData.V( :, :, CrankAngleNo_plot, k );
    temp_PIVem_SpeedMap = abs( complex( piv_u_em, piv_v_em ) );
    
    temp_PIV_mask = ~isnan( temp_PIVem_SpeedMap );
    temp_PIV_mask = double( temp_PIV_mask );
    temp_PIV_mask( temp_PIV_mask==0 ) = NaN;
    temp_ccm_u = ccm_u .* temp_PIV_mask;
    temp_ccm_v = ccm_v .* temp_PIV_mask;
    
    temp_RI_piv_ccm = Find_Relevance_Index( complex( piv_u_em, piv_v_em ), complex( temp_ccm_u, temp_ccm_v ), 'normal_num' );
    ri_ccm(i,1) = CrankAngle_plot;
    ri_ccm(i,2) = temp_RI_piv_ccm;
    
%     str = string(sprintf('%d',k));
%     app = append('Cycle ', str, ' complete')

end
