figureprop.Clim = [ 0 50 ];

ColourQuiver_SB( PODApprox.X, PODApprox.Y, PODApprox.U, PODApprox.V, figureprop );
title('PODApprox U and V')

%%
ColourQuiver_SB( PODApprox.X, PODApprox.Y, ccm_u, ccm_v, figureprop );
title('ccm_u and v')

%%
ColourQuiver_SB( PODApprox.X, PODApprox.Y, temp_ccm_u, temp_ccm_v, figureprop );
title('masked ccm_u and v')

%%
ColourQuiver_SB( PODApprox.X, PODApprox.Y, temp_piv_u, temp_piv_v, figureprop );
title('masked data cycle 300')

%% plot A
ColourQuiver_SB( MaskedData.X, MaskedData.Y, MaskedData.U( :, :, 11, 56 ), MaskedData.V( :, :, 11, 56 ), figureprop );
title('piv cycle 56')
Find_Relevance_Index( complex( MaskedData.U( :, :, 11, 56 ), MaskedData.V( :, :, 11, 56 ) ), complex( temp_ccm_u, temp_ccm_v ), 'normal_num' )

temp_ccm_u = ccmdata.(ccm_cad).v;
temp_ccm_v = ccmdata.(ccm_cad).w;

temp_PIV_mask = ~isnan( temp_PIVem_SpeedMap );
temp_PIV_mask = double( temp_PIV_mask );
temp_PIV_mask( temp_PIV_mask==0 ) = NaN;
temp_ccm_u = temp_ccm_u .* temp_PIV_mask;
temp_ccm_v = temp_ccm_v .* temp_PIV_mask;

Find_Relevance_Index( complex( MaskedData.U( :, :, 11, 56 ), MaskedData.V( :, :, 11, 56 ) ), complex( temp_ccm_u, temp_ccm_v ), 'normal_num' );

% trusted

%%
ColourQuiver_SB( PODApprox.X, PODApprox.Y, piv_u_em, piv_v_em, figureprop );
% trusted

%%

%% plot B
ColourQuiver_SB( MaskedData.X, MaskedData.Y, PODVel.(temp_index).u, PODVel.(temp_index).v, figureprop );
title('PODVel POD299')
Find_Relevance_Index( complex( PODVel.(temp_index).u, PODVel.(temp_index).v ), complex( temp_ccm_u, temp_ccm_v ), 'normal_num' )

Find_Relevance_Index( complex( MaskedData.U( :, :, 11, 56 ), MaskedData.V( :, :, 11, 56 ) ), complex( temp_ccm_u, temp_ccm_v ), 'normal_num' )


% plot looks the same as plot A. Therefore error in find_relevance_index
% code in RIctpPOD.

% Is temp_ccm_u the same? Might be an issue with what I'm calling
% temp_ccm_u across different codes.
