function [ hf1, ha1, hc1 ] = Plot_ColourMapandQuiver( Figure_No, Data_x, Data_y, Data_u, Data_v, Data_mag,...
                                xPcolorShift, yPcolorShift, NorScale,...
                                XLabel_String, YLabel_String, CLabel_String, Plot_FontSize, Plot_FontName,...
                                XLim_Array, YLim_Array, Clim_max, SpacingReduceFactor, isQuiver, BWvector_threshold )

hf1 = figure( Figure_No );

% define figure background colour
hf1.Color = [ 1 1 1 ];

colormap( parula(64) );

clf
hold on
box on
x_pcolor = Data_x - xPcolorShift;
y_pcolor = Data_y - yPcolorShift;
mag_pcolor = Data_mag/NorScale;

hplots = pcolor( x_pcolor, y_pcolor, mag_pcolor );
hplots.LineStyle = 'none';

% shading interp 

ha1 = gca;
ha1.YDir = 'normal';
ha1.CLim = [ 0 Clim_max ];



hc1 = colorbar;
hc1.Label.String = CLabel_String;
% hc1.Location = 'northoutside';


axis equal
xlim( XLim_Array )
ylim( YLim_Array )


xlabel( XLabel_String )
ylabel( YLabel_String )

if isQuiver
    % create new axes for vectors
    hAxesVec = axes;
    
    % set vector axes to transparent
    set(hAxesVec,'Color','none')
    
    % set colormap for  vectors
    colormap(hAxesVec,fliplr(gray));
    
    % plot vectors using background to define contrast colour
    % (modified file exchange)
    
    my_ncquiverref_contrast(    Data_x( SpacingReduceFactor:SpacingReduceFactor:(end-SpacingReduceFactor), SpacingReduceFactor:SpacingReduceFactor:(end-SpacingReduceFactor) ),...
                                Data_y( SpacingReduceFactor:SpacingReduceFactor:(end-SpacingReduceFactor), SpacingReduceFactor:SpacingReduceFactor:(end-SpacingReduceFactor) ),...
                                Data_u( SpacingReduceFactor:SpacingReduceFactor:(end-SpacingReduceFactor), SpacingReduceFactor:SpacingReduceFactor:(end-SpacingReduceFactor) )/NorScale,...
                                Data_v( SpacingReduceFactor:SpacingReduceFactor:(end-SpacingReduceFactor), SpacingReduceFactor:SpacingReduceFactor:(end-SpacingReduceFactor) )/NorScale,...
                                '',1,0,'col',BWvector_threshold*Clim_max,...
                                Data_mag( SpacingReduceFactor:SpacingReduceFactor:(end-SpacingReduceFactor), SpacingReduceFactor:SpacingReduceFactor:(end-SpacingReduceFactor) )/NorScale,...
                                'nonquiver')
    
    caxis(hAxesVec,[0 1]);
    
    axis equal
    
    set(findall(hf1,'-property','FontSize'),'FontSize',Plot_FontSize)
    set(findall(hf1,'-property','FontName'),'FontName', Plot_FontName)
    set(findall(hf1,'-property','interpreter'),'interpreter', 'latex')
    
    % set to overlap original axes again due to colorbar change
    hAxesVec.Position = ha1.Position;
    hAxesVec.XLim = ha1.XLim;
    hAxesVec.YLim = ha1.YLim;
    hAxesVec.XTick = ha1.XTick;
    hAxesVec.YTick = ha1.YTick;
    
    linkprop([ha1,hAxesVec],{'Position','XLim','YLim','XTick','YTick','XColor','YColor'});
    
else
    set(findall(hf1,'-property','FontSize'),'FontSize',Plot_FontSize)
    set(findall(hf1,'-property','FontName'),'FontName', Plot_FontName)
    set(findall(hf1,'-property','interpreter'),'interpreter', 'latex')
end

% ha1.Position = [ 0 0 1 1 ];
% ha1.DataAspectRatio = [ 1 1 1 ];
ha1.PlotBoxAspectRatio = [ diff(XLim_Array) diff(YLim_Array) 1 ];
hf1.Position(1) = 200;
hf1.Position(2) = 200;
hf1.Position(4) = hf1.Position(3) * diff(YLim_Array) / diff(XLim_Array);
% 
% ha1.XColor = 'none';
% ha1.YColor = 'none';


end