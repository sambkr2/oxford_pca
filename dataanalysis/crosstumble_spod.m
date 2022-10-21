%% Load
load('../dataformatting/x20190517_CrossTumble_CR12p5_T1_C33_DVA_Motored_Processed_all.mat')

%% Get original data
CrankAngle_plot = -250;
[ ~, CrankAngleNo_plot ] = ismember( CrankAngle_plot, MaskedData.CrankAngle );

% Get spatial vectors
x = MaskedData.X;
y = MaskedData.Y;

% Get velocity data at the specified crank angle
u = squeeze(MaskedData.U(:,:,CrankAngleNo_plot,:));
w = squeeze(MaskedData.V(:,:,CrankAngleNo_plot,:));

%% Prepare data for SPOD use
nCycle = size(u,3);

% Need to find ones
originalIndex = find(~isnan(MaskedData.U(:,:,CrankAngleNo_plot,1)) == 1);
nValidPoints = length(originalIndex);

% Take valid points and put them into a column
SPODData.X = x(originalIndex);
SPODData.Y = y(originalIndex);

% Initialise velocity data
SPODData.U = zeros( nValidPoints, nCycle );
SPODData.V = zeros( nValidPoints, nCycle );

% For each cycle, fill in the SPOData velocities
for cycle = 1:nCycle
    temp_u = u(:,:,cycle);
    temp_w = w(:,:,cycle);
    SPODData.U( :,cycle ) = temp_u(originalIndex);
    SPODData.V( :,cycle ) = temp_w(originalIndex);
end

%% SPOD
complex_velo = [complex( SPODData.U, SPODData.V )];
velo = [real(complex_velo); imag(complex_velo)];
veloT = velo';

% rpm
rpm = 1500;

% time between snapshots. *2 factor accounts for 4 stroke engine
dt = 1/(rpm/60)*2;

% SPOD
[L2,P2,f2] = spod(veloT,[],[],[],dt);

%% Reconstruction
% Number of velocity locations
nVel = size(velo,1)/2;

% First choose the SPOD modes and frequency to plot
sMode = [1];
    
for fplot = 1:9
    %
    data = P2(fplot,:,sMode)*L2(fplot,sMode);
    data = data';
    
    % Retrieve velocities
    spodU = real(data(1:nVel));
    spodV = real(data(nVel+1:end));
    
    % Reshape into PIV grid format
    pivy = size(CFDData.Data.y_PIVGrid,1);
    pivz = size(CFDData.Data.z_PIVGrid,2);
    splotu = NaN(pivy,pivz);
    splotv = NaN(pivy,pivz);
    splotu(originalIndex) = spodU;
    splotv(originalIndex) = spodV;
    
    % plot
    figureprop.axes_lim = [ -25 25 -20 2 ];
    figureprop.xlabel = '{\it x} (mm)';
    figureprop.ylabel = '{\it z} (mm)';
    % figureprop.velocity_normalisation = 5;
    figureprop.sparse_vector = 2;
    figureprop.Clim = [ 0 50 ];
    
    figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, splotu*1000, splotv*1000, figureprop );    
    t = {['SPOD mode 1', num2str(f2(fplot)), ' Hz']};
    title(t)
    set(gcf,'position',[440 500 560 280])
end
