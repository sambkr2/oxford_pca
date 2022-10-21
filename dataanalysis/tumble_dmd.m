%% Load
load('../dataformatting/x20180706_Tumble_CR12p5_T1_C33_DVA_Motored_Processed_all_masked.mat')

%% Get variables
CrankAngle_plot = -285;
[ ~, CrankAngleNo_plot ] = ismember( CrankAngle_plot, MaskedData.CrankAngle );

% Get spatial vectors
x = MaskedData.X;
y = MaskedData.Y;

% Get velocity data at the specified crank angle
u = squeeze(MaskedData.U(:,:,CrankAngleNo_plot,:));
w = squeeze(MaskedData.V(:,:,CrankAngleNo_plot,:));

%% Prepare data for DMD use
nCycle = size(u,3);

% Need to find ones
originalIndex = find(~isnan(MaskedData.U(:,:,CrankAngleNo_plot,1)) == 1);
nValidPoints = length(originalIndex);

% Take valid points and put them into a column
DMDData.X = x(originalIndex);
DMDData.Y = y(originalIndex);

% Initialise velocity data
DMDData.U = zeros( nValidPoints, nCycle );
DMDData.V = zeros( nValidPoints, nCycle );

% For each cycle, fill in the DMDData velocities
for cycle = 1:nCycle
    temp_u = u(:,:,cycle);
    temp_w = w(:,:,cycle);
    DMDData.U( :,cycle ) = temp_u(originalIndex);
    DMDData.V( :,cycle ) = temp_w(originalIndex);
end

%% DMD
velo = [ DMDData.U; DMDData.V ];
X1 = velo(:,1:end-1);
X2 = velo(:,2:end);

% Choose cut-off
% r = 299;
r = GDmode;

% DMD
[P2, L2, b2] = DMD(X1,X2,r);

%% Plot DMD spectrum
% rpm
rpm = 1500;

% time between snapshots. *2 factor accounts for 4 stroke engine (seconds
% per cycle)
dt = 1/(rpm/60)*2;

% get DMD frequency
dmdFreq = log(diag(L2))/dt/2/pi;

% get amplitude
amp = abs(b2(:));

% scatter raw
% figure
% semilogy(imag(dmdFreq), amp,'o','MarkerFaceColor','blue')
% xlim([0 1])

% scatter neat
figure1 = figure( 'Color', [ 1 1 1 ] );
axes1 = axes( 'Parent', figure1, 'FontName', 'Times', 'FontSize', 20,'linewidth',1);
hold( axes1, 'on' )
box( axes1, 'on')
semilogy(imag(dmdFreq), amp,'o','MarkerFaceColor','blue')
xlim([0 6.5])
xlabel('Frequency (Hz)')
ylabel('DMD mode amplitude')
title({'De-noised spectrum:','tumble plane, -285 CAD'})
% title({'Raw spectrum:','tumble plane, -285 CAD'})

%% export fig
name = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/TP/tp_dmdspectrum_clean_',ccm_cad,'.png'];
% name = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/TP/tp_dmdspectrum_raw_',ccm_cad,'.png'];
exportgraphics(gcf,name,'resolution',600)

%% Plot settings
pivgrid = CFDData.Data.x_PIVGrid;

figureprop.axes_lim = [ -25 25 -30 10 ];
figureprop.xlabel = '{\it x} (mm)';
figureprop.ylabel = '{\it z} (mm)';
figureprop.velocity_normalisation = 5;
figureprop.sparse_vector = 2;
figureprop.Clim = [ 0 50 ];

%% Plot first mode
index = 1;
first = P2(:,index)*L2(index,index)*b2(index,index);

% Separate U and V velocities
Ucol = real(first(1:nValidPoints));
Vcol = real(first(nValidPoints+1:end));

% Reshape into PIV grid format
pivx = size(pivgrid,1);
pivz = size(pivgrid,2);
dmdplotu = NaN(pivx,pivz);
dmdplotv = NaN(pivx,pivz);
dmdplotu(originalIndex) = Ucol;
dmdplotv(originalIndex) = Vcol;

figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, dmdplotu, dmdplotv, figureprop );    
ylim([-30 10])
t = {['DMD mode ', num2str(index), ', 0 Hz']};
title(t)

%% Plot full matrix
data = P2*L2*b2;

% Separate U and V velocities
locs = size(data);
Ucol = real(data(1:locs/2));
Vcol = real(data(locs/2+1:end));

% Reshape into PIV grid format
pivx = size(pivgrid,1);
pivz = size(pivgrid,2);
dmdplotu = NaN(pivx,pivz);
dmdplotv = NaN(pivx,pivz);
dmdplotu(originalIndex) = Ucol;
dmdplotv(originalIndex) = Vcol;

figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, dmdplotu, dmdplotv, figureprop );    
ylim([-30 10])
t = {['All modes reconstructed']};
title(t)

%% Find largest amplitudes
% sort amplitudes in size order
dvals = sort(amp, 'descend');

% number of values to look at
nvals = 15;
largest = dvals(1:nvals);

% find index in amp/b2 where the largest values are
ind = zeros(nvals,1);
for k = 1:nvals
    ind(k,1) = find(amp == largest(k));
end

% attach frequency to amp inde
ampf = [ind, dmdFreq(ind)];

% find index of positive frequencies
posampf_ind = imag(ampf(:,2)) >= 0;

% retain positive frequency indices
posf_ind = ampf(posampf_ind);

%% Inspect modes at highest amplitudes

for kk = 1:length(posf_ind)
    
    % DMD mode * DMD eigenvalue * amplitude
    recon = P2(:,posf_ind(kk))*L2(posf_ind(kk),posf_ind(kk))*b2(posf_ind(kk),1);
    
    % Separate U and V velocities
    locs = size(recon);
    Ucol = real(recon(1:locs/2));
    Vcol = real(recon(locs/2+1:end));
    
    % Reshape into PIV grid format
    pivx = size(pivgrid,1);
    pivz = size(pivgrid,2);
    dmdplotu = NaN(pivx,pivz);
    dmdplotv = NaN(pivx,pivz);
    dmdplotu(originalIndex) = Ucol;
    dmdplotv(originalIndex) = Vcol;
    
    figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, dmdplotu, dmdplotv, figureprop );    
    ylim([-30 10])
    t = {['DMD mode ', num2str(posf_ind(kk)), ', ', num2str(round(imag(dmdFreq(posf_ind(kk))))),' Hz']};
    title(t)

end

%% WARNING
% Careful not to overwrite images!!

%% Inspect at 0 Hz freq
% near-zero frequencies
zeroind = imag(dmdFreq) <= 0.1 & imag(dmdFreq) >= 0;
zeroindf = find(zeroind == 1);

for jj = 1:length(zeroindf)
    
    % DMD mode * DMD eigenvalue * amplitude
%     recon = P2(:,posf_ind(kk))*L2(posf_ind(kk),posf_ind(kk))*b2(posf_ind(kk),1);
%     recon = P2(:,zeroindf(jj))*L2(zeroindf(jj),zeroindf(jj))*amp(zeroindf(jj),1);
    recon = P2(:,zeroindf(jj))*-amp(zeroindf(jj),1);
    
    % Separate U and V velocities
    locs = size(recon);
    Ucol = real(recon(1:locs/2));
    Vcol = real(recon(locs/2+1:end));
    
    % Reshape into PIV grid format
    pivx = size(pivgrid,1);
    pivz = size(pivgrid,2);
    dmdplotu = NaN(pivx,pivz);
    dmdplotv = NaN(pivx,pivz);
    dmdplotu(originalIndex) = Ucol;
    dmdplotv(originalIndex) = Vcol;
    
    figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, dmdplotu, dmdplotv, figureprop );    
%     t = {['DMD mode ', num2str(zeroindf(jj)), ', ', num2str(round(imag(dmdFreq(zeroindf(jj))))),' Hz']};
%     title(t)

%     set(gcf, 'position', [440 377 560 290])
    dmdname = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/TP/DMD',num2str(zeroindf(jj)),'.png'];
    exportgraphics(gcf,dmdname,'resolution',600)

end

%% Trying different r values
r11 = 5;
[P11, L11, b11] = DMD(X1,X2,r11);
data11 = P11*L11*b11;

% Separate U and V velocities
locs = size(data11);
Ucol = real(data11(1:locs/2));
Vcol = real(data11(locs/2+1:end));

% Reshape into PIV grid format
pivx = size(pivgrid,1);
pivz = size(pivgrid,2);
dmdplotu = NaN(pivx,pivz);
dmdplotv = NaN(pivx,pivz);
dmdplotu(originalIndex) = Ucol;
dmdplotv(originalIndex) = Vcol;

figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, dmdplotu, dmdplotv, figureprop );    
ylim([-30 10])
t = {['All modes reconstructed, r = ',num2str(r11)]};
title(t)

%%
data111 = P11(:,177)*L11(177,177)*b11(177,1);

% Separate U and V velocities
locs = size(data111);
Ucol = real(data111(1:locs/2));
Vcol = real(data111(locs/2+1:end));

% Reshape into PIV grid format
pivx = size(pivgrid,1);
pivz = size(pivgrid,2);
dmdplotu = NaN(pivx,pivz);
dmdplotv = NaN(pivx,pivz);
dmdplotu(originalIndex) = Ucol;
dmdplotv(originalIndex) = Vcol;

figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, dmdplotu, dmdplotv, figureprop );    
ylim([-30 10])
t = {['All modes reconstructed, cutoff 11']};
title(t)

%% Plot unit circle
figure
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta))
hold on
grid on
scatter(real(diag(L2)),imag(diag(L2)))
axis([-1.1 1.1 -1.1 1.1])

%% Plotting real and imaginary components of DMD modes separately
dmdmode = P2(:,292);

% Separate U and V velocities
locs = size(dmdmode);
Ucol = dmdmode(1:locs/2);
Vcol = dmdmode(locs/2+1:end);

% Reshape into PIV grid format
pivx = size(pivgrid,1);
pivz = size(pivgrid,2);
dmdplotu = complex(NaN(pivx,pivz),NaN(pivx,pivz));
dmdplotv = complex(NaN(pivx,pivz),NaN(pivx,pivz));
dmdplotu(originalIndex) = Ucol;
dmdplotv(originalIndex) = Vcol;

figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, real(dmdplotu)*1000, real(dmdplotv)*1000, figureprop );
figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, imag(dmdplotu)*1000, imag(dmdplotv)*1000, figureprop );

