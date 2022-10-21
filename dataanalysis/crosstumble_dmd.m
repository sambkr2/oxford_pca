%% Load
load('../dataformatting/x20190517_CrossTumble_CR12p5_T1_C33_DVA_Motored_Processed_all.mat')

%% Get variables
CrankAngle_plot = -285;
[ ~, CrankAngleNo_plot ] = ismember( CrankAngle_plot, MaskedData.CrankAngle );

% Get spatial vectors
x = MaskedData.X;
y = MaskedData.Y;

% Get velocity data at the specified crank angle
u = squeeze(MaskedData.U(:,:,CrankAngleNo_plot,:));
w = squeeze(MaskedData.V(:,:,CrankAngleNo_plot,:));

t_cad = num2str(CrankAngle_plot);
ccm_cad = strrep(t_cad,'-','m');

%% Prepare data for DMD use
% nCycle = size(u,3);
nCycle = 100;

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
counter = 1;
for cycle = 201:300
    temp_u = u(:,:,cycle);
    temp_w = w(:,:,cycle);
    DMDData.U( :,counter ) = temp_u(originalIndex);
    DMDData.V( :,counter ) = temp_w(originalIndex);
    counter  = counter +1;
end

%% SVD to find GD cut-off
temp_velo_data = complex(DMDData.U, DMDData.V);
temp_PODResult = Perform_POD_SB( temp_velo_data, 'Centered', 'Direct' );
svs = diag(temp_PODResult.svdS);
beta = size(u,2)/size(u,1);
tau = optimal_SVHT_coef(beta,0) * median(svs); % find cut-off tau
modes = svs(svs>tau);
GDmode = length(modes); % Gavish Donoho threshold mode

%% DMD
velo = [ DMDData.U; DMDData.V ];
X1 = velo(:,1:end-1);
X2 = velo(:,2:end);

% Choose cut-off
% r = 299;
% r = 50;
r = GDmode;
% r = 99;

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
ampr = abs(real(b2(:)));

% % scatter raw
% figure
% semilogy(imag(dmdFreq), amp,'o','MarkerFaceColor','blue')
% xlim([0 1])

% scatter neat
figure1 = figure( 'Color', [ 1 1 1 ] );
axes1 = axes( 'Parent', figure1, 'FontName', 'Times', 'FontSize', 24,'linewidth',1.5);
hold( axes1, 'on' )
box( axes1, 'on')
grid on
semilogy(imag(dmdFreq), ampr,'o','MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410],'MarkerSize',8)
% semilogy(imag(dmdFreq(11)), ampr(11),'o','MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerEdgeColor',[0.9290 0.6940 0.1250],'MarkerSize',8)
% semilogy(imag(dmdFreq(14)), ampr(14),'o','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize',8)
% semilogy(imag(dmdFreq(292)), ampr(292),'o','MarkerFaceColor',[0.494, 0.1840, 0.5560],'MarkerEdgeColor',[0.494, 0.1840, 0.5560],'MarkerSize',8)
% semilogy(imag(dmdFreq(298)), ampr(298),'o','MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerEdgeColor',[0.4660 0.6740 0.1880],'MarkerSize',8)
xlim([0 6.5])
xticks([0:1:6])
xlabel('Frequency (Hz)')
ylabel('DMD mode amplitude')
% title({'Raw spectrum:','cross-tumble plane, -285 CAD'})
% title({'De-noised spectrum:','cross-tumble plane, -285 CAD'})

%% export fig
name = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/CTP/ctp_dmdspectrum_raw_',ccm_cad,'.png'];
% name = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/CTP/ctp_dmdspectrum_clean_',ccm_cad,'.png'];
exportgraphics(gcf,name,'resolution',600)

%% plot parameters
pivgrid = CFDData.Data.y_PIVGrid;

figureprop.axes_lim = [ -25 25 -20 2 ];
figureprop.xlabel = '{\it y} (mm)';
figureprop.ylabel = '{\it z} (mm)';
% figureprop.ylabel = [];
figureprop.velocity_normalisation = 5;
figureprop.sparse_vector = 2;
figureprop.Clim = [ 0 50 ];

zeroind = imag(dmdFreq) <= 0.2 & imag(dmdFreq) >= 0;
zeroindf = find(zeroind == 1);

%% Plot first DMD mode

index = 8;
% first = P2(:,index)*L2(index,index)*b2(index,1);
% first = P2(:,index)*L2(index,index)*-amp(index,1);
first = P2(:,index)*ampr(index,1);

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
% ylim([-30 10])
% t = {['DMD mode ', num2str(index), ', 0 Hz']};
% title(t)
% set(gcf, 'position', [440 377 560 290])
% set(gcf, 'units','centimeters')
set(gcf, 'position', [357 379 762 497])
set(gcf,'OuterPosition',[357 379 762 497])
set(gca,'Position',[0.1225 0.1100 0.7302 0.8150])
% set(colorbar,'position' ,[0.8635 0.2244 0.0210 0.4858])
% Position: [444 507 560 331]


%% Plot specific DMD modes

% index = [11 12];
% index = [292];
% first = P2(:,index)*L2(index,index)*b2(index,1);
% first = P2(:,index)*L2(index,index)*-amp(index,1);
% first = P2(:,index)*amp(index,1);
% first = P2(:,index)*ampr(index,1);

% eleven and twelve
el = P2(:,11)*ampr(11,1);
tw = P2(:,14)*ampr(14,1);
av = (el+tw)/2;
first = av;

% Separate U and V velocities (real)
Ucol = real(first(1:nValidPoints));
Vcol = real(first(nValidPoints+1:end));

% Separate U and V velocities (imag)
Ucoli = imag(first(1:nValidPoints));
Vcoli = imag(first(nValidPoints+1:end));

% Reshape into PIV grid format (real)
pivx = size(pivgrid,1);
pivz = size(pivgrid,2);
dmdplotu = NaN(pivx,pivz);
dmdplotv = NaN(pivx,pivz);
dmdplotu(originalIndex) = Ucol;
dmdplotv(originalIndex) = Vcol;

% Reshape into PIV grid format (imag)
dmdplotui = NaN(pivx,pivz);
dmdplotvi = NaN(pivx,pivz);
dmdplotui(originalIndex) = Ucoli;
dmdplotvi(originalIndex) = Vcoli;

% For combined modes 11 and 12 plot
figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, dmdplotu, dmdplotv, figureprop );    

% For modes 296 ad 292 individually
% figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, dmdplotu, dmdplotv, figureprop );    

% ylim([-30 10])
% t = {['DMD mode ', num2str(index), ', 0 Hz']};
% title(t)
% set(gcf, 'position', [440 377 560 290])
% set(gcf, 'units','centimeters')
% set(gcf,'OuterPosition',[1 1 20 14])
set(gcf, 'position', [357 379 762 497])
set(gcf,'OuterPosition',[357 379 762 497])
set(gca,'Position',[0.1225 0.1100 0.7302 0.8150])

%%
dmdname = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/CTP/DMD1114.png'];
exportgraphics(gcf,dmdname,'resolution',600)

%% Find largest frequencies
% sort amplitudes in size order
% dvals = sort(amp(:), 'descend');
dvals = sort(ampr, 'descend');

% number of values to look at
nvals = 12;
largest = dvals(1:nvals);

% find index in amp/b2 where the largest values are
ind = zeros(nvals,1);
for k = 1:nvals
%     ind(k,1) = find(amp == largest(k));
    ind(k,1) = find(ampr == largest(k));
end

% attach frequency to amp index
ampf = [ind, dmdFreq(ind)];

% find index of positive frequencies
posampf_ind = imag(ampf(:,2)) >= 0;

% retain positive frequency indices
posf_ind = ampf(posampf_ind);

%% Inspect modes at highest amplitudes

for kk = 1:length(posf_ind)
    
    % DMD mode * DMD eigenvalue * amplitude
%     recon = P2(:,posf_ind(kk))*L2(posf_ind(kk),posf_ind(kk))*b2(posf_ind(kk),1);
    recon = P2(:,posf_ind(kk))*ampr(posf_ind(kk),1);
    
    % Separate U and V velocities (real)
    locs = size(recon);
    Ucol = real(recon(1:locs/2));
    Vcol = real(recon(locs/2+1:end));

    % Separate U and V velocities (imag)
%     locs = size(recon);
%     Ucol = imag(recon(1:locs/2));
%     Vcol = imag(recon(locs/2+1:end));
    
    % Reshape into PIV grid format (real)
    pivx = size(pivgrid,1);
    pivz = size(pivgrid,2);
    dmdplotu = NaN(pivx,pivz);
    dmdplotv = NaN(pivx,pivz);
    dmdplotu(originalIndex) = Ucol;
    dmdplotv(originalIndex) = Vcol;

    % Reshape into PIV grid format (imag)
%     dmdplotui = NaN(pivx,pivz);
%     dmdplotvi = NaN(pivx,pivz);
%     dmdplotui(originalIndex) = Ucoli;
%     dmdplotvi(originalIndex) = Vcoli;
    
    figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, dmdplotu, dmdplotv, figureprop );    
    
%     figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, dmdplotui, dmdplotvi, figureprop );    
    
%     t = {['DMD mode ', num2str(posf_ind(kk)), ', ', num2str(round(imag(dmdFreq(posf_ind(kk))))),' Hz']};
%     title(t)

end

%% Plot mode shapes (pairs)


%% WARNING
% Careful not to overwrite images!!

%% Inspect at 0 Hz freq
% near-zero frequencies
zeroind = imag(dmdFreq) <= 0.2 & imag(dmdFreq) >= 0;
zeroindf = find(zeroind == 1);

for jj = 1:length(zeroindf)
    
    % DMD mode * DMD eigenvalue * amplitude
%     recon = P2(:,posf_ind(kk))*L2(posf_ind(kk),posf_ind(kk))*b2(posf_ind(kk),1);
%     recon = P2(:,zeroindf(jj))*L2(zeroindf(jj),zeroindf(jj))*amp(zeroindf(jj),1);
%     recon = P2(:,zeroindf(jj))*-amp(zeroindf(jj),1);
    recon = P2(:,zeroindf(jj))*ampr(zeroindf(jj),1);

    
    % Separate U and V velocities (real)
    locs = size(recon);
    Ucol = real(recon(1:locs/2));
    Vcol = real(recon(locs/2+1:end));

    % Separate U and V velocities (imag)
    locs = size(recon);
    Ucoli = imag(recon(1:locs/2));
    Vcoli = imag(recon(locs/2+1:end));
    
    % Reshape into PIV grid format
    pivx = size(pivgrid,1);
    pivz = size(pivgrid,2);
    dmdplotu = NaN(pivx,pivz);
    dmdplotv = NaN(pivx,pivz);
    dmdplotu(originalIndex) = Ucol;
    dmdplotv(originalIndex) = Vcol;

    % Reshape into PIV grid format (imag)
    dmdplotui = NaN(pivx,pivz);
    dmdplotvi = NaN(pivx,pivz);
    dmdplotui(originalIndex) = Ucoli;
    dmdplotvi(originalIndex) = Vcoli;
    
    figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, dmdplotu, dmdplotv, figureprop );  
    set(gcf, 'position', [440 377 560 290])
%     t = {['DMD mode ', num2str(zeroindf(jj)), ', ', num2str(round(imag(dmdFreq(zeroindf(jj))))),' Hz']};
%     title(t)
    figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, dmdplotui, dmdplotvi, figureprop );    

    set(gcf, 'position', [440 377 560 290])

end

%%
dmdname = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/CTP/DMD',num2str(zeroindf(jj)),'.png'];
exportgraphics(gcf,dmdname,'resolution',600)

%%
Pa = P2(:,292);
Pb = P2(:,11);
Pc = P2(:,65);
Pd = P2(:,296);
Pe = P2(:,73);
P = [Pa Pb Pc Pd Pe];
L = zeros(5,5);
L(1,1) = L2(292,292);
L(2,2) = L2(11,11);
L(3,3) = L2(65,65);
L(4,4) = L2(296,296);
L(5,5) = L2(73,73);
b = zeros(4,1);
b(1,1) = b2(292,1);
b(2,1) = b2(11,1);
b(3,1) = b2(65,1);
b(4,1) = b2(296,1);
b(5,1) = b2(73,1);
data = P*L*b;

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
t = {['All modes reconstructed']};
title(t)


%% Trying different r values
r11 = 50;
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
set(gca, 'FontName','times', 'FontSize', 26,'linewidth',1.5)
scatter(real(diag(L2)),imag(diag(L2)),80,[0, 0.4470, 0.7410],'filled')
scatter(real(diag(L2(298,298))),imag(diag(L2(298,298))),120,[0.466, 0.6740, 0.1880],'filled','square')
scatter(real(diag(L2(292,292))),imag(diag(L2(292,292))),120,[0.494, 0.1840, 0.5560],'filled','square')
scatter(real(diag(L2(11,11))),imag(diag(L2(11,11))),120,[0.9290 0.6940 0.1250],'filled','square')
scatter(real(diag(L2(14,14))),imag(diag(L2(14,14))),120,[0.8500 0.3250 0.0980],'filled','square')
axis([-1.1 1.1 0 1.1])
set(gcf,'position',[207 377 872 420])
xticks([-1:0.2:1])
xlabel('Real')
ylabel('Imaginary')

%% export fig
name = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/unitcircle.png'];
exportgraphics(gcf,name,'resolution',600)

%% now forward/backward dmd

fblambdas = [];

% compute forward DMD
[U, Sf, Vf] = svd(Y1, 'econ');
f_Atilde = U' * Y2 * Vf / Sf;

% compute backward DMD
[U, S, V] = svd(Y2, 'econ');
b_Atilde = U' * Y1 * V / S;

% estimate forward/backward DMD
Atilde = (f_Atilde * inv(b_Atilde)) ^ 0.5;
[W, D] = eig(Atilde);
lam = diag(D);

fblambdas = [fblambdas; lam];

phi = pinv(fblambdas) * Y2 * V * inv(S) * W;

%%
% rpm
rpm = 1500;

% time between snapshots. *2 factor accounts for 4 stroke engine (seconds
% per cycle)
dt = 1/(rpm/60)*2;

% get DMD frequency
dmdFreq = log(fblambdas)/dt/2/pi;

% get amplitude
amp = abs(b);

% % scatter raw
% figure
% semilogy(imag(dmdFreq), amp,'o','MarkerFaceColor','blue')
% xlim([0 1])

% scatter neat
figure1 = figure( 'Color', [ 1 1 1 ] );
axes1 = axes( 'Parent', figure1, 'FontName', 'Times', 'FontSize', 20,'linewidth',1);
hold( axes1, 'on' )
box( axes1, 'on')
semilogy(imag(dmdFreq), amp,'filled','MarkerFaceColor','blue')
% xlim([0 6.5])
xlabel('Frequency (Hz)')
ylabel('DMD mode amplitude')
% title({'Raw spectrum:','cross-tumble plane, -285 CAD'})
title({'De-noised spectrum:','cross-tumble plane, -285 CAD'})