%% Load
% load('../../PIV_matlab/example_PlotsAndMetrics/x20190517_CrossTumble_CR12p5_T1_C33_DVA_Motored.mat')
% load('../dataformatting/x20190517_CrossTumble_CR12p5_T1_C33_DVA_Motored.mat')
load('../dataformatting/x20180706_Tumble_CR12p5_T1_C33_DVA_Motored_Processed_all_masked.mat')

%% Get original data
CrankAngle_plot = -285;
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

% Zero velocity data
SPODData.U = zeros( nValidPoints, nCycle );
SPODData.V = zeros( nValidPoints, nCycle );

% For each cycle, fill in the SPOData velocities
for cycle = 1:nCycle
    temp_u = u(:,:,cycle);
    temp_w = w(:,:,cycle);
    SPODData.U( :,cycle ) = temp_u(originalIndex);
    SPODData.V( :,cycle ) = temp_w(originalIndex);
end

% Tomorrow: What to do with SPODData? Turn into complex number?

%% SPOD
% Turn into complex number and transpose, so time in first dimension for
% SPOD algorithm
complex_velo = [complex( SPODData.U, SPODData.V )];
complex_veloT = complex_velo';

% SPOD meaning:
% The columns of L contain the modal energy spectra. 
% P contains the SPOD modes whose spatial dimensions are identical to those 
% of X. The first index of P is the frequency and the last one the mode 
% number ranked in descending order by modal energy. 
% F is the frequency vector. 

% rpm
rpm = 1500;

% time between snapshots. *2 factor accounts for 4 stroke engine
dt = 1/(rpm/60)*2;

% SPOD
[L2,P2,f2] = spod(complex_veloT,[],[],[],dt);

% Why do I get negative frquencies? 
% But first, just try to plot the SPOD modes by calculating SPODApprox
% go through the algorithm step by step and plot

%% Calc SPOD Approx
originalrows = size(u,1);
originalcols = size(u,2);

PODcumsum = PODResult.Mode( :, 1:nModes ) * PODResult.Coeff( 1:nModes, CycleNo ) + repmat( PODResult.EnsembleMean, 1, length( CycleNo ) );

SPODApprox.U = nan( originalrows, originalcols, 300 );
SPODApprox.V = nan( originalrows, originalcols, 300 );
for jj = 1 : length( CycleNo )
    SPODApprox.U(:,:,1,jj) = Convert_PODFormat( real( PODcumsum(:,jj) ), 'POD2Original', PODResult.nRowsInOriginalGrid, PODResult.nColsInOriginalGrid, PODResult.IndexInOriginalGrid );
    SPODApprox.V(:,:,1,jj) = Convert_PODFormat( imag( PODcumsum(:,jj) ), 'POD2Original', PODResult.nRowsInOriginalGrid, PODResult.nColsInOriginalGrid, PODResult.IndexInOriginalGrid );
end

%%
% Second, we visualize the 1st and 2nd SPOD modes at three frequencies.
figure
count = 1;
for fi = [10 15 20]
    for mi = [1 2]
        subplot(3,2,count)
        contourf(x,r,real(squeeze(P(fi,:,:,mi))),11,'edgecolor','none'), axis equal tight, caxis(max(abs(caxis))*[-1 1])
        xlabel('x'), ylabel('r'), title(['$f=' num2str(f(fi),'%.2f$') ', mode ' num2str(mi) ', $\lambda=' num2str(L(fi,mi),'%.2g$')])
        xlim([0 10]); ylim([0 2])
        count = count + 1;
    end
end

%% For each cycle, turn velocities into column vectors
% Desired output
A = NaN(dim,300);
B = NaN(dim,300);
C = NaN(dim*2,300);
E = NaN(dim,300);

% For each cycle, turn velocities into column vectors
for cycle = 1:300
    A(1:dim,cycle) = reshape(v(:,:,cycle),[],1);
    a = A(1:dim,cycle);
    B(1:dim,cycle) = reshape(w(:,:,cycle),[],1);
    b = B(1:dim,cycle);
    E(:,cycle) = complex(a,b);
end

C(1:dim,:) = A;
C(dim+1:dim*2,:) = B;

% Need first dimension as time for SPOD
D = C';
F = E';

%% Reshape E back into grid
G = NaN(83,51,300);
for cycle = 1:300
    G(:,:,cycle) = reshape(E(:,cycle),83,51);
end

H = permute(G,[3 1 2]);

%% Get time in the first dimension for SPOD
vt = permute(v,[3 1 2]);
wt = permute(w,[3 1 2]);

%% Animate the first 100 snapshots of the pressure field.
figure('name','Pressure of the symmetric component of a turbulent jet')
for ti=1:10
    quiver(MaskedData.X,MaskedData.Y,real(complex_velo(ti,:)),imag(complex_velo(ti,:)))
%     axis equal tight, shading interp, caxis([4.43 4.48])
    xlabel('x'), ylabel('r')
    pause(0.05)
    drawnow
end

%% SPOD
