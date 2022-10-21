%%
maxindices = zeros(300,1);
for kk = 1:300
    str = num2str(kk);
    str2 = ['cycle',str];
    maximum = max(ri2.(str2)(:,2));
    maxind = find(ri2.(str2)(:,2) == maximum);
    maxpodmode(kk,1) = maxind - 1;
end
nnz(maxpodmode)

%% Scatter plot
figure3 = figure( 'Color', [ 1 1 1 ] );
axes3 = axes( 'Parent', figure3, 'FontSize', 20,'linewidth',1);
hold( axes3, 'on' )
box( axes3, 'on')
scatter(1:300,maxpodmode(:,1),'filled','o')
% xlim([0 10])
xlabel('Cycle number [-]')
ylabel('POD mode [-]')
t = {['Cross-tumble CFD against PIV POD modes'], ['Optimum mode for each cycle']};
title(t,'FontSize',20)
fontname(figure3,"times")

%% Histogram
figure4 = figure( 'Color', [ 1 1 1 ] );
axes4 = axes( 'Parent', figure4, 'FontSize', 24, 'fontname', 'times','linewidth',1);
hold( axes4, 'on' )
box( axes4, 'on')
histogram(maxpodmode)
xlabel('Optimum POD mode')
ylabel('Number of cycles')
t = {['Frequency of optimum POD modes']};
% title(t,'FontSize',20)
% fontname(gcf,"times")
xticks([0:1:12])
xlim([-0.5 12.5])
set(gcf,'Position',[439 377 608 420])

%%
pivname = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/CTP/CTP_histogram.png'];
exportgraphics(gcf,pivname,'resolution',600)
