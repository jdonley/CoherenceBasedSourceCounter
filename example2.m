clear;
%% Run example.m first to generate the three .mat files listed below
%% After the .mat files have been gernated run example2.m

S{1}='Nsources_vs_SynthParams.ReverbTime_MSC.mat';
S{2}='Nsources_vs_SynthParams.SNR_MSC.mat';
S{3}='Nsources_vs_IntraNodeDist_MSC.mat';

%% Plot Success Rate (SR) results

markLst = {'+','o','.','^','x','s','d','*','v','>','<','p','h'};
lineLst = {': ','--','- ','-.'};
Cols = [0.1,0.1,0.8;
        0.8,0,0;];
methodtxt = {'MSC'};

Nmethods = numel(methodtxt);
fH = figure(2);
fH.Name = ['Success Rate - Comparison'];
fH.Color = 'w';
scl = 1;
ha = tightPlots(...
    1, 3, ...
    (88.9/10 + 6.35/10 + 88.9/10)*scl, ...
    [1 1], ...
    [0.5 0.2]*scl, ...
    [1.5 0.1]*scl, ...
    [1.5 0.1]*scl, ...
    'centimeters');

for a = 1:3
    load(S{a});
    Nsrc = numel(vals{1});
    
    axes(ha(a));
    ax = gca;
    ax.ColorOrder = [0, 0, 0];
    set(0,'DefaultAxesColorOrder', [0,0,0]);
    
    ax.LineStyleOrder = ...
        cell2mat([reshape(markLst(repmat(1:Nsrc,Nmethods,1))',[],1), ...
        reshape(lineLst(repmat(1:Nmethods,Nsrc,1)),[],1)]);
    set(0,'DefaultAxesLineStyleOrder', cell2mat([reshape(markLst(repmat(1:Nsrc,Nmethods,1))',[],1), ...
        reshape(lineLst(repmat(1:Nmethods,Nsrc,1)),[],1)]));
    
    pltNVpairs = {...
        'MarkerSize', 4,...
        'LineWidth', 0.5};
    
    plot(ax, vals{2},SR  ,'color',Cols(1,:), pltNVpairs{:}); hold on
    
    
    
    limbuff = 0.05; %10 percent
    xrange = [min(vals{2}) max(vals{2})];
    yrange = [0 100];
    % xlabel(variables{2}); xlim(xrange+[-1 1]*(limbuff*diff(xrange)));
    tx=xlabel(varnames{2}, 'interpreter','latex'); xlim(xrange+[-1 1]*(limbuff*diff(xrange)));
    ty=ylabel('Success Rate (SR) (\%)', 'interpreter','latex'); ylim(yrange+[-1 1]*(limbuff*diff(yrange)));
    
    titletxt = 'Comparisons';
    % title(titletxt);
    hold off
    
    if a ~= 1
        ax.YLabel=[];
        ax.YTickLabel=[];
    end
    
    ax.TickDir = 'both';
    ax.FontSize = 9;
    ax.FontName = 'Cambria';
    tx.FontName = 'Times';
    ty.FontName = 'Times';
    
end


hold on;
h = zeros(Nsrc+Nmethods, 1);
h(1:Nsrc) = plot(ax, NaN(Nsrc) ); set(h(1:Nsrc),'LineStyle', 'none');
ax.ColorOrder = [Cols(1,:); Cols(2,:)];
h(Nsrc+1:end) = plot(ax, NaN(Nmethods) ); set(h(Nsrc+1:end),'Marker', 'none');
arrayfun(@(i) set(h(Nsrc+i),'LineStyle',lineLst{i}),1:Nmethods);
srcTxt = [' ' varnames{1}];
hL = legend(h, ...
    [mat2cell([num2str([vals{1}]'), repmat(srcTxt,Nsrc,1)],ones(Nsrc,1),numel(srcTxt)+1); ...
    methodtxt'], ...
    'Location', 'southeast');

hL.FontSize = 9;
hL.FontName = 'Times';

hold off;

ax.Units = 'centimeters';
hL.Units = ax.Units;
hL.Position(1:2) = ax.Position(1:2) + ax.Position(3:4).*[1 0.5] - [-0.1, hL.Position(4)/2];

fH.Units = ax.Units;
fH.Position(3) = hL.Position(1) + hL.Position(3);
% fH.Units = 'centimeters';
% fH.Position(3) = 10;
% fH.Position(4) = fH.Position(3)*3/4;

%% Save results and save plot
print(['results' '.png'],'-dpng','-r600');
