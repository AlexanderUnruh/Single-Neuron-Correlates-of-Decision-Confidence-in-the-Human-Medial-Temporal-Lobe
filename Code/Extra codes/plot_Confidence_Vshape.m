function [subplothandle] = plot_Confidence_Vshape(mfr1, mfr2,  val1, val2,  xlabel_in,  fig1, posSub,  ylabelon, possible_colors, legendtype)

                

dbstop if error
alpha= sqrt(2 * 0.05);
 

axfontsize = 12;

set(0, 'currentFigure',fig1)



set(gcf, 'PaperPositionMode', 'auto');




% possible_colors= [0   0  0.85;...
%     0.3 0.45 1;...
%     1 0.5,0.4;...
%     0.99   0   0;...
%     0.99   0   0;...
%     1 0.5,0.4;...
%       0.3 0.45 1;...
%     0   0  0.85];
% 
% 
% possible_colors = [[0 0.4470 0.7410];...
%     [0.8500 0.3250 0.0980];...
%     [0.9290 0.5940 0.1250];...
%     [0.4940 0.1840 0.5560];...
%     [0.4940 0.1840 0.5560];...
%     [0.9290 0.5940 0.1250];...
%     [0.8500 0.3250 0.0980];...
%     [0 0.4470 0.7410]];


possible_colors = [possible_colors; flipud(possible_colors)];

% compute firing rate by condition


mfr3 = [mfr1(1:end-1,:);mfr2(1:end-1,:)];
val3 = [val1(1:end-1,:);val2(1:end-1,:)];

mfr4 = [mfr1(2:end,:); mfr2(2:end,:)];
val4 = [val1(1:end-1,1); val2(1:end-1,1)];


% define windows of analysis

arange3 = [621:1:1201];
arange4 = [39:1:500];


% arange3 = [841:1:1201];
% arange4 = [140:1:500];



% only if significant extract the correlation
vecrange3 =  ([ zeros(621-1,1); ones([1202-621],1); zeros(1241-1202,1)]);
vecrange4 =  ([ zeros(39-1,1); ones([501-39],1); zeros(1241-501,1)]);
[c11, p11] = corr(mfr1, abs(val1), 'type', 'Spearman');
[c12, p12] = corr(mfr2, abs(val2), 'type', 'Spearman');

c1= abs((c11+c12)/2);
p1 =  (p11<alpha & p12<alpha & (sign(c11) == sign(c12)));


% next trial
[c21, p21] = corr(mfr1(2:end,:), abs(val1(1:end-1)), 'type', 'Spearman');
[c22, p22] = corr(mfr2(2:end,:), abs(val2(1:end-1)), 'type', 'Spearman');

c2= abs((c21+c22)/2);
p2 =  (p21<alpha & p22<alpha & (sign(c21) == sign(c22)));


arange3 = find(vecrange3 & p1); % current trial 
arange4 = find(vecrange4 & p2);  % next trial 






% extract mfr for plot 


mean3 = nanmean(mfr3(:, arange3),2);
mean4 = nanmean(mfr4(:, arange4),2);
meanAll = nanmean([mfr3(:, arange3) , mfr4(:, arange4)],2);


gg=1;
for b = [100:-25:25, -25:-25:-100]  % Print
    
    plotMFR3 (gg,:)=  nanmean(mfr3(val3 == b,:),1);
    plotMFR4 (gg,:)=  nanmean(mfr4(val4 == b,:),1);
    plotMFRse3 (gg,:)= nanstd(mfr3(val3 == b,:),1)/sqrt(sum(val3 == b));
    plotMFRse4 (gg,:)= nanstd(mfr4(val4 == b,:),1)/sqrt(sum(val4 == b));
    
    toplot(gg,:)  =  nanmean(meanAll(val4 == b,:),1);
    toplotSE(gg,:) = nanstd(meanAll(val4 == b,:),1)/sqrt(sum(val4 == b));
    
    gg=gg+1;
    
end


toplot1 = mean(plotMFR3(:, arange3),2) + mean(plotMFR4(:, arange4),2);
toplotSE1 = mean(plotMFRse3(:, arange3),2) + mean(plotMFRse4(:, arange4),2);

  
 
%  subplot('Position',  [0.1300    0.1100    0.7750    0.7850]) % do not delete, old position plot
subplothandle =subplot('Position',  posSub);


hold on
set(gca,'FontSize',axfontsize,'FontWeight', 'bold');



for gg=1:8
    eb = errorbar([gg],toplot(gg), toplotSE(gg), toplotSE(gg)  ,'.', 'linewidth', 1.5, 'color', possible_colors(gg,:));
    
    plot([gg], toplot(gg), '.', 'markersize',20, 'color', possible_colors(gg,:));
end

axis tight
xlim([0.5, 8.5])
line([4.5, 4.5], [0, subplothandle.YLim(2)], 'linewidth', 1.5, 'linestyle', '--', 'color', [0.6,0.6,0.6])

if ylabelon == 1
	ylabel('Firing Rate \pm SEM (spikes/s)')
end


xlabel(xlabel_in)
if legendtype == 1
    set(gca, 'Xtick', [1:8], 'Xticklabel', {'100','75','50','25',  '25', '50', '75', '100'})
elseif legendtype == 2
    set(gca, 'Xtick', [1:8], 'Xticklabel', {'S','mid-S','mid-L','L', 'L', 'mid-L', 'mid-S', 'S'})
end



