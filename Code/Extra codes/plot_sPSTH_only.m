function[psth, tickspos, leg]= plot_sPSTH_only(mfr1, mfr2, val1,  val2,  x_ax_lim, endline,  endline2, xlabel_in, variable, fig1, posSub, rangePlot, ylabelon, doRho, possible_colors, legendtype)


dbstop if error



axfontsize = 12;

set(0, 'currentFigure',fig1)

%annotation('textbox',[0.06 0.958 0.2 0.05],'string',variable,'Units','normalized','FontName','Helvetica','FontSize',24,'EdgeColor','none')
% annotation('textbox',[0.1 0.958 0.1 0.05],'string',variable,'Units','normalized','FontName','Helvetica','FontSize',20,'EdgeColor','none')



set(gcf, 'PaperPositionMode', 'auto');
% compute firing rate by condition


alphafill=0.3;


mfr3 = [mfr1;mfr2];
val3 = [val1;val2];


gg=1;

for b = [100:-25:25]  % Print
    
    plotMFR1 (gg,:)= nanmean(mfr1(val1 == b,:),1);
    plotMFR2 (gg,:)= nanmean(mfr2(val2 == b,:),1);
    plotMFR3 (gg,:)= nanmean(mfr3(val3 == b,:),1);
    plotMFRse3 (gg,:)= nanstd(mfr3(val3 == b,:),1)/sqrt(sum(val3 == b));
    gg=gg+1;
end


convseup   =  plotMFR3 + plotMFRse3;
convsedown =  plotMFR3 - plotMFRse3;  % SE
%       convseup   = plotMFR3 +(1.96*plotMFRse3);
%       convsedown = plotMFR3 -(1.96*plotMFRse3); %         % CI 95
sefill = [convseup, convsedown(:,end:-1:1)];



%% compute correlation coefficient





% subplot('Position',  [0.1300    0.1100    0.7750    0.7850]) % do not delete, old position plot
psth = subplot('Position',  posSub);


hold on
set(gca,'FontSize',axfontsize,'FontWeight', 'bold');


rangePlot = rangePlot + endline;


plot(rangePlot,plotMFR3(1,:), 'linewidth', 1.5, 'color', possible_colors(1,:))
plot(rangePlot,plotMFR3(2,:), 'linewidth', 1.5, 'color', possible_colors(2,:))
plot(rangePlot,plotMFR3(3,:), 'linewidth', 1.5, 'color', possible_colors(3,:))
plot(rangePlot,plotMFR3(4,:), 'linewidth', 1.5, 'color', possible_colors(4,:))




max2use = max(max(sefill));
max2use = max2use ;

if endline2 ==1
    
    area([0, 2000], [max2use,max2use], 0, 'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
    
elseif endline2 == 2
    
    area([   0, 1000] , [max2use, max2use], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
    area([1200, 2200] , [max2use, max2use], 0,  'Facealpha', 0.5, 'Facecolor', [0.3,0.3,0.3], 'Edgecolor', 'none')
    area([2400, 3400] , [max2use, max2use], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
    
end



fill([rangePlot, rangePlot(end:-1:1)], sefill(1,:) ,possible_colors(1,:), 'FaceAlpha', alphafill, 'EdgeColor', 'none'); hold on;
fill([rangePlot, rangePlot(end:-1:1)], sefill(2,:) ,possible_colors(2,:), 'FaceAlpha', alphafill, 'EdgeColor', 'none'); hold on;
fill([rangePlot, rangePlot(end:-1:1)], sefill(3,:) ,possible_colors(3,:), 'FaceAlpha', alphafill, 'EdgeColor', 'none'); hold on;
fill([rangePlot, rangePlot(end:-1:1)], sefill(4,:) ,possible_colors(4,:), 'FaceAlpha', alphafill, 'EdgeColor', 'none'); hold on;

plot(rangePlot,plotMFR3(1,:), 'linewidth', 2, 'color', possible_colors(1,:))
plot(rangePlot,plotMFR3(2,:), 'linewidth', 2, 'color', possible_colors(2,:))
plot(rangePlot,plotMFR3(3,:), 'linewidth', 2, 'color', possible_colors(3,:))
plot(rangePlot,plotMFR3(4,:), 'linewidth', 2, 'color', possible_colors(4,:))

axis tight


ylim([0,max2use])
xlim(x_ax_lim + endline)

if ylabelon == 1
    ylabel('Firing Rate \pm SEM (spikes/s)')
    if doRho == 1 || doRho == 2
        
        if legendtype == 1
            leg= legend('100', '75', '50', '25');
        elseif legendtype == 2
            leg= legend('short', 'mid-short', 'mid-long', 'long');
        end
        
        leg.Position = [ 0.1044    0.6269    0.0712    0.1794];
        leg.FontSize = 11;
        leg.FontWeight = 'bold';
    end
else
    if doRho ==0
        
        if legendtype == 1
            leg= legend('100', '75', '50', '25');
        elseif legendtype == 2
            leg= legend('short', 'mid-short', 'mid-long', 'long');
        end
        leg.Position = [0.8744    0.6269    0.0712    0.1794];
        leg.FontSize = 11;
        leg.FontWeight = 'bold';
    end
end


xlabel(xlabel_in)
set(gca, 'Xtick', [-7000:1000:7000], 'Xticklabel', {'-7','-6','-5','-4','-3', '-2', '-1', '0', '1', '2', '3', '4', '5', '6', '7'})



if doRho == 1
    
    [c1, p1] = corr(mfr1,val1, 'type', 'Spearman');
    [c2, p2] = corr(mfr2,val2, 'type', 'Spearman');
    c3= abs((c1+c2)/2);
    
    
    colormap((hot))
    
    
    
    
    alpha  =  sqrt(0.05*2);
    p2plot =  ~(p1<alpha & p2<alpha & (sign(c1) == sign(c2)));
    
    %subplot('Position',  [0.1300    0.8950    0.7750    0.05])  % do not delete, old postition plot
    subplot('Position', [posSub(1), posSub(2)+posSub(4)+0.02 ,posSub(3), 0.05])
    imagesc(rangePlot, 1:5, p2plot', [0 1]);
    xlim(x_ax_lim)
    axis off
    
    
    % %     subplot('Position', [posSub(1), posSub(2)+posSub(4)+0.001+0.055 ,posSub(3), 0.05])  % this only if mean rho plotted above
    subplot('Position', [posSub(1), posSub(2)-0.18, posSub(3), 0.05])
    imagesc(rangePlot, 1:5, c3', [0 0.6]);
    
    
    xlim(x_ax_lim)
    axis off
    shiftdown = 0.18;
    
    if ylabelon ==0
        
        
        heightPos3 = posSub(2)-shiftdown ;
        
        poscbh = [0.953,   heightPos3,    0.0080,  0.12 -  heightPos3];
        cbh = colorbar;
        cbh.FontSize= 11;
        cbh.Position= poscbh;
        cbh.FontWeight= 'bold';
        cbh.Ticks=[0, 0.3, 0.6];
    else
        
        asdf=  annotation('textbox',[0.000,   posSub(2)-shiftdown+0.004 ,  0.2000  ,  0.0500],'string','|<Rho>|','Units','normalized','FontName','Helvetica','FontSize',11, 'FontWeight','bold', 'EdgeColor','none');
        pvalanno=  annotation('textbox',[0.00500,   posSub(2)+posSub(4)+0.02+0.004  ,  0.2000  ,  0.0500],'string','Signif.','Units','normalized','FontName','Helvetica','FontSize',11, 'FontWeight','bold', 'EdgeColor','none');
    end
end






if doRho == 2
    
    [c1, p1] = corr(mfr1,val1, 'type', 'Spearman');
    [c2, p2] = corr(mfr2,val2, 'type', 'Spearman');
    c3= abs((c1+c2)/2);
    
    
    colormap((hot))
    
    
    
    
    alpha  =  sqrt(0.05*2);
    p2plot =  ~(p1<alpha & p2<alpha & (sign(c1) == sign(c2)));
    
    %subplot('Position',  [0.1300    0.8950    0.7750    0.05])  % do not delete, old postition plot
    tickspos = subplot('Position', [posSub(1), posSub(2)+posSub(4)+0.001 ,posSub(3), 0.05]);
    imagesc(rangePlot, 1:5, p2plot', [0 1]);
    xlim(x_ax_lim)
    axis off
 
end

















