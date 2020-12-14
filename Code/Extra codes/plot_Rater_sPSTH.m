function[psth]= plot_Rater_sPSTH(mfr1, mfr2, val1,  val2,  x_ax_lim, endline,  endline2, xlabel_in, variable, fig1, posSub, rangePlot, ylabelon, doRho, possible_colors)


dbstop if error


asdfasdf
axfontsize = 12;

set(0, 'currentFigure',fig1)

%annotation('textbox',[0.06 0.958 0.2 0.05],'string',variable,'Units','normalized','FontName','Helvetica','FontSize',24,'EdgeColor','none')
% annotation('textbox',[0.1 0.958 0.1 0.05],'string',variable,'Units','normalized','FontName','Helvetica','FontSize',20,'EdgeColor','none')



set(gcf, 'PaperPositionMode', 'auto');


% extract rasters rate by condition



alphafill=0.3;


mfr3 = [mfr1;mfr2];
val3 = [val1;val2];

clear plotMFR*
gg=1;
for b = [100:-25:25]  % Print
    
    plotMFR1 {gg,:}= [mfr1(val1 == b,:)];
    plotMFR2 {gg,:}= [mfr2(val2 == b,:)];
    plotMFR3 {gg,:}= [mfr3(val3 == b,:)];
    
    gg=gg+1;
end





marksize=8;
% subplot('Position',  [0.1300    0.1100    0.7750    0.7850]) % do not delete, old position plot
psth = subplot('Position',  posSub);


hold on
set(gca,'FontSize',axfontsize,'FontWeight', 'bold');




% max2use = max(max(sefill));
max2use = 381 ;

if endline2 ==1
    
    area([0, 2000], [max2use,max2use], 0, 'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
    
elseif endline2 == 2
    
    area([   0, 1000] , [max2use, max2use], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none', 'ShowBaseLine', 'off')
    area([1200, 2200] , [max2use, max2use], 0,  'Facealpha', 0.5, 'Facecolor', [0.3,0.3,0.3], 'Edgecolor', 'none', 'ShowBaseLine', 'off')
    area([2400, 3400] , [max2use, max2use], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none', 'ShowBaseLine', 'off')
    
end


rangePlot = rangePlot + endline;
gg=1;
for trialn = 1:numel(plotMFR3{4})
    try
        plot(plotMFR3{4}{trialn}, gg,'.', 'markersize', marksize,'color', possible_colors(4,:))
    catch
    end
    gg=gg+1;
end

for trialn = 1:numel(plotMFR3{3})
    try
        plot(plotMFR3{3}{trialn}, gg, '.', 'markersize', marksize, 'color', possible_colors(3,:))
    catch
    end
    gg=gg+1;
end


for trialn = 1:numel(plotMFR3{2})
    try
        plot(plotMFR3{2}{trialn}, gg, '.', 'markersize', marksize, 'color', possible_colors(2,:))
    catch
    end
    gg=gg+1;
end



for trialn = 1:numel(plotMFR3{1})
    try
        plot(plotMFR3{1}{trialn}, gg, '.', 'markersize', marksize, 'color', possible_colors(1,:))
    catch
    end
    gg=gg+1;
end


axis tight


ylim([0,max2use])
xlim(x_ax_lim + endline)

if ylabelon == 1
    
    if doRho == 1 || doRho == 2
        ylabel('Trial Number')
        leg= legend('100', '75', '50', '25');
        leg.Position = [ 0.1044    0.6269    0.0712    0.1794];
        leg.FontSize = 11;
        leg.FontWeight = 'bold';
    end
else
    if doRho ==0
        leg= legend('100', '75', '50', '25');
        leg.Position = [0.8744    0.6269    0.0712    0.1794];
        leg.FontSize = 11;
        leg.FontWeight = 'bold';
    end
end


xlabel(xlabel_in)
set(gca, 'Xtick', [-7000:1000:7000], 'Xticklabel', {'-7','-6','-5','-4','-3', '-2', '-1', '0', '1', '2', '3', '4', '5', '6', '7'})



if d\oRho == 1
    
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
    subplot('Position', [posSub(1), posSub(2)+posSub(4)+0.001 ,posSub(3), 0.05])
    imagesc(rangePlot, 1:5, p2plot', [0 1]);
    xlim(x_ax_lim)
    axis off
    
end

















