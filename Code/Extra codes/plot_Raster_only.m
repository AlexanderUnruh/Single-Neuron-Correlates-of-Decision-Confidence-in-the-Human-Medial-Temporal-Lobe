function[psth]= plot_Raster_only(mfr1, mfr2, val1,  val2,  x_ax_lim, endline,  endline2, xlabel_in, variable, fig1, posSub, rangePlot, ylabelon, doRho, possible_colors, marksize)


dbstop if error


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






% subplot('Position',  [0.1300    0.1100    0.7750    0.7850]) % do not delete, old position plot
psth = subplot('Position',  posSub);


hold on
set(gca,'FontSize',axfontsize,'FontWeight', 'bold');




% max2use = max(max(sefill));
max2use = size(mfr3,1) +4;

if endline2 ==1
    
    area([0, 2000], [max2use,max2use], -4, 'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none', 'ShowBaseLine', 'off')
    
elseif endline2 == 2
    
    area([   0, 1000] , [max2use, max2use], -4,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none', 'ShowBaseLine', 'off')
    area([1200, 2200] , [max2use, max2use], -4,  'Facealpha', 0.5, 'Facecolor', [0.3,0.3,0.3], 'Edgecolor', 'none', 'ShowBaseLine', 'off')
    area([2400, 3400] , [max2use, max2use], -4,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none', 'ShowBaseLine', 'off')
    
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


axis off

















