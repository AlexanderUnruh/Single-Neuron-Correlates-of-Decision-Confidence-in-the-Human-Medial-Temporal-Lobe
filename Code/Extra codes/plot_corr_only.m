function[psth]= plot_corr_only(mfr1, mfr2, val1, val2, val3, val4,  x_ax_lim, endline, endline2,  xlabel_in,variable, fig1, posSub, rangePlot, ylabelon)



dbstop if error
axfontsize = 12;

set(0, 'currentFigure',fig1)

%annotation('textbox',[0.06 0.958 0.2 0.05],'string',variable,'Units','normalized','FontName','Helvetica','FontSize',24,'EdgeColor','none')
% annotation('textbox',[0.1 0.958 0.1 0.05],'string',variable,'Units','normalized','FontName','Helvetica','FontSize',20,'EdgeColor','none')


set(gcf, 'PaperPositionMode', 'auto');





%% compute correlation coefficient



[c1, p1] = corr(mfr1,val1, 'type', 'Spearman');
[c2, p2] = corr(mfr2,val2, 'type', 'Spearman');


c11= abs((c1+c2)/2);




%subplot('Position',  [0.1300    0.1100    0.7750    0.7850]) % do not delete, old position plot
subplot('Position',  posSub)
hold on
set(gca,'FontSize',axfontsize,'FontWeight', 'bold');




if endline2 ==1
    
    area([0, 2000], [1,1], -1, 'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')

elseif endline2 == 2
    
    area([   0, 1000] , [1, 1], -1,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
    area([1200, 2200] , [1, 1], -1,  'Facealpha', 0.5, 'Facecolor', [0.3,0.3,0.3], 'Edgecolor', 'none')
    area([2400, 3400] , [1, 1], -1,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
    
end



rangePlot = rangePlot + endline;
plot(rangePlot,c11, 'color', [0, 0.4470, 0.7410], 'linewidth', 1.5)



axis tight
ylim([0, 0.6])
xlim(x_ax_lim + endline)

if ylabelon == 1
	ylabel('|Mean Rho|')
end





[c3, p3] = corr(mfr1,val3, 'type', 'Spearman');
[c4, p4] = corr(mfr2,val4, 'type', 'Spearman');


c21= abs((c3+c4)/2);

rangePlot = [-4999 : 10: 7400] + endline;
plot(rangePlot,c21, 'color', [0.9290 0.5940 0.1250], 'linewidth', 1.5)



%% A value

colormap((hot))
xlabel(xlabel_in)
set(gca, 'Xtick', [-7000:1000:7000], 'Xticklabel', {'-7','-6','-5','-4','-3', '-2', '-1', '0', '1', '2', '3', '4', '5', '6', '7'})


alpha  =  sqrt(0.05*2);
p2plot =  ~(p1<alpha & p2<alpha & (sign(c1) == sign(c2)));

%subplot('Position',  [0.1300    0.8950    0.7750    0.05])  % do not delete, old postition plot
ax11 = subplot('Position', [posSub(1), posSub(2)+posSub(4)+0.001 ,posSub(3), 0.05]);
imagesc(rangePlot, 1:5, p2plot', [0 1]);
colormap(ax11, [0, 0.4470, 0.7410; 1,1,1])
xlim(x_ax_lim)
axis off





%% B value 

xlabel(xlabel_in)
set(gca, 'Xtick', [-7000:1000:7000], 'Xticklabel', {'-7','-6','-5','-4','-3', '-2', '-1', '0', '1', '2', '3', '4', '5', '6', '7'})


alpha  =  sqrt(0.05*2);
p2plot =  ~(p3<alpha & p4<alpha & (sign(c3) == sign(c4)));

%subplot('Position',  [0.1300    0.8950    0.7750    0.05])  % do not delete, old postition plot
ax12 = subplot('Position', [posSub(1), posSub(2)+posSub(4)+0.051 ,posSub(3), 0.05]);
imagesc(rangePlot, 1:5, p2plot', [0 1]);
colormap(ax12, [0.9290 0.5940 0.1250; 1,1,1])
xlim(x_ax_lim)
axis off

asfdasdfasdf

