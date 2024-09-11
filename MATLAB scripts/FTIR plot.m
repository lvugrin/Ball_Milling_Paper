%% import data

a = importdata('18c6ether.txt') ;   % load the first dataset 

%% plotting

figure,
ax1 = subplot(2,1,1);   
plot(ax1,a(:,1),a(:,2),'LineWidth',2.5), hold on        
box off,                                           
xlim([400 4000]);                                         
ylabel('\it T \rm/ a.u.', 'FontSize', 14);                
               
axgCp=gca;                                                   
axgCp.TickDir = 'out';  axgCp.XMinorTick='off'; axgCp.YMinorTick='off'; axgCp.TickLength=[0.010, 0.0250];      
axgCp.LineWidth=0.2;                                                
axgCp.FontSize=28;
axgCp.XDir='reverse';                                                 
dummyh = line('Linestyle', 'none', 'Marker', 'none', 'Color', 'none');                       
legend(dummyh,'a)','Location','northwestoutside');                                        
legend ('boxoff'), hold on

ax2 = subplot(2,1,2);   
plot(ax2,b(:,1),b(:,2),'LineWidth',2.5),hold on
box off 
xlim([400 4000]);
ylabel('\it T \rm/ a.u.', 'FontSize', 14);
xlabel('wave number / cm^-^1', 'FontSize', 14); 
axgCp=gca;
axgCp.TickDir = 'out';  axgCp.XMinorTick='off'; axgCp.YMinorTick='off'; axgCp.TickLength=[0.010, 0.0250];
axgCp.LineWidth=0.2; 
axgCp.FontSize=28;
axgCp.XDir='reverse';
dummyh = line('Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh,'b)','Location','northwestoutside');
legend ('boxoff'); hold on

ax3 = subplot(2,1,2);    
plot(ax3,c(:,1),c(:,2),'LineWidth',2.5),hold on
box off 
xlim([400 4000]);
ylabel('\it T \rm/ a.u.', 'FontSize', 14);
xlabel('wave number / cm^-^1', 'FontSize', 14); 
axgCp=gca;
axgCp.TickDir = 'out';  axgCp.XMinorTick='off'; axgCp.YMinorTick='off'; axgCp.TickLength=[0.010, 0.0250];
axgCp.LineWidth=0.2; 
axgCp.FontSize=28;
axgCp.XDir='reverse';
dummyh = line('Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh,'b)','Location','northwestoutside');
legend ('boxoff'); hold on


linkaxes([ax1,ax2],'x');
xticklabels(ax1,{});  

% end;