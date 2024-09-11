%% PXRD - importing raw data 

A = importdata('18c6-rofyiv.txt');

%% Plot from xy file

figure,
plot(forplottingwithphases(:,1),forplottingwithphases(:,2),  'linewidth', 1.0), hold on  %experimental
plot(forplottingwithphases(:,1),forplottingwithphases(:,3), 'linewidth',1.0), hold on   %calculated
plot(forplottingwithphases(:,1),forplottingwithphases(:,4),  'linewidth',0.00005), hold on %bckg
plot(forplottingwithphases(:,1),forplottingwithphases(:,5), 'linewidth',1.0), hold on %monohydrate complex
plot(forplottingwithphases(:,1),forplottingwithphases(:,7) 'linewidth',1.0), hold on %anhydrous
plot(forplottingwithphases(:,1),forplottingwithphases(:,6), 'linewidth',1.0), hold off %kcl

legend('Experimental', 'Calculated','Difference', '[18c6K]ClxH_2_O','[18c6K]Cl','KCl','LineWidth', 4);
legend('boxoff')
%xlim([5 60])
xlabel('2Theta / °');
ylabel('Intensity / counts')
ax = gca;
set(ax, 'Box','off', 'tickdir','out', 'xminortick','on','linewidth',1)
% axis tight