%% Importing data
a=Monitoring_spectra_without_jar;
sp=a(280:900,2:end);
wn=wnHS(280:900,1);

%% Removing Cosmic Rays
%sp(818,:) = (sp(817,:)+sp(819,:))./2;
%sp(1228,:) = (sp(1227,:)+sp(1229,:))./2;
% sp=Cosmic(sp,155, 205);
figure,plot(wn,sp);

%% baseline editing
[spals, als] = ALSbaseline(sp, 30000, 0.01, 30);
if (min(min(spals))) < 0, NNALS = NonNegativity(spals); end;
figure, plot(NNALS);

%% NORM
%NNALS = [NNALS1 ; NNALS2 ; NNALS3];
N2sp = Norm2(NNALS);

%% 2d plot
x1=100;
x2=size(N2sp,1); %change 
t1=1;
t2=size(N2sp,2); %full reaction time

%Vrijeme = Times(1:end,1)./60;

figure,
surf(Times(t1:t2,1)./60, wn(x1:x2,1), N2sp(x1:x2,t1:t2).^(1.1), 'FaceColor','interp','FaceLighting','gouraud','EdgeColor','none' );
colormap(jet(30000));
% colormap(hot);
shading interp;
view(90,-90);
xlabel('Time / min');
ylabel('Raman shift / cm^{-1}')
% title('LV-1VI')
set(gca, 'yDir','reverse')
axis tight


%% To improve intensity of Raman bands
x1 = 100;
% x2 = 600;
x2=size(N2sp,1);

t1=1;
%t2=size(N2sp,2);
t2=50; 

% Extract relevant data
time_range = Times(t1:t2, 1) ./ 60;
wn_range = wn(x1:x2, 1);
spectra_subset = N2sp(x1:x2, t1:t2);

% Logarithmic transformation to the spectra data
spectra_subset = log(1 + spectra_subset);

% 2d plot
figure;
surf(time_range, wn_range, spectra_subset, 'FaceColor', 'interp', 'EdgeColor', 'none');
colormap(parula); % testing different colormaps like 'copper' or 'gray', hot

% Enhancing contrast by manually adjusting the factor as needed
caxis([min(spectra_subset(:)), max(spectra_subset(:))] * 1.2); 

shading interp;
view(90, -90);
xlabel('Time / min');
ylabel('Raman shift / cm^{-1}');
set(gca, 'yDir', 'reverse');
xlim([min(time_range), max(time_range)]);
ylim([min(wn_range), max(wn_range)]);
