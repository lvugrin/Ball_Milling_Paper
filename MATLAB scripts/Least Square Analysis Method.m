% Data analysis
a = New_spectra_without_jar; %Monitoring_spectra_without_jar
sp = a(280:900, :);
wn = wnHS(280:900, 1);

% Removing Cosmic Rays
sp = Cosmic(sp, 205, 415);
sp = Cosmic(sp, 440, 440);

figure, plot(wn, sp);

% Baseline correction
[spals, als] = ALSbaseline(sp, 50000, 0.01, 30);
if (min(min(spals))) < 0
    NNALS = NonNegativity(spals);
end

figure, 
plot(wn, sp(:, 5)), hold on, 
plot(wn, NNALS(:, 5)), hold on,
plot(wn, als(:, 5)), hold off

% Normalisation
N2sp = Norm2(NNALS);

%%
% Time conversion
Times_minutes = Times ./ 60;

% Time definition over regions
interval1 = [1,2]; 
interval2 = [3,50]; 
interval3 = [50 , 90]; 
region1_indices = find(Times_minutes >= interval1(1) & Times_minutes <= interval1(2));
region2_indices = find(Times_minutes >= interval2(1) & Times_minutes <= interval2(2));
region3_indices = find(Times_minutes >= interval3(1) & Times_minutes <= interval3(2));

% Averaging spectra region
average_spectrum_region1 = mean(N2sp(:, region1_indices), 2);
average_spectrum_region2 = mean(N2sp(:, region2_indices), 2);
average_spectrum_region3 = mean(N2sp(:, region3_indices), 2);

% Matrix of components
component_spectra = [average_spectrum_region1, average_spectrum_region2, average_spectrum_region3]; % Dimenzije 621x3
X = component_spectra; % Dimensions X: 621x3

% Dimension checking
disp('Dimenzije X:');
disp(size(X)); % Estimated [621, 3]
disp('Dimenzije N2sp:');
disp(size(N2sp)); % Estimated [621, 338]

% Transpose
X_transpose = X'; % Dimenzije X_transpose: 3x621

% Dimension check
disp('Dimenzije X_transpose:');
disp(size(X_transpose)); % Očekuje se [3, 621]

% LSM formula
beta = (X_transpose * X) \ (X_transpose * N2sp);

% Beta factor
disp('Dimenzije beta:');
disp(size(beta)); % Estimated[3, 338]

%  Savitzky-Golay filter , smoothing section
%smooth_beta = sgolayfilt(beta', 3, 11)'; 

windowSize = 50; 
smooth_beta = movmean(beta, windowSize, 2); 

% Normalizaiton of beta factor
normalized_beta = normalize(smooth_beta, 2, 'range');

% Plotting
figure;
plot(Times_minutes, normalized_beta(1, :),   'DisplayName','Starting', 'LineWidth', 1.5);
hold on;
plot(Times_minutes,  normalized_beta(2, :), 'DisplayName', 'Melting', 'LineWidth', 1.5);
plot(Times_minutes, normalized_beta(3, :),   'DisplayName','Final', 'LineWidth', 1.5);
xlabel('Time / min');
ylabel('Coefficients');
legend ('boxoff');
title('Profile of Components');
hold off;

% Display
disp('Coefficients over time:');
disp(smooth_beta');
