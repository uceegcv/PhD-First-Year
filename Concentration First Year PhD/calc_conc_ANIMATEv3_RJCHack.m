%% Calculate concentration from raw CYRIL data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters to check before running script%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optode_dist = 2.3433; %cm
DPF = 4.99; %4.99 baby head, 6.26 adult head, 4.16 adult arm, 5.51 adult leg (from Duncan 1994)
raw = xlsread('Intensity_LUMO_cyril.xlsx'); % Have manually copied raw data file of interest into Matlab and renamed to 'raw'
time = raw(:,1);
raw = raw(:,2:6);

%% Start calculation
% Read in raw CYRIL files
% to be written as necessary %

% Read in calibration file
%wl_cal = xlsread('Intensity_LUMO_cyril.xlsx');


wl_cal = [720 760 800 850 890];

wl = [720 760 800 850 890]; % Wavelengths used

% Extinction coefficients
% Wavelength, HbO2, HHb, CCO, wavelength dependency
%[epsilon_labview_770to906] = specific_extinction_coeffs_770to906;
%%

ext_coeffs = tissue_specific_extinction_coefficient_650to1000;
ext_coeffs_lambda = 650:1000;

path_dep = DPF_Lambda_Dependency_740to915;
dpf_lambda = 740:915;

for i = 1:length(wl)
    wavNow = wl(i);
    [~,ind] = min(abs(ext_coeffs_lambda-wavNow));
    ext_coeffs_crop(i,:) = ext_coeffs(ind,3:5);
    
    [~,ind] = min(abs(dpf_lambda-wavNow));
    dpfCorrFact(i,:) = path_dep(ind,2);
end

% ind_720 = 71;
% ext_720 = ext_coeffs(71,3:5);
% 
% ind_760 = 111;
% ext_760 = ext_coeffs(ind_760,3:5);
% dpf_760 = path_dep(21,2);
% 
% ind_800 = 151
% ext_800 = ext_coeffs(ind_800,3:5);
% dpf_800 = path_dep(61,2);
% 
% ind_850 = 201
% ext_850 = ext_coeffs(ind_850,3:5);
% dpf_850 = path_dep(111,2);
% 
% ind_890 = 241
% ext_890 = ext_coeffs(ind_890,3:5);
% dpf_890 = path_dep(151,2);
% ext_coeffs_crop = [ext_720;ext_760;ext_800;ext_850;ext_890];
% HbO2, HHb, CCO
% need to look at this
%dpfs = [1.05;dpf_760;dpf_800;dpf_850;dpf_890];


%% Invert extinctions
ext_coeffs_inv = pinv(ext_coeffs_crop);

%% Calculate concentration

[x,y] = size(raw);

% Attenuation
for ii = 1:x % row
    for jj = 1:y % column
        atten_WD(ii,jj) = log10(raw(1,jj)/raw(ii,jj))/(optode_dist*DPF*dpfCorrFact(jj)); % log 10
    end
end
%%
for ii = 1:x
    conc(:,ii) =  ext_coeffs_inv * atten_WD(ii,:)';
end
%%
conc = conc';

figure
FigH = figure('Position', get(0, 'Screensize'));
F    = getframe(FigH);
imwrite(F.cdata, 'concentrations.png', 'png')
plot(time,conc,'LineWidth',2)
xlabel('Time (s)')
ylabel('Change in Concentration (M)')
legend({'HbO','HbR','CCO'},'Location','Best')
%title('Average Counts at Detector')
%xlim([600 1000])
ax = gca;
ax.FontSize = 25;
ax.ColorOrder = [1 0 0; 0 0 1; 0 1 0];
saveas(gcf,'concentrations.png')

% Output  calculated concentration to excel files
% to be written %