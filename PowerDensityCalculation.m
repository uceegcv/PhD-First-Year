%% Script ot cacluate the power denisty of a broadband source
% need a power meter and a spectrometer sharing the same bandwidth so it can
% work.


%% 720 Data Import
% Find file names in directory
dir_name = cd;
%dir_name = '\gina_data\tile_1';
files = dir('*720*.txt');
%files = dir("gina_data\tile_5\LED_890nm_tile5_HR2B10321__6__05-57-11-243.txt");

% write to file
splitName = split(dir_name, "\");
name = splitName{end};

fileName = strcat(name,'.txt');

% Import all files in folder
for ii = 1:length(files)
    %data{ii,1} = importdata([dir_name,'\', files(ii,1).name]);
    data{ii,1} = importdata([files(ii,1).name]);
end

% extract lambda vector
Lambda=data{1,1}.data(:,1);

% background spectra
% Select intensity data for all background files
for ii = 1:10 % number of background files
    data_bg{ii,1} = data{ii+10,1}.data;
end

% Take average of background files intensity data
for jj = 1:length(data_bg{ii,1}(:,1))
    for ii = 1:length(data_bg)
        data_bg_mean{jj,:}(ii,1) = (data_bg{ii,1}(jj,2));
    end
end

for ii = 1:length(data_bg_mean)
    data_bg_mean2(ii,1) = mean(data_bg_mean{ii,1}(:,1));
end

% Light source spectra
% Select intensity data for all light source files
for ii = 1:10 % number of background files
    data_HL{ii,1} = data{ii,1}.data;
end

% figure
% figure(1)
% plot(Lambda,data_HL{1,1}(:,2))
% xlabel('Wavelength (nm)')
% ylabel('Detected counts')

%
% [pks,locs,w,p] = findpeaks(data_HL{1,1}(:,2),Lambda,'MinPeakHeight',500,'MinPeakWidth',15,'Annotate','extents');
% findpeaks(data_HL{1,1}(:,2),Lambda,'MinPeakHeight',500,'MinPeakWidth',15,'Annotate','extents')
% title('Signal Peak Widths')

% Take average of light files intensity data
for jj = 1:length(data_HL{ii,1}(:,1))
    for ii = 1:length(data_HL)
        data_HL_mean{jj,:}(ii,1) = (data_HL{ii,1}(jj,2));
    end
end

for ii = 1:length(data_HL_mean)
    data_HL_mean2(ii,1) = mean(data_HL_mean{ii,1}(:,1));
end


% preprocessing
% subtract mean spectrum and mean background
mean_spec_720 = data_HL_mean2 - data_bg_mean2;

% sum all intensities
mean_spec_int_720 = sum(mean_spec_720);

% plot mean spectra
% figure(2)
% plot(Lambda,mean_spec_720)
% xlabel('Wavelength (nm)')
% ylabel('Counts')
% title('Average Counts at Detector')
% xlim([600 950])
% ax = gca;
% ax.FontSize = 25;
%% 760
files_760 = dir('*760*.txt');
%files = dir("gina_data\tile_5\LED_890nm_tile5_HR2B10321__6__05-57-11-243.txt");

% write to file
splitName = split(dir_name, "\");
name = splitName{end};

fileName = strcat(name,'.txt');

% Import all files in folder
for ii = 1:length(files)
    %data{ii,1} = importdata([dir_name,'\', files(ii,1).name]);
    data_760{ii,1} = importdata([files_760(ii,1).name]);
end

% extract lambda vector
Lambda=data_760{1,1}.data(:,1);

% background spectra
% Select intensity data for all background files
for ii = 1:10 % number of background files
    data_bg_760{ii,1} = data_760{ii+10,1}.data;
end

% Take average of background files intensity data
for jj = 1:length(data_bg_760{ii,1}(:,1))
    for ii = 1:length(data_bg_760)
        data_bg_760_mean{jj,:}(ii,1) = (data_bg_760{ii,1}(jj,2));
    end
end

for ii = 1:length(data_bg_760_mean)
    data_bg_760_mean2(ii,1) = mean(data_bg_760_mean{ii,1}(:,1));
end

% Light source spectra
% Select intensity data for all light source files
for ii = 1:10 % number of background files
    data_HL_760{ii,1} = data_760{ii,1}.data;
end

% Take average of light files intensity data
for jj = 1:length(data_HL_760{ii,1}(:,1))
    for ii = 1:length(data_HL_760)
        data_HL_760_mean{jj,:}(ii,1) = (data_HL_760{ii,1}(jj,2));
    end
end

for ii = 1:length(data_HL_760_mean)
    data_HL_760_mean2(ii,1) = mean(data_HL_760_mean{ii,1}(:,1));
end

% preprocessing
% subtract mean spectrum and mean background
mean_spec_760 = data_HL_760_mean2 - data_bg_760_mean2;

% sum all intensities
mean_spec_int_760 = sum(mean_spec_760);

% plot mean spectra
% figure(3)
% plot(Lambda,mean_spec_760)
% xlabel('Wavelength (nm)')
% ylabel('Counts')
% title('Average Counts at Detector')
% xlim([600 950])
% ax = gca;
% ax.FontSize = 25;

%% 800

files_800 = dir('*800*.txt');
%files = dir("gina_data\tile_5\LED_890nm_tile5_HR2B10321__6__05-57-11-243.txt");

% write to file
splitName = split(dir_name, "\");
name = splitName{end};

fileName = strcat(name,'.txt');

% Import all files in folder
for ii = 1:length(files)
    %data{ii,1} = importdata([dir_name,'\', files(ii,1).name]);
    data_800{ii,1} = importdata([files_800(ii,1).name]);
end

% extract lambda vector
Lambda=data_800{1,1}.data(:,1);

% background spectra
% Select intensity data for all background files
for ii = 1:10 % number of background files
    data_bg_800{ii,1} = data_800{ii+10,1}.data;
end

% Take average of background files intensity data
for jj = 1:length(data_bg_800{ii,1}(:,1))
    for ii = 1:length(data_bg_800)
        data_bg_800_mean{jj,:}(ii,1) = (data_bg_800{ii,1}(jj,2));
    end
end

for ii = 1:length(data_bg_800_mean)
    data_bg_800_mean2(ii,1) = mean(data_bg_800_mean{ii,1}(:,1));
end

% Light source spectra
% Select intensity data for all light source files
for ii = 1:10 % number of background files
    data_HL_800{ii,1} = data_800{ii,1}.data;
end

% Take average of light files intensity data
for jj = 1:length(data_HL_800{ii,1}(:,1))
    for ii = 1:length(data_HL_800)
        data_HL_800_mean{jj,:}(ii,1) = (data_HL_800{ii,1}(jj,2));
    end
end

for ii = 1:length(data_HL_800_mean)
    data_HL_800_mean2(ii,1) = mean(data_HL_800_mean{ii,1}(:,1));
end

% preprocessing
% subtract mean spectrum and mean background
mean_spec_800 = data_HL_800_mean2 - data_bg_800_mean2;

% sum all intensities
mean_spec_int_800 = sum(mean_spec_800);

% plot mean spectra
% figure(3)
% plot(Lambda,mean_spec_800)
% xlabel('Wavelength (nm)')
% ylabel('Counts')
% title('Average Counts at Detector')
% xlim([600 950])
% ax = gca;
% ax.FontSize = 25;

%% 850

files_850 = dir('*850*.txt');
%files = dir("gina_data\tile_5\LED_890nm_tile5_HR2B10321__6__05-57-11-243.txt");

% write to file
splitName = split(dir_name, "\");
name = splitName{end};

fileName = strcat(name,'.txt');

% Import all files in folder
for ii = 1:length(files)
    %data{ii,1} = importdata([dir_name,'\', files(ii,1).name]);
    data_850{ii,1} = importdata([files_850(ii,1).name]);
end

% extract lambda vector
Lambda=data_850{1,1}.data(:,1);

% background spectra
% Select intensity data for all background files
for ii = 1:10 % number of background files
    data_bg_850{ii,1} = data_850{ii+10,1}.data;
end

% Take average of background files intensity data
for jj = 1:length(data_bg_850{ii,1}(:,1))
    for ii = 1:length(data_bg_850)
        data_bg_850_mean{jj,:}(ii,1) = (data_bg_850{ii,1}(jj,2));
    end
end

for ii = 1:length(data_bg_850_mean)
    data_bg_850_mean2(ii,1) = mean(data_bg_850_mean{ii,1}(:,1));
end

% Light source spectra
% Select intensity data for all light source files
for ii = 1:10 % number of background files
    data_HL_850{ii,1} = data_850{ii,1}.data;
end

% Take average of light files intensity data
for jj = 1:length(data_HL_850{ii,1}(:,1))
    for ii = 1:length(data_HL_850)
        data_HL_850_mean{jj,:}(ii,1) = (data_HL_850{ii,1}(jj,2));
    end
end

for ii = 1:length(data_HL_850_mean)
    data_HL_850_mean2(ii,1) = mean(data_HL_850_mean{ii,1}(:,1));
end

% preprocessing
% subtract mean spectrum and mean background
mean_spec_850 = data_HL_850_mean2 - data_bg_850_mean2;

% sum all intensities
mean_spec_int_850 = sum(mean_spec_850);

% plot mean spectra
% figure(3)
% plot(Lambda,mean_spec_850)
% xlabel('Wavelength (nm)')
% ylabel('Counts')
% title('Average Counts at Detector')
% xlim([600 950])
% ax = gca;
% ax.FontSize = 25;


%% 890

files_890 = dir('*890*.txt');
%files = dir("gina_data\tile_5\LED_890nm_tile5_HR2B10321__6__05-57-11-243.txt");

% write to file
splitName = split(dir_name, "\");
name = splitName{end};

fileName = strcat(name,'.txt');

% Import all files in folder
for ii = 1:length(files)
    %data{ii,1} = importdata([dir_name,'\', files(ii,1).name]);
    data_890{ii,1} = importdata([files_890(ii,1).name]);
end

% extract lambda vector
Lambda=data_890{1,1}.data(:,1);

% background spectra
% Select intensity data for all background files
for ii = 1:10 % number of background files
    data_bg_890{ii,1} = data_890{ii+10,1}.data;
end

% Take average of background files intensity data
for jj = 1:length(data_bg_890{ii,1}(:,1))
    for ii = 1:length(data_bg_890)
        data_bg_890_mean{jj,:}(ii,1) = (data_bg_890{ii,1}(jj,2));
    end
end

for ii = 1:length(data_bg_890_mean)
    data_bg_890_mean2(ii,1) = mean(data_bg_890_mean{ii,1}(:,1));
end

% Light source spectra
% Select intensity data for all light source files
for ii = 1:10 % number of background files
    data_HL_890{ii,1} = data_890{ii,1}.data;
end

% Take average of light files intensity data
for jj = 1:length(data_HL_890{ii,1}(:,1))
    for ii = 1:length(data_HL_890)
        data_HL_890_mean{jj,:}(ii,1) = (data_HL_890{ii,1}(jj,2));
    end
end

for ii = 1:length(data_HL_890_mean)
    data_HL_890_mean2(ii,1) = mean(data_HL_890_mean{ii,1}(:,1));
end

% preprocessing
% subtract mean spectrum and mean background
mean_spec_890 = data_HL_890_mean2 - data_bg_890_mean2;

% sum all intensities
mean_spec_int_890 = sum(mean_spec_890);

% plot mean spectra
% figure(3)
% plot(Lambda,mean_spec_890)
% xlabel('Wavelength (nm)')
% ylabel('Counts')
% title('Average Counts at Detector')
% xlim([600 950])
% ax = gca;
% ax.FontSize = 25;

%%

all_LED_plot_lambda = [Lambda mean_spec_720 mean_spec_760 mean_spec_800 mean_spec_850 mean_spec_890];
%%

all_LED_plot = [mean_spec_720 mean_spec_760 mean_spec_800 mean_spec_850 mean_spec_890];
all_LED_plot_vector = [mean_spec_720;mean_spec_760;mean_spec_800;mean_spec_850;mean_spec_890];

figure(7)
FigH = figure('Position', get(0, 'Screensize'));
F    = getframe(FigH);
imwrite(F.cdata, 'LED_spectra_plot2.png', 'png')
plot(Lambda,all_LED_plot,'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Counts')
title('Average Counts at Detector')
xlim([600 1000])
ax = gca;
ax.FontSize = 25;
%legend({'720nm - FWHM = 23nm','760nm - FWHM = 26nm','800nm - FWHM = 29nm','850nm - FWHM = 32nm','890nm - FWHM = 47nm'})
[pks1,locs1,w1,p1] = findpeaks(mean_spec_720, 'MinPeakHeight',500,'MinPeakWidth',15,'WidthReference','halfheight');      % Identify Prominent Peaks ...
[pks2,locs2] = findpeaks(mean_spec_760, 'MinPeakHeight',500,'MinPeakWidth',15);
[pks3,locs3] = findpeaks(mean_spec_800, 'MinPeakHeight',500,'MinPeakWidth',15);
[pks4,locs4] = findpeaks(mean_spec_850, 'MinPeakHeight',500,'MinPeakWidth',15);
[pks5,locs5] = findpeaks(mean_spec_890, 'MinPeakHeight',500,'MinPeakWidth',15);
peak_720 = mean_spec_720(locs1);
xval_720peak = Lambda(locs1);
halfheight_720 = pks1/2;
text(xval_720peak,peak_720,'\leftarrow 728.2','FontSize',20);
%widths_720 = find(mean_spec_720 == halfheight_720);
%[~,~,idx]=unique(round(abs(a-n)),'stable');

peak_760 = mean_spec_760(locs2);
xval_760peak = Lambda(locs2);
halfheight_760 = pks2/2;
text(xval_760peak,peak_760,'\leftarrow 760.3','FontSize',20);

peak_800 = mean_spec_800(locs3);
xval_800peak = Lambda(locs3);
halfheight_800 = pks3/2;
text(xval_800peak,peak_800,'\leftarrow 802.7','FontSize',20);

peak_850 = mean_spec_850(locs4);
xval_850peak = Lambda(locs4);
halfheight_850 = pks4/2;
text(xval_850peak,peak_850,'\leftarrow 847.2','FontSize',20);

peak_890 = mean_spec_890(locs5);
xval_890peak = Lambda(locs5);
halfheight_890 = pks5/2;
text(xval_890peak,peak_890,'\leftarrow 884.0','FontSize',20);

hold on
line([xval_720peak xval_720peak], [0 peak_720],'LineWidth',2,'Color','black','LineStyle','--');
line([xval_760peak xval_760peak], [0 peak_760],'LineWidth',2,'Color','black','LineStyle','--');
line([xval_800peak xval_800peak], [0 peak_800],'LineWidth',2,'Color','black','LineStyle','--');
line([xval_850peak xval_850peak], [0 peak_850],'LineWidth',2,'Color','black','LineStyle','--');
line([xval_890peak xval_890peak], [0 peak_890],'LineWidth',2,'Color','black','LineStyle','--');
f=get(gca,'Children');
legend([f(15),f(14),f(13),f(12),f(11)],'720nm - FWHM = 23.9nm','760nm - FWHM = 30.3nm','800nm - FWHM = 34.1nm','850nm - FWHM = 38.6nm','890nm - FWHM = 54.5nm','Location','NorthWest')
% ylv = ylim;
% plot([1;1]*Lambda(locs1), ylv(:)*ones(size(pks1')))                                      % ... & Plot Vertical Lines Through Them
saveas(gcf,'LED_spectra_plot2.png')
hold off

% [pks,locs,w,p] = findpeaks(all_LED_plot,Lambda,'MinPeakHeight',500,'MinPeakWidth',15,'Annotate','extents');
% findpeaks(all_LED_plot,Lambda,'MinPeakHeight',500,'MinPeakWidth',15,'Annotate','extents')
% title('Signal Peak Widths')
%%

% Weight spectrum
int_weight = mean_spec(:,:)./mean_spec_int_720;

% plot weighted spectra
figure(3)
plot(Lambda,int_weight)
xlabel lambda
ylabel counts

% power repartition as function of wl (in units of A)
power_weight = int_weight * 0.000001125; % value from your power measurement in Amps
% plot weighted power
figure(4)
plot(Lambda,power_weight)
xlabel lambda
ylabel A

% interpolate to wavelengths
wl = data{1,1}.data(:,1);
wl_int = 400:5:1100;
power_weight_int = interp1(wl, power_weight, wl_int, 'spline');

%% import responsivity of powermeter
% responsivity in units A/W
% replace by the responsivity of your detector
resp = importdata('powermeter.txt',',');
resp_data = resp(1:end,2);

wl_int = wl_int.';
% plot calibration curve
figure(5)
plot(wl_int,resp_data)
xlabel 'wavelength (nm)'
ylabel 'Responsivity A/W'

%% power in W/nm
% get power in Watts as function of wavelength = power repartition /
% responsivity
power_watts = power_weight_int' ./ resp_data;
figure(6)
plot(wl_int,power_watts)
xlabel 'wavelength (nm)'
ylabel 'Power (W)'
TotalPower=sum(power_watts(1:121));
%TotalPower2=trapz(wl_int(1:121),power_watts(1:121)); test to see if sum if
%good enough

%save data

% fid = fopen(fileName,'w');
% %X = cat(2,rot90(wl_int),power_watts)
% X = [wl_int(:) power_watts(:)];
% writematrix(X, fileName);


%% inverion
% this was just to check if the value measured by powermeater in Watts was
% consistant with calibration
% PowerAt600=0.00016743./resp_data(41)
% PowerAt800=0.00016743./resp_data(81)
% PowerAt900=0.00016743./resp_data(101)
%% different boundries
% %power repartition as function of wl (in units of A)
% power_weight_threshold = int_weight * 0.00016743;
% % plot weighted power
% figure
% plot(Lambda,power_weight_threshold)
% xlabel lambda
% ylabel A
% 
% % interpolate to wavelengths
% wl = data{1,1}.data(:,1);
% wl_int = 400:5:1100;
% power_weight_int = interp1(wl, power_weight_threshold, wl_int, 'spline');
% 
% % responsivity in units A/W
% resp = importdata('C:\Users\rmapfla\Desktop\Power Measurements HL2000\CalibrationS130C_500mw.xlsx');
% resp_data = resp.data(:,2);
% % plot calibration curve
% figure
% plot(wl_int,resp_data)
% xlabel lambda
% ylabel A/W
% 
% % get power in Watts as function of wavelength = power repartition /
% % responsivity
% power_watts_threshold = power_weight_int' ./ resp_data(53:136);
% figure
% plot(wl_int(53:136),power_watts_threshold)
% xlabel lambda
% ylabel W
% TotalPower=sum(power_watts_threshold)






