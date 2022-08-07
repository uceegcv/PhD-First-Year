%% Complete brain intensity calculator using Rob's paper
list = readmatrix('wav_650_900.xlsx');

tissue_name = ["skin","brain","breast","bone","soft tissue","fibrous tissue","fatty tissue","csf"];
tissue_a = [46 24.2 16.8 22.9 18.9 27.1 18.4 12];
tissue_b = [1.421 1.611 1.055 0.716 1.286 1.627 0.672 0.8];

% Vector of wavelengths from 650 to 900nm in steps of 1nm
%lambda = linspace(650e-9,900e-9,251);
lambda = list(:,9);
%%
%l = 1;
for l = 1:length(lambda)
    for v = 1:length(tissue_a)
        tissue_reduced_mus(v,l) = tissue_a(v) * ((lambda(l)/(500))^(-tissue_b(v)));
    end
    %l = l+1;
end

%% Extinction Spectra

% extinction coefficient in cm^-1 mM^-1
full_HbR = list(:,2); 
full_HbO = list(:,3);
full_CCO = interp1(list((1:66),4),list((1:66),5),list(:,9));

% absorption coefficient in cm^-1
full_water = interp1(list((1:30),7),list((1:30),8),list(:,9));
full_lipid = list(:,10);

absorption_HbR = (full_HbR*0.0024)/2.303;
absorption_HbO = (full_HbO*0.0056)/2.303;
absorption_CCO = (full_CCO*0.0003)/2.303;
absorption_water = full_water/10;
absorption_lipid = full_lipid/10;

%% Absorption Coefficient

%values here from Jem's 2020 paper
conc_HbO = 0.056; %in mMol
conc_HbR = 0.024;

B = 0.012; %in mm^-1

W = 0.8; %constants
L = 0.116;

k = 1;
for f = 1:length(lambda)
    mua_brain(f) = absorption_HbO(f) + absorption_HbR(f) + absorption_water(f)*W + absorption_lipid(f)*L + B;
    k = k+1;
end

mua_brain = mua_brain.';
%%
figure(1)
plot(lambda,mua_brain, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Absorption Coefficient (mm^-1)')
title('Graph showing absorption coefficient of brain tissue against wavelength')
ax = gca;
ax.FontSize = 25;

%%

reduced_scat_coeff_brain = tissue_reduced_mus(2,:);
sprime_brain_mm = reduced_scat_coeff_brain/10;
sprime_brain_mm = sprime_brain_mm.';

reduced_scat_coeff_skin = tissue_reduced_mus(1,:);
sprime_skin_mm = reduced_scat_coeff_skin/10;
sprime_skin_mm = sprime_skin_mm.';

reduced_scat_coeff_bone = tissue_reduced_mus(4,:);
sprime_bone_mm = reduced_scat_coeff_bone/10;
sprime_bone_mm = sprime_bone_mm.';

reduced_scat_coeff_st = tissue_reduced_mus(5,:);
sprime_st_mm = reduced_scat_coeff_st/10;
sprime_st_mm = sprime_st_mm.';

reduced_scat_coeff_csf = tissue_reduced_mus(8,:);
sprime_csf_mm = reduced_scat_coeff_csf/10;
sprime_csf_mm = sprime_csf_mm.';
%%
for f = 1:length(lambda)
    mueff_brain(f) = sqrt(3*mua_brain(f)*(mua_brain(f)+sprime_brain_mm(f)));
    %k = k+1;
end

mueff_brain = mueff_brain.';

for f = 1:length(lambda)
    pathlength_brain(f) = 1/mueff_brain(f);
    %k = k+1;
end

pathlength_brain = pathlength_brain.';
%%
master_matrix = [lambda,full_HbO,full_HbR,full_CCO,full_water,full_lipid,absorption_HbO,absorption_HbR,absorption_CCO,absorption_water,absorption_lipid,sprime_brain_mm,mua_brain,mueff_brain,pathlength_brain];
%%
figure(2)
plot(lambda,pathlength_brain, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Mean free pathlength (mm)')
title('Graph showing mean free pathlength of infrared light in brain tissue')
ax = gca;
ax.FontSize = 25;

%%

absorptions_list = [absorption_HbR,absorption_HbO,absorption_CCO,absorption_water,absorption_lipid];
figure(3)
plot(lambda,absorptions_list, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Absorption Coefficient (mm^-1)')
title('Graph showing absorption coefficient of different chromophores in the brain')
legend({'HbR','HbO','CCO','water','lipid'})
ax = gca;
ax.FontSize = 25;

%%
% absorptions_list_2 = [Gemma_exts(:,4),Gemma_exts(:,3),Gemma_exts(:,5),Gemma_exts(:,2),Gemma_exts(:,9)];
% figure(3)
% plot(Gemma_exts(:,1),absorptions_list_2, 'LineWidth',3)
% xlabel('Wavelength (nm)')
% ylabel('Absorption Coefficient (mm^-1 mM^-1)')
% title('Graph showing absorption coefficient of different chromophores in the brain')
% legend({'HbR','HbO','CCO','water','lipid'})
% ax = gca;
% ax.FontSize = 25;

%%
figure(4)
plot(lambda,sprime_brain_mm, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Reduced scattering coefficient (mm^-1)')
title('Graph showing relationship between reduced scattering coefficient of the brain and wavelength')
ax = gca;
ax.FontSize = 23;
%%

figure(5)
plot(lambda,sprime_bone_mm, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Reduced scattering coefficient (mm^-1)')
title('Graph showing relationship between reduced scattering coefficient of bone and wavelength')
ax = gca;
ax.FontSize = 25;
% 
figure(6)
plot(lambda,sprime_skin_mm, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Reduced scattering coefficient (mm^-1)')
title('Graph showing relationship between reduced scattering coefficient of skin and wavelength')
ax = gca;
ax.FontSize = 25;
% 
% figure(7)
% plot(lambda,sprime_st_mm, 'LineWidth',3)
% xlabel('Wavelength (nm)')
% ylabel('Reduced scattering coefficient (mm^-1)')
% title('Graph showing relationship between reduced scattering coefficient of soft tissue and wavelength')
% ax = gca;
% ax.FontSize = 23;

%% Skin 

for f = 1:length(lambda)
    mua_skin_baseline(f) = (7.84*10^8)*(lambda(f)^(-3.255));
end

mua_skin_baseline = mua_skin_baseline.';
mua_skin_baseline = mua_skin_baseline/10;

for f = 1:length(lambda)
    mua_melanin(f) = (6.6*10^11)*(lambda(f)^(-3.33));
end

mua_melanin = mua_melanin.';
mua_melanin = mua_melanin/10;

prompt = 'Select skin colour (light, medium and dark): ';
answer1 = input(prompt,'s');

if strcmp(answer1,'light')
    mel_conc = 0.04;
elseif strcmp(answer1,'medium')
    mel_conc = 0.13;
elseif strcmp(answer1,'dark')
    mel_conc = 0.3;
end

absorption_HbR_skin = (full_HbR*0.000065)/2.303;
absorption_HbO_skin = (full_HbO*0.000065)/2.303;
SO2_skin = 0.71;
water_conc_skin = 0.6;

% for f = 1:length(lambda)
%     total_mua_skin(f) = (mua_melanin(f)*mel_conc) + (absorption_HbO_skin(f)*SO2_skin) + (absorption_HbR_skin(f)*(1-SO2_skin)) + (absorption_water(f)*water_conc_skin);
% end
for f = 1:length(lambda)
    mua_epidermis(f) = (mua_melanin(f)*mel_conc) + (absorption_HbO(f)*SO2_skin) + (absorption_HbR(f)*(1-SO2_skin)) + (absorption_water(f)*water_conc_skin);
end

for f = 1:length(lambda)
    mua_dermis(f) = (absorption_HbO(f)*SO2_skin) + (absorption_HbR(f)*(1-SO2_skin)) + (absorption_water(f)*water_conc_skin);
end

% for f = 1:length(lambda)
%     total_mua_dermis(f) = (mua_melanin(f)*mel_conc) + (absorption_HbO(f)*SO2_skin) + (absorption_HbR(f)*(1-SO2_skin)) + (absorption_water(f)*water_conc_skin);
% end
mua_epidermis = mua_epidermis.';
mua_dermis = mua_dermis.';

for f = 1:length(lambda)
    mueff_epidermis(f) = sqrt(3*mua_epidermis(f)*(mua_epidermis(f)+sprime_skin_mm(f)));
end

for f = 1:length(lambda)
    mueff_dermis(f) = sqrt(3*mua_dermis(f)*(mua_dermis(f)+sprime_skin_mm(f)));
end

mueff_epidermis = mueff_epidermis.'; 
mueff_dermis = mueff_dermis.';

for f = 1:length(lambda)
    pathlength_epidermis(f) = 1/mueff_epidermis(f);
end

pathlength_epidermis = pathlength_epidermis.';

for f = 1:length(lambda)
    pathlength_dermis(f) = 1/mueff_dermis(f);
end

pathlength_dermis = pathlength_dermis.';

figure(7)
plot(lambda,pathlength_epidermis, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Mean free pathlength (mm)')
title(sprintf('Graph showing mean free pathlength of infrared light in %s skin epidermis tissue', answer1))
ax = gca;
ax.FontSize = 25;

figure(8)
plot(lambda,pathlength_dermis, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Mean free pathlength (mm)')
title('Graph showing mean free pathlength of infrared light in dermis tissue')
ax = gca;
ax.FontSize = 25;
%%
figure(9)
plot(lambda,mua_melanin, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Absorption (mm^-1)')
title('Graph showing absorption coefficient of melanin with wavelength')
ax = gca;
ax.FontSize = 25;
%%
thickness_epidermis = 0.05;
thickness_dermis = 7.5;
initial_intensity = 1;

for f = 1:length(lambda)
    intensity_epidermis(f) = initial_intensity*(exp(-mueff_epidermis(f)*thickness_epidermis));
end

intensity_epidermis = intensity_epidermis.';

for f = 1:length(lambda)
    intensity_dermis(f) = initial_intensity*(exp(-mueff_dermis(f)*thickness_dermis));
end

intensity_dermis = intensity_dermis.';

figure(10)
plot(lambda,intensity_epidermis, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Relative output intensity')
title(['Graph showing relative output light intensity against wavelength for epidermis of ',answer1, ' skin'])
ax = gca;
ax.FontSize = 25;

figure(11)
plot(lambda,intensity_dermis, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Relative output intensity')
title('Graph showing relative output light intensity against wavelength for dermis')
ax = gca;
ax.FontSize = 25;

%% Bone

skull_thickness = 3;
% in cm^-1
absorption_bone = interp1(list((1:13),11),list((1:13),12),list(:,9));
% in mm^-1
absorption_bone = absorption_bone/10;

for f = 1:length(lambda)
    mueff_bone(f) = sqrt(3*absorption_bone(f)*(absorption_bone(f)+sprime_bone_mm(f)));
end

mueff_bone = mueff_bone.';

for f = 1:length(lambda)
    pathlength_bone(f) = 1/mueff_bone(f);
end

pathlength_bone = pathlength_bone.';

for f = 1:length(lambda)
    intensity_bone(f) = initial_intensity*(exp(-mueff_bone(f)*skull_thickness));
end

intensity_bone = intensity_bone.';

%% Dura Matter

dura_thickness = 0.7;
reduced_scat_coeff_ft = tissue_reduced_mus(6,:);
sprime_ft_mm = reduced_scat_coeff_ft/10;
sprime_ft_mm = sprime_ft_mm.';

figure(12)
plot(lambda,sprime_ft_mm, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Reduced scattering coefficient (mm^-1)')
title('Graph showing relationship between reduced scattering coefficient of fibrous tissue and wavelength')
ax = gca;
ax.FontSize = 23;

for f = 1:length(lambda)
    mua_ft(f) = 0.03;
end

mua_ft = mua_ft.';

for f = 1:length(lambda)
    mueff_ft(f) = sqrt(3*mua_ft(f)*(mua_ft(f)+sprime_ft_mm(f)));
end

mueff_ft = mueff_ft.';

for f = 1:length(lambda)
    pathlength_ft(f) = 1/mueff_ft(f);
end

pathlength_ft = pathlength_ft.';

for f = 1:length(lambda)
    intensity_ft(f) = initial_intensity*(exp(-mueff_ft(f)*dura_thickness));
end

intensity_ft = intensity_ft.';

%% CSF

csf_thickness = 3;
% in cm^-1
absorption_csf = absorption_water;
% in mm^-1
%absorption_bone = absorption_bone/10;

for f = 1:length(lambda)
    mueff_csf(f) = sqrt(3*absorption_csf(f)*(absorption_csf(f)+sprime_csf_mm(f)));
end

mueff_csf = mueff_csf.';

for f = 1:length(lambda)
    pathlength_csf(f) = 1/mueff_csf(f);
end

pathlength_csf = pathlength_csf.';

for f = 1:length(lambda)
    intensity_csf(f) = initial_intensity*(exp(-mueff_csf(f)*csf_thickness));
end

intensity_csf = intensity_csf.';

%% Total Intensity

for f = 1:length(lambda)
    intensity_before_brain(f) = initial_intensity*intensity_epidermis(f)*intensity_dermis(f)*intensity_bone(f)*intensity_ft(f)*intensity_csf(f);
end

intensity_before_brain = intensity_before_brain.';

figure(13)
plot(lambda,intensity_before_brain, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Relative Intensity')
title(['Graph showing relative intensity of light reaching the brain at different wavelengths for ', answer1,' skin'])
ax = gca;
ax.FontSize = 23;

%%
figure(14)
plot(lambda,absorptions_list, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Absorption Coefficient (mm^-1)')
xlim([650 900])
%title('Graph showing absorption coefficient of different chromophores in the brain')
legend({'HbR','HbO','CCO','water','lipid'})
ax = gca;
ax.FontSize = 25;
%%
zoomed_absorptions = [absorption_lipid,absorption_CCO];
figure(15)
plot(lambda,absorption_lipid, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Absorption Coefficient (mm^-1)')
xlim([650 900])
%title('Graph showing absorption coefficient of different chromophores in the brain')
legend({'lipid'})
ax = gca;
ax.FontSize = 25;

figure(16)
plot(lambda,absorption_CCO, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Absorption Coefficient (mm^-1)')
xlim([650 900])
%title('Graph showing absorption coefficient of different chromophores in the brain')
legend({'CCO'})
ax = gca;
ax.FontSize = 25;
%%
HbO_HbR_CCO = [full_HbO,full_HbR,full_CCO];
figure(17)
plot(lambda,HbO_HbR_CCO, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Extinction Coefficient (cm^-1 mM^-1)')
%xlim([650 900])
title('Graph showing extinction coefficients of HbO, HbR and CCCO')
legend({'HbO','HbR','CCO'})
ax = gca;
ax.FontSize = 25;

%%
water_brain = absorption_water*W;
lipid_brain = absorption_lipid*L;
absorptions_brain = [absorption_HbO,absorption_HbR,absorption_CCO,water_brain,lipid_brain];

figure(18)
plot(lambda,absorptions_brain, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Absorption Coefficient (mm^-1)')
%xlim([650 900])
title('Graph showing absorption coefficient of different chromophores in the brain')
legend({'HbO','HbR','CCO','water','lipid'},'Location','NorthWest')
ax = gca;
ax.FontSize = 25;

%%
water_brain = absorption_water*W;
lipid_brain = absorption_lipid*L;
absorptions_brain = [absorption_HbO,absorption_HbR,absorption_CCO,water_brain,lipid_brain];

figure(18)
plot(lambda,absorptions_brain, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Absorption Coefficient (mm^-1)')
%xlim([650 900])
title('Graph showing absorption coefficient of different chromophores in the brain')
legend({'HbO','HbR','CCO','water','lipid'},'Location','NorthWest')
ax = gca;
ax.FontSize = 25;

%%

% 
% for f = 1:length(lambda)
%     mua_dermis(f) = (absorption_HbO(f)*SO2_skin) + (absorption_HbR(f)*(1-SO2_skin)) + (absorption_water(f)*water_conc_skin);
% end

absorption_HbO_epidermis = absorption_HbO_skin*SO2_skin;
absorption_HbR_epidermis = absorption_HbR_skin*(1-SO2_skin);
absorption_melanin_epidermis = mua_melanin*mel_conc;
absorption_water_epidermis = absorption_water*water_conc_skin;


absorptions_epidermis = [absorption_HbO_epidermis,absorption_HbR_epidermis,absorption_melanin_epidermis,absorption_water_epidermis];

figure(19)
plot(lambda,absorptions_epidermis, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Absorption Coefficient (mm^-1)')
%xlim([650 900])
title(['Graph showing absorption coefficient of different chromophores in the epidermis of ', answer1, ' skin'])
legend({'HbO','HbR','melanin','water'},'Location','NorthEast')
ax = gca;
ax.FontSize = 23;


absorption_HbO_dermis = absorption_HbO_skin*SO2_skin;
absorption_HbR_dermis = absorption_HbR_skin*(1-SO2_skin);
absorption_water_dermis = absorption_water*water_conc_skin;

absorptions_dermis = [absorption_HbO_epidermis,absorption_HbR_epidermis,absorption_water_epidermis];

figure(20)
plot(lambda,absorptions_dermis, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Absorption Coefficient (mm^-1)')
%xlim([650 900])
title('Graph showing absorption coefficient of different chromophores in the dermis')
legend({'HbO','HbR','water'},'Location','NorthWest')
ax = gca;
ax.FontSize = 25;

%%
list_for_pdf = [lambda,full_HbO,full_HbR,full_CCO,absorption_water,absorption_lipid,absorption_HbO,absorption_HbR,absorption_CCO];
list_for_pdf2 = [lambda,mua_brain,mua_epidermis,mua_dermis,absorption_bone,mua_ft,intensity_before_brain];

just_water_lipid = [absorption_water,absorption_lipid];
figure(21)
plot(lambda,just_water_lipid, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Absorption Coefficient (mm^-1)')
xlim([700 820])
title('Graph showing absorption coefficient of water and lipid')
legend({'water','lipid'},'Location','NorthWest')
ax = gca;
ax.FontSize = 25;

%%
figure(22)
plot(lambda,intensity_bone, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Relative Intensity')
title('Graph showing relative intensity of light passing through the skull')
ax = gca;
ax.FontSize = 23;

figure(23)
plot(lambda,intensity_ft, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Relative Intensity')
title('Graph showing relative intensity of light passing through extra-cerebral tissues')
ax = gca;
ax.FontSize = 23;
%%
figure(24)
plot(lambda,intensity_csf, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Relative output intensity')
title('Graph showing relative output light intensity against wavelength for csf')
ax = gca;
ax.FontSize = 25;
%%
figure(25)
plot(lambda,sprime_csf_mm, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Reduced scattering coefficient (mm^-1)')
title('Graph showing relationship between reduced scattering coefficient of csf and wavelength')
ax = gca;
ax.FontSize = 25;