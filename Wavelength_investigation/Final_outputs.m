%% Tissue reduced scattering coefficient finder
list = readmatrix('wav_650_900.xlsx');

tissue_name = ["skin","brain","breast","bone","soft tissue","fibrous tissue","fatty tissue","CSF"];
tissue_a = [46 24.2 16.8 22.9 18.9 27.1 18.4 12];
tissue_b = [1.421 1.611 1.055 0.716 1.286 1.627 0.672 0.8];

% Vector of wavelengths from 650 to 900nm in steps of 1nm
%lambda = linspace(650e-9,900e-9,251);
lambda = list(:,9);
%%
% Equation from Jacques review paper
for l = 1:length(lambda)
    for v = 1:length(tissue_a)
        tissue_reduced_mus(v,l) = tissue_a(v) * ((lambda(l)/(500))^(-tissue_b(v)));
    end
end
%% Extinction Spectra

% extinction coefficient in cm^-1 mM^-1
full_HbR = list(:,2); 
full_HbO = list(:,3);
full_CCO = interp1(list((1:66),4),list((1:66),5),list(:,9));

% absorption coefficient in cm^-1
full_water = interp1(list((1:30),7),list((1:30),8),list(:,9));
full_lipid = list(:,10);

% Convert molar extinction coefficients to absorption coefficients in mm^-1
% These are specifically for the brain (multiplication terms are for 
% brain concentrations) - 2.3 is log(e) factor
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
% Equation from Jem 2020 paper to find total absorption from all
% chromophores
k = 1;
for f = 1:length(lambda)
    mua_brain(f) = absorption_HbO(f) + absorption_HbR(f) + absorption_water(f)*W + absorption_lipid(f)*L + B;
    k = k+1;
end

mua_brain = mua_brain.';

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
%%
% Equation to find mu effective
for f = 1:length(lambda)
    mueff_brain(f) = sqrt(3*mua_brain(f)*(mua_brain(f)+sprime_brain_mm(f)));
    %k = k+1;
end

mueff_brain = mueff_brain.';
% Pathlength is inverse of mu effective
for f = 1:length(lambda)
    pathlength_brain(f) = 1/mueff_brain(f);
    %k = k+1;
end

pathlength_brain = pathlength_brain.';
%%
master_matrix = [lambda,full_HbO,full_HbR,full_CCO,full_water,full_lipid,absorption_HbO,absorption_HbR,absorption_CCO,absorption_water,absorption_lipid,sprime_brain_mm,mua_brain,mueff_brain,pathlength_brain];

%% User input age and gender of participant

prompt = 'How old is participant? (32 weeks, 34 weeks, term baby or adult): ';
answer2 = input(prompt,'s');

if strcmp(answer2,'32 weeks')
    agekey = 1;
elseif strcmp(answer2,'34 weeks')
    agekey = 2;
elseif strcmp(answer2,'term baby')
    agekey = 3;
elseif strcmp(answer2,'adult')
    agekey = 4;
end
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
% Effect of skin colour changes results
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

%% Bone

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

%% Dura Matter

reduced_scat_coeff_ft = tissue_reduced_mus(6,:);
sprime_ft_mm = reduced_scat_coeff_ft/10;
sprime_ft_mm = sprime_ft_mm.';

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

%% CSF

% csf_thickness = 3;
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
%% Working on it

thickness_adult_scalp = [7.43 6.3 8.98 9.08 9.45 7.54 8.02 7.5 7.04 6.27 6.11];
thickness_adult_scalp = thickness_adult_scalp.';
thickness_adult_skull = [5.01 5.88 6.29 2.72 2.97 4.32 4.13 4.5 5.63 5.62 5.07];
thickness_adult_skull = thickness_adult_skull.';
thickness_adult_csf = [2.07 2.07 3.22 2.07 2.19 2.52 7.12 3.51 4.2 1.97 2.52];
thickness_adult_csf = thickness_adult_csf.';
thickness_regions = ["common","Oz","Pz","T7","T8","C3","Cz","C4","Fz","Fp1","Fp2"];
thickness_regions = thickness_regions.';
thickness_adult = [thickness_adult_scalp,thickness_adult_skull,thickness_adult_csf];
%%
initial_intensity = 1;
thickness_epidermis = 0.05;

if agekey == 4
    thicknesses = thickness_adult;
end

%%


for f = 1:length(lambda)
    intensity_epidermis(f) = initial_intensity*(exp(-mueff_epidermis(f)*thickness_epidermis));
end

intensity_epidermis = intensity_epidermis.';

for f = 1:length(lambda)
    intensity_dermis(f) = initial_intensity*(exp(-mueff_dermis(f)*thickness_dermis));
end

intensity_dermis = intensity_dermis.';

for f = 1:length(lambda)
    intensity_bone(f) = initial_intensity*(exp(-mueff_bone(f)*skull_thickness));
end

intensity_bone = intensity_bone.';

% for f = 1:length(lambda)
%     intensity_ft(f) = initial_intensity*(exp(-mueff_ft(f)*dura_thickness));
% end
% 
% intensity_ft = intensity_ft.';
for f = 1:length(lambda)
    intensity_csf(f) = initial_intensity*(exp(-mueff_csf(f)*csf_thickness));
end

intensity_csf = intensity_csf.';

for f = 1:length(lambda)
    intensity_before_brain(f) = initial_intensity*intensity_epidermis(f)*intensity_dermis(f)*intensity_bone(f)*intensity_ft(f);
end

intensity_before_brain = intensity_before_brain.';


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



%% Reduced scattering coefficients all in one

reduced_scat_coeff_csf = tissue_reduced_mus(8,:);
sprime_csf_mm = reduced_scat_coeff_csf/10;
sprime_csf_mm = sprime_csf_mm.';

reduced_scats = [sprime_brain_mm,sprime_bone_mm,sprime_skin_mm,sprime_st_mm,sprime_csf_mm];
figure(24)
plot(lambda,reduced_scats, 'LineWidth',4)
ax = gca;
%ax.YColor = 'blue';
leftcolors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 1 0 1];
colororder(leftcolors)
xlabel('Wavelength (nm)') 
ylabel('Reduced Scattering Coefficient (mm^-^1)')
legend({'Brain', 'Bone', 'Skin', 'Soft tissue', 'CSF'},'Location','NorthEast','FontSize',23)
ax.FontSize = 23;
%xlim([785 915])


for f = 1:length(lambda)
    intensity_2(f) = intensity_epidermis(f)*intensity_dermis(f);
end
intensity_2 = intensity_2.';
for f = 1:length(lambda)
    intensity_3(f) = intensity_2(f)*intensity_bone(f);
end
intensity_3 = intensity_3.';
for f = 1:length(lambda)
    intensity_4(f) = intensity_3(f)*intensity_ft(f);
end
intensity_4 = intensity_4.';

intensity_accumulate = [intensity_epidermis,intensity_2,intensity_3,intensity_4];
figure(25)
plot(lambda,intensity_accumulate, 'LineWidth',4)
ax = gca;
%ax.YColor = 'blue';
leftcolors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 1 0 1];
colororder(leftcolors)
xlabel('Wavelength (nm)') 
ylabel('Relative Intensity')
legend({'Epidermis', 'Epidermis + Dermis', 'Epidermis + Dermis + Skull', 'Epidermis + Dermis + Skull + Dura'},'Location','NorthEast','FontSize',23)
ax.FontSize = 23;

%%

