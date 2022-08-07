%% Tissue reduced scattering coefficient finder
list = readmatrix('wav_650_900.xlsx');

tissue_name = ["skin","brain","breast","bone","soft tissue","fibrous tissue","fatty tissue"];
tissue_a = [46 24.2 16.8 22.9 18.9 27.1 18.4];
tissue_b = [1.421 1.611 1.055 0.716 1.286 1.627 0.672];

% Vector of wavelengths from 650 to 900nm in steps of 1nm
%lambda = linspace(650e-9,900e-9,251);

Gemma_exts = tissue_specific_extinction_coefficient_650to1000;
lambda = Gemma_exts(:,1);
%%
%l = 1;
for l = 1:length(lambda)
    for v = 1:length(tissue_a)
        tissue_reduced_mus(v,l) = tissue_a(v) * ((lambda(l)/(500))^(-tissue_b(v)));
    end
    %l = l+1;
end

% %%
% %tissue_matrix_names = [tissue_name;tissue_a;tissue_b;tissue_reduced_mus];
% 
% 
% %%Plot values of reduced scattering coefficient at user specified
% %%wavelenght
% 
% prompt1 = 'What is the desired tissue type? (skin, brain, breast, bone, soft tissue, fibrous tissue, fatty tissue): ';
% answer1 = input(prompt1,'s');
% 
% if strcmp(answer1,'skin')
%     scat_val = tissue_reduced_mus(1,:);
% elseif strcmp(answer1,'brain')
%     scat_val = tissue_reduced_mus(2,:);
% elseif strcmp(answer1,'breast')
%     scat_val = tissue_reduced_mus(3,:);
% elseif strcmp(answer1,'bone')
%     scat_val = tissue_reduced_mus(4,:);
% elseif strcmp(answer1,'soft tissue')
%     scat_val = tissue_reduced_mus(5,:);
% elseif strcmp(answer1,'fibrous tissue')
%     scat_val = tissue_reduced_mus(6,:);
% elseif strcmp(answer1,'fatty tissue')
%     scat_val = tissue_reduced_mus(7,:);
% end
% 
% %%
% 
% figure(1)
% plot(lambda,scat_val,'marker','.','MarkerSize',15)
% xlabel('Wavelength (nm)')
% ylabel('Reduced scattering coefficient (cm^-1)')
% title('Graph of reduced scattering coefficient against wavelength for',answer1)

%% Extinction Spectra
% list = readmatrix('extinctions.xlsx');
%list = readmatrix('wav_970.xlsx');
%list = readmatrix('vals_650.xlsx');
%list = readmatrix('wav_650_900.xlsx');

lambda = Gemma_exts(:,1);

full_HbR = Gemma_exts(:,4);
full_HbO = Gemma_exts(:,3);
full_CCO = Gemma_exts(:,5);
full_water = Gemma_exts(:,2);
%full_water = interp1(list((1:30),7),list((1:30),8),lambda);
full_lipid = Gemma_exts(:,9);
full_HbR = full_HbR/10000;
full_HbO = full_HbO/10000;
full_CCO = full_CCO/10000;
full_water = full_water/10000;
full_lipid = full_lipid/10000;

absorption_HbR = full_HbR*2.303;
absorption_HbO = full_HbO*2.303;
absorption_CCO = full_CCO*2.303;
absorption_water = full_water*2.303;
absorption_lipid = full_lipid*2.303;
% %%
% combos = nchoosek(1:length(list),5);
% 
% for l = 1:length(combos)
%     A{l} = [list(combos(l,1),2) list(combos(l,1),3) list(combos(l,1),4); list(combos(l,2),2) list(combos(l,2),3) list(combos(l,2),4); list(combos(l,3),2) list(combos(l,3),3) list(combos(l,3),4); list(combos(l,4),2) list(combos(l,4),3) list(combos(l,4),4); list(combos(l,5),2) list(combos(l,5),3) list(combos(l,5),4)];
% end
% 
% grad_HbO = gradient(list(:,3));
% grad_HbR = gradient(list(:,2));
% grad_CCO = gradient(list(:,4));
% 
% grad_plot = [grad_HbO,grad_HbR,grad_CCO];
% 
% for l = 1:length(list)
%     total_grad(l,1) = abs(grad_HbO(l)) + abs(grad_HbR(l)) + abs(grad_CCO(l));
% end
% 
% slope_list = [total_grad,list(:,1)];
% 
% for l = 1:length(list)
%     dist_HbR_HbO(l,1) = abs(list(l,3) - list(l,2));
%     dist_HbR_CCO(l,1) = abs(list(l,2) - list(l,4));
%     dist_CCO_HbO(l,1) = abs(list(l,3) - list(l,4));
%     total_dist(l,1) = dist_HbR_HbO(l,1) + dist_HbR_CCO(l,1) + dist_CCO_HbO(l,1);
% end
% 
% for i = 1:length(A)
%     condition_number(i,:) = cond(A{i});
% end
% 
% % for i = 1:length(A)
% %     determinant(i,:) = det(A{i});
% % end
% 
% N = min(condition_number);
% found = find(condition_number==N);
% 
% indices = combos(found,:);
% best_wavelengths = list(indices,1);
% %det_min = determinant(found);
% 
% for k = 1:length(A)
%     wavelength_inds{k} = combos(k,:);
%     wavelength_list{k} = list(wavelength_inds{k},1);
% end
% wav_list2 = cell2mat(wavelength_list);
% wav_official = wav_list2.';
% 
% xvals = 1:length(condition_number);
% xvals2 = reshape(xvals,[length(condition_number),1]);
% 
% master_list = [condition_number,wav_official];
% %master_list_sorted = sort(master_list,1);
% 
% [~,idx] = sort(master_list(:,1)); % sort just the first column
% master_list_sorted = master_list(idx,:);
% 
% spectra_HbR = list(:,2);
% spectra_HbO = list(:,3);
% spectra_CCO = list(:,4);

%% Absorption Coefficient

%values here from Jem's 2020 paper
conc_HbO = 0.056; %in mMol
conc_HbR = 0.024;

B = 0.012; %in mm^-1

W = 0.8; %constants
L = 0.116;

k = 1;
for f = 1:length(lambda)
    mua_brain(f) = full_HbO(f)*conc_HbO + full_HbR(f)*conc_HbR + absorption_water(f)*W + absorption_lipid(f)*L + B;
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
plot(lambda,pathlength_brain, 'LineWidth',2)
xlabel('Wavelength (nm)')
ylabel('Mean free pathlength (mm)')
title('Graph showing mean free pathlength of infrared light in brain tissue')
ax = gca;
ax.FontSize = 16;

%%

absorptions_list = [absorption_HbR,absorption_HbO,absorption_CCO,absorption_water,absorption_lipid];
figure(3)
plot(lambda,absorptions_list, 'LineWidth',3)
xlabel('Wavelength (nm)')
ylabel('Absorption Coefficient (mm^-1 mM^-1)')
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

% %%
% figure(4)
% plot(lambda,sprime_brain_mm, 'LineWidth',3)
% xlabel('Wavelength (nm)')
% ylabel('Reduced scattering coefficient (mm^-1)')
% title('Graph showing relationship between reduced scattering coefficient of the brain and wavelength')
% ax = gca;
% ax.FontSize = 23;
% %%
% 
% figure(5)
% plot(lambda,sprime_bone_mm, 'LineWidth',3)
% xlabel('Wavelength (nm)')
% ylabel('Reduced scattering coefficient (mm^-1)')
% title('Graph showing relationship between reduced scattering coefficient of bone and wavelength')
% ax = gca;
% ax.FontSize = 25;
% 
% figure(6)
% plot(lambda,sprime_skin_mm, 'LineWidth',3)
% xlabel('Wavelength (nm)')
% ylabel('Reduced scattering coefficient (mm^-1)')
% title('Graph showing relationship between reduced scattering coefficient of skin and wavelength')
% ax = gca;
% ax.FontSize = 25;
% 
% figure(7)
% plot(lambda,sprime_st_mm, 'LineWidth',3)
% xlabel('Wavelength (nm)')
% ylabel('Reduced scattering coefficient (mm^-1)')
% title('Graph showing relationship between reduced scattering coefficient of soft tissue and wavelength')
% ax = gca;
% ax.FontSize = 23;