%% Determinant Finder 3 wavelengths

% list = readmatrix('extinctions.xlsx');
%list = readmatrix('wav_970.xlsx');
%list = readmatrix('vals_650.xlsx');
list = readmatrix('wav_650_900.xlsx');

combos = nchoosek(1:length(list),5);
%%

for l = 1:length(combos)
    A{l} = [list(combos(l,1),2) list(combos(l,1),3) list(combos(l,1),4); list(combos(l,2),2) list(combos(l,2),3) list(combos(l,2),4); list(combos(l,3),2) list(combos(l,3),3) list(combos(l,3),4); list(combos(l,4),2) list(combos(l,4),3) list(combos(l,4),4); list(combos(l,5),2) list(combos(l,5),3) list(combos(l,5),4)];
end

%%

for i = 1:length(A)
    condition_number(i,:) = cond(A{i});
end

% for i = 1:length(A)
%     determinant(i,:) = det(A{i});
% end

%%
N = min(condition_number);
found = find(condition_number==N);

indices = combos(found,:);
best_wavelengths = list(indices,1);
%det_min = determinant(found);

%%
for k = 1:length(A)
    wavelength_inds{k} = combos(k,:);
    wavelength_list{k} = list(wavelength_inds{k},1);
end
wav_list2 = cell2mat(wavelength_list);
wav_official = wav_list2.';
% %%
% find_wavelengths = find(combos(:,1)==13 & combos(:,2)==21 & combos(:,3)==39);
% condition_3lambda = condition_number(find_wavelengths);
% det_3lamb = determinant(find_wavelengths);
% %%
% find_wavelengths2 = find(combos(:,1)==8 & combos(:,2)==23 & combos(:,3)==31);
% condition_3lambda2 = condition_number(find_wavelengths2);
% det_3lamb2 = determinant(find_wavelengths2);
%%
xvals = 1:length(condition_number);
xvals2 = reshape(xvals,[length(condition_number),1]);
%%
master_list = [condition_number,wav_official];
%master_list_sorted = sort(master_list,1);

[~,idx] = sort(master_list(:,1)); % sort just the first column
master_list_sorted = master_list(idx,:);
%% Uncomment if want to plot first 50 values
% conds_toplot = master_list_sorted(1:50,1);
% wav1_plot = master_list_sorted(1:50,2);
% wav2_plot = master_list_sorted(1:50,3);
% wav3_plot = master_list_sorted(1:50,4);
% wav4_plot = master_list_sorted(1:50,5);
% wav5_plot = master_list_sorted(1:50,6);
% wavelengths_plot = [wav1_plot,wav2_plot,wav3_plot,wav4_plot,wav5_plot];

%% Uncomment if want to plot first 500 values
conds_toplot = master_list_sorted(1:500,1);
wav1_plot = master_list_sorted(1:500,2);
wav2_plot = master_list_sorted(1:500,3);
wav3_plot = master_list_sorted(1:500,4);
wav4_plot = master_list_sorted(1:500,5);
wav5_plot = master_list_sorted(1:500,6);
wavelengths_plot = [wav1_plot,wav2_plot,wav3_plot,wav4_plot,wav5_plot];

%% Uncomment if want to plot all values (minus the extremes)
% conds_toplot = master_list_sorted(1:length(combos)-30,1);
% wav1_plot = master_list_sorted(1:length(combos)-30,2);
% wav2_plot = master_list_sorted(1:length(combos)-30,3);
% wav3_plot = master_list_sorted(1:length(combos)-30,4);
% wav4_plot = master_list_sorted(1:length(combos)-30,5);
% wav5_plot = master_list_sorted(1:length(combos)-30,6);
% wavelengths_plot = [wav1_plot,wav2_plot,wav3_plot,wav4_plot,wav5_plot];
new_plot = [list(:,2),list(:,3),list(:,4)];
%% Uncomment if want to plot all values (including the extremes)
% conds_toplot = master_list_sorted(1:length(combos),1);
% wav1_plot = master_list_sorted(1:length(combos),3);
% wav2_plot = master_list_sorted(1:length(combos),4);
% wav3_plot = master_list_sorted(1:length(combos),5);
% wavelengths_plot = [wav1_plot,wav2_plot,wav3_plot];
%%
newcolors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840];
colororder(newcolors)
figure(1)
%h = yyaxis(wavelengths_plot,conds_toplot,list(:,1),new_plot,'plot');
yyaxis left 
plot(wavelengths_plot,conds_toplot,'LineStyle', 'none' ,'marker','.', 'MarkerSize',30)
ax = gca;
ax.YColor = 'blue';
leftcolors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 1 0 1];
colororder(leftcolors)
xlabel('Wavelength (nm)') 
ylabel('Condition Number')
xlim([650 905])

yyaxis right
plot(list(:,1),new_plot, 'LineWidth',1.5)
rightcolors = [0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840];
colororder(rightcolors)
ylabel('Extinction Coefficient')

%'LineStyle', 'none' ,'marker','.', 'MarkerSize',30);
%xlim([0 200000])
%ylim([650 1000])
title({'Graph showing different five-wavelength combinations and their corresponding condition number','alongside extinction coefficients of chromophores'},'FontSize',16)
legend({'Wavelength 1', 'Wavelength 2', 'Wavelength 3', 'Wavelength 4', 'Wavelength 5', 'HbO', 'HbR', 'oxCCO'},'Location','southwest','FontSize',16)

ax.FontSize = 16;

%hold on


%figure(2)
%plot(list(:,1),new_plot)
% stem(conds_toplot,wavelengths_plot, 'LineWidth',1)

%%
grad_HbO = gradient(list(:,3));
grad_HbR = gradient(list(:,2));
grad_CCO = gradient(list(:,4));
grads = [grad_HbO,grad_HbR,grad_CCO];

figure(2)
%h = yyaxis(wavelengths_plot,conds_toplot,list(:,1),new_plot,'plot');
yyaxis left 
plot(wavelengths_plot,conds_toplot,'LineStyle', 'none' ,'marker','.', 'MarkerSize',30)
ax = gca;
ax.YColor = 'blue';
leftcolors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 1 0 1];
colororder(leftcolors)
xlabel('Wavelength (nm)') 
ylabel('Condition Number')
xlim([650 905])

yyaxis right
plot(list(:,1),grads, 'LineWidth',1.5)
rightcolors = [0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840];
colororder(rightcolors)
ylabel('Slope of Extinction Coefficient')

%'LineStyle', 'none' ,'marker','.', 'MarkerSize',30);
%xlim([0 200000])
%ylim([650 1000])
title({'Graph showing different five-wavelength combinations and their corresponding condition number','alongside the slope of the extinction coefficients of chromophores'},'FontSize',16)
legend({'Wavelength 1', 'Wavelength 2', 'Wavelength 3', 'Wavelength 4', 'Wavelength 5', 'HbO', 'HbR', 'oxCCO'},'Location','southwest','FontSize',16)

ax.FontSize = 16;


%%

for l = 1:length(list)
    total_grad(l,1) = (grad_HbO(l)) + (grad_HbR(l)) + (grad_CCO(l));
end

slope_list = [total_grad,list(:,1)];

figure(3)
%h = yyaxis(wavelengths_plot,conds_toplot,list(:,1),new_plot,'plot');
yyaxis left 
plot(wavelengths_plot,conds_toplot,'LineStyle', 'none' ,'marker','.', 'MarkerSize',30)
ax = gca;
ax.YColor = 'blue';
leftcolors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 1 0 1];
colororder(leftcolors)
xlabel('Wavelength (nm)') 
ylabel('Condition Number')
xlim([650 905])

yyaxis right
plot(list(:,1),total_grad, 'LineWidth',1.5)
rightcolors = [0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840];
colororder(rightcolors)
ylabel('Slope of Extinction Coefficient')

%'LineStyle', 'none' ,'marker','.', 'MarkerSize',30);
%xlim([0 200000])
%ylim([650 1000])
title({'Graph showing different five-wavelength combinations and their corresponding condition number','alongside the sum of the slopes of the extinction coefficients of chromophores'},'FontSize',16)
legend({'Wavelength 1', 'Wavelength 2', 'Wavelength 3', 'Wavelength 4', 'Wavelength 5', 'Sum of Sloped', 'HbR', 'oxCCO'},'Location','southwest','FontSize',16)

ax.FontSize = 16;


%%

% new_plot = [list(:,2),list(:,3),list(:,4)];
% figure(2)
% plot(list(:,1),new_plot)

%% Slope investigations
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
% figure(3)
% plot(list(:,1),total_grad)
% 
% figure(4)
% plot(list(:,1),grad_plot)