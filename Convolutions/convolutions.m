list = readmatrix('Conv_values.xlsx');

pin_720 = 0:0.1:3.99;

pin_720 = pin_720.';

artif_row = 0:0.01:3.99;
artif_row = artif_row.';

new_list = [list, artif_row];

%% 720nm pin

[Lia,Locb] = ismember(list(:,1),list(:,4));

[row,col,V] = find(Lia);

for f = 1:length(row)
    condvals_720(f) = list(row(f),3);
end

condvals_720 = condvals_720.';

takezero_720 = condvals_720(2:end);
otherhalf_720 = flip(takezero_720);

total_720 = [otherhalf_720 ; condvals_720];
%total_720 = total_720.';
%% 760nm pin

[Lia,Locb] = ismember(list(:,1),list(:,5));

[row,col,V] = find(Lia);

for f = 1:length(row)
    condvals_760(f) = list(row(f),3);
end

condvals_760 = condvals_760.';
%% 800nm pin

[Lia,Locb] = ismember(list(:,1),list(:,6));

[row,col,V] = find(Lia);

for f = 1:length(row)
    condvals_800(f) = list(row(f),3);
end

condvals_800 = condvals_800.';
%% 850nm pin

[Lia,Locb] = ismember(list(:,1),list(:,7));

[row,col,V] = find(Lia);

for f = 1:length(row)
    condvals_850(f) = list(row(f),3);
end

condvals_850 = condvals_850.';
%% 890nm pin

[Lia,Locb] = ismember(list(:,1),list(:,8));

[row,col,V] = find(Lia);

for f = 1:length(row)
    condvals_890(f) = list(row(f),3);
end

condvals_890 = condvals_890.';
