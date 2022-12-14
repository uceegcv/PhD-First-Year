function e=Wray_tissue_specific_extinction_coefficient_650to1042
% Specific	extinction	coeffs	taken from Wray 1987 
% S. Wray, M. Cope, D. T. Delpy, J. S. Wyatt, and E. O. Reynolds, 
%“Characterization of the near infrared absorption spectra of cytochrome aa3 and haemoglobin for the non-invasive monitoring of cerebral oxygenation.,” 
%Biochim. Biophys. Acta, vol. 933, pp. 184–192, 1988.
% 
% specific extinction coefficients in	OD/M/cm

%__________________________________________________________________________
% Specific Extinction is described using base 10 logarithm units and Specific Absorption is
% described using natural logarithm units. The difference between them is a
% scaling factor of ln10:specific absorption coefficient = specific
% extinction coefficient x 2.3025851.
%__________________________________________________________________________

% ---------------------------------------------------------------------------------------------------------
%Wavelengths(nm) HbO2 (OD/M/cm)	HHb (OD/M/cm) 
%
% ---------------------------------------------------------------------------------------------------------
c=[650 506 3743
    653 479 3644
    657 459 3546
    660 445 3442
    663 435 3326
    667 430 3193
    670 427 3043
    674 426 2879
    677 425 2713
    681 423 2550
    684 420 2392
    687 417 2243
    691 415 2108
    694 415 1990
    698 416 1887
    701 421 1798
    705 427 1720
    708 435 1647
    711 444 1579
    715 454 1513
    718 466 1450
    721 478 1392
    725 490 1343
    728 503 1307
    732 517 1286
    735 532 1286
    738 546 1307
    742 561 1349
    745 576 1412
    748 592 1490 
    752 608 1574
    755 624 1644
    758 641 1678
    762 658 1660
    765 675 1590
    768 692 1485
    772 710 1365
    775 728 1250
    778 745 1149
    782 763 1066
    785 781 999
    788 799 948
    792 817 908
    795 835 877
    798 852 852
    801 869 832
    805 886 816
    808 903 804
    811 920 796
    814 936 789
    818 952 785
    821 968 782
    824 983 779
    828 998 778
    831 1013 778
    834 1028 777
    837 1043 777
    840 1057 777
    844 1070 778
    847 1084 780
    850 1097 781
    853 1110 784
    857 1122 788
    860 1134 792
    863 1146 797
    866 1157 803
    869 1167 810
    872 1177 817
    876 1187 824
    879 1196 833
    882 1204 841
    885 1212 850
    888 1219 858
    891 1226 866
    895 1232 873
    898 1238 880
    901 1243 885
    904 1248 890
    907 1252 892
    910 1255 893
    913 1258 892
    916 1260 889
    920 1261 883
    923 1262 875
    926 1263 865
    929 1263 852
    932 1263 836
    935 1262 819
    938 1259 799
    941 1258 778
    944 1255 755
    947 1251 729
    950 1248 702
    953 1244 674
    956 1239 645
    959 1233 616
    962 1226 587
    965 1220 558
    968 1213 529
    972 1205 501
    975 1197 473
    978 1189 445
    981 1180 419
    984 1170 394
    986 1160 370
    989 1149 347
    992 1138 325
    995 1126 305
    998 1114 285
    1001 1103 267
    1004 1089 250
    1007 1076 234
    1010 1062 219
    1013 1047 206
    1016 1033 193
    1019 1018 182
    1022 1003 171
    1025 988 162
    1028 972 153
    1031 956 145
    1034 940 138
    1036 923 131
    1039 905 125
    1042 889 120];    

%Extrapolate   
e = zeros([3 length([650:1042])]);
e(1,:) = [650:1042]';
e(2,:) = interp1(c(:,1),c(:,2),[650:1042],'linear','extrap');         
e(3,:) = interp1(c(:,1),c(:,3),[650:1042],'linear','extrap');    
e = e';                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          