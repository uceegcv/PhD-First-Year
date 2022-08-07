%pipelineExampe_RJC

DPFs = [6 6 6 6 6];
path_dep = DPF_Lambda_Dependency_740to915;
chanToPlot = 1;

%Convert LUMO data to .nirs
LUMOFileName = ['/Users/georginaleadley/Documents/PhD Data/5_wav_test2_Liam27-5-2022_14-22.LUMO'];
nirs = DOTHUB_LUMO2nirs_5wav(LUMOFileName);

%convert intenstiy to DOD
dod = hmrIntensity2OD(nirs.d);

%convert to concentration
conc = hmrOD2Conc5wav(dod,nirs.SD,DPFs);

%plot concentrations
figure;
plot(nirs.t,squeeze(conc(:,:,chanToPlot)));


