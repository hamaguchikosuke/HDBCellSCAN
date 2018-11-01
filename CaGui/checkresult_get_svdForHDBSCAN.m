% how to improve svd data with small number of movie frames?
opsFile='N:\home\ImagingData\DualLickMice\Ai148VGLUT1Cre_1138\20180725_d3\1_500um\regops_Ai148VGLUT1Cre_1138_20180725_d3_plane1_ch1.mat';
load(opsFile,'ops');
%%
TiffFiles{1}='N:\home\ImagingData\DualLickMice\Ai148VGLUT1Cre_1138\20180725_d3\1_500um\plane1_ch1\x4movie\20180725_d3_1_500um_Ai148VGLUT1Cre_1138_2P_plane1_ch1_x4_001.tif';
TiffFiles{2}='N:\home\ImagingData\DualLickMice\Ai148VGLUT1Cre_1138\20180725_d3\1_500um\plane1_ch1\x4movie\20180725_d3_1_500um_Ai148VGLUT1Cre_1138_2P_plane1_ch1_x4_002.tif'

ops.writeSVDroi = 0; % no need to write svd file
[ops, U, Sv, V, Fs, sdmov] = get_svdForHDBSCAN(ops,TiffFiles);