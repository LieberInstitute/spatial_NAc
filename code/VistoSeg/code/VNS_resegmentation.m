%fname = 'nuclei.mat';
load(fname)
BW = bwareafilt(mask_dark_blue,[20 2000]); % 'n' will the cluster number associated with nuclei from 2nd run
save('*nuclei_final.mat','BW')
