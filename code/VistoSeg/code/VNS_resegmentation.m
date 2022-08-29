function VNS(fname,N)
%fname = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Lieber_Institute_OTS-20-7748_rush_posterior_A1.tif
%N = 5;

tic
disp('Importing capture area')
Img1 = imread(fname); %import image
toc

tic
disp('Performing rgb to Lab color space conversion')
lab_he = rgb2lab(Img1); % convert from rgb color space to Lab color space
toc
ab = lab_he(:,:,2:3); % extract a*b color space from Lab
ab = im2single(ab);
tic
disp('Applying Kmeans')
pixel_labels = imsegkmeans(ab,N,'NumAttempts',3); % apply Kmeans
toc

tic
disp('saving outputs')
parfor i = 1:N
mask{i} = pixel_labels==i;
cluster{i} = he .* uint8(mask{i});
imwrite(cluster{i},[fname(1:end-4),'_cluster',num2str(i),'.tif'])
end

save([fname(1:end-4),'_mask.mat'],'mask','-v7.3')
save([fname(1:end-4),'_cluster.mat'],'cluster','-v7.3')
toc 

load('*mask.mat')
BW = bwarefilt(mask{n},[20 2000]); % 'n' will the cluster number associated with nuclei from 2nd run
save('*nuclei_final.mat','BW')