clc;
clear all
close all

eps=0.00001;
lh=5;
tekrar=5000;
pic1=imread('dataset\lena.tif');
pic= im2gray(pic1);
[m,n]=size(pic);
uilk=imnoise(pic, 'salt & pepper', 0.2);
pic=double(pic);
uilk=double(uilk);
yeni=uilk;
for i=2:m-1
    for j=2:n-1
        c=0; top=0;
if ((uilk(i,j)==0 || uilk(i,j)==255)) 
    if uilk(i-1,j)~=0 && uilk(i-1,j)~=255 c=c+1; top=top+uilk(i-1,j); else; end
 if uilk(i+1,j)~=0  && uilk(i+1,j)~=255 c=c+1; top=top+uilk(i+1,j); else; end
  if uilk(i,j-1)~=0  && uilk(i,j-1)~=255 c=c+1; top=top+uilk(i,j-1); else; end
if uilk(i,j+1)~=0  && uilk(i,j+1)~=255 c=c+1; top=top+uilk(i,j+1); else; end
if c~=0; yeni(i,j)=top/c; else; end
        
end 
    end
end
yeni2=f_TVFF_1(pic, yeni,m,n,5,eps);
% yeni2=f_TVFF_2(pic, yeni,m,n,5,eps);
for i=2:m-1
    for j=2:n-1
if ((uilk(i,j)==0 || uilk(i,j)==255)) 
   yeni2(i,j)=(yeni2(i-1,j)+yeni2(i+1,j)+yeni2(i,j-1)+yeni2(i,j+1))/4;  
end 
    end
end
yeni2=uint8(yeni2);
pic=uint8(pic);
[peaksnr, snr] = psnr(yeni2,pic);
[ssimIndex2, ~] = ssim(yeni2,pic);
score = multissim(yeni2,pic);
fprintf('psnr: %.4f  ssim: %.4f ms-ssim : %.4f\n',peaksnr,ssimIndex2,score);
figure; 
subplot(1,3,1); imshow(pic); title('Orginal');
subplot(1,3,2); imshow(uint8(uilk)); title('Noisy');
subplot(1,3,3); imshow(yeni2); title('Denoised');

