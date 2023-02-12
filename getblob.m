function [final] = getblob(image)

%checking if input parameter is 2D Image or path to image

% T1= [0.1, 0, 0; 0, 0.1, 0; 0, 0, 1];
% T2 = [10, 0, 0; 0, 10, 0; 0, 0, 1];
% t1 = maketform('affine',T1);
% t2 = maketform('affine',T2);

oriim = image;
h=ones(7,7);
	               
                    
                    
                                 
                    
                    
   %  histlev = [5, 6, 7];
   %   histlev = [256];
   histlev = [256];

%	s=floor(meanstd/2)




%procesing
[x y z]=size(image);  %getting image size



if(z == 1)
image(:,:,2)=image(:,:,1);
image(:,:,3)=image(:,:,1);
end






%TRY CONVERT RGB TO HSV HERE
% 
  hsvimg = rgb2hsv(image);  %rgb2hsv
% 
  hsvimg = im2uint8(rescale(hsvimg));
  
  labimg = rgb2lab(image);  %rgb2lab
% 
  labimg = im2uint8(rescale(labimg));
  
 % image = hsvimg;






image1 = bitand(bitand(image,hsvimg), bitand(image,labimg)); %close skull stripping with levels 2
%image = bitor(bitand(image,hsvimg),bitand(image,labimg));
%image = bitxor(bitor(image,hsvimg), bitor(image,labimg)); %for tumor level 3
%image = bitxor(bitxor(bitand(image,hsvimg), bitand(image,labimg)),bitor(bitand(image,hsvimg), bitand(image,labimg))); %for skull stripping
%image = bitand(bitxor(bitand(image,hsvimg), bitand(image,labimg)),bitor(bitand(image,hsvimg), bitand(image,labimg))); %for skull stripping
%image = bitor(bitxor(bitxor(image,hsvimg), bitxor(image,labimg)),bitxor(bitand(image,hsvimg), bitand(image,labimg))); %for tumor
%image = bitand(bitxor(image,hsvimg), bitxor(image,labimg)); %mrisegmentation
%image = bitand(bitor(image,hsvimg), bitor(image,labimg)); %mrisegmentation
image2 = bitxor(bitand(image,hsvimg), bitand(image,labimg)); %close skull stripping with levels 2

%image1=averagefilt(image1); image2=averagefilt(image2);


image= double(double(image1)+double(image2));

%image= image2-image1;
image=rescale(image);
image=im2uint8(image);
%image=imadjust(image,[0.05 0.3 0.15; 1 1 1],[]);
image=imadjust(image,[0.6 0.1 0.1; 1 1 1],[]);
gr1=rgb2gray(image1);
%roughc=rescale(gr1);
%gr1=im2uint8(roughc);
% gr1 = roughclustering(gr1,2,x,y);
gr2=rgb2gray(image2);
%roughc=rescale(gr2);
%gr2=im2uint8(roughc);
% gr2 = roughclustering(gr2,2,x,y);

grf=double(gr1)+double(gr2);

grf=rescale(grf);
grf=im2uint8(grf);
grf = roughclustering(grf,11,x,y);




    se2 = strel('disk',0);


grf = imerode(grf,se2,'full');


 se1 = strel('disk',0);


grf = imdilate(grf,se1,'full');

[x,y,z] = size(grf);


grf=averagefilt(grf);





roughc=image;

roughc=rescale(roughc);
roughc=im2uint8(roughc);

for i=1:length(histlev)

    roughc(:,:,1) = roughclustering(image(:,:,1),histlev(i),x,y);
    roughc(:,:,2) = roughclustering(image(:,:,2),histlev(i),x,y);
    roughc(:,:,3) = roughclustering(image(:,:,3),histlev(i),x,y);
end  



%     se2 = strel('disk',0);
% 
% 
% roughc = imerode(roughc,se2,'full');
% 
%     se1 = strel('disk',);
% 
% 
% roughc = imdilate(roughc,se1,'full');

[x,y,z] = size(roughc);

image=roughc;
%image=averagefilt(roughc);

 % image=histeq(image,128);    
% 
%     %%figure('WindowStyle','docked');imshow(image);

  
    grf=rgb2gray(image);
    
%    min(min(grf))
%    imhist(histeq(grf))
    
%     %%figure
% %%figure('WindowStyle','docked');imshow(grf,[]);
    
    
      grf = histeq(grf,256);  %worked
   
    

%    
%         bone=imbinarize(grf,0.8);
% %%figure('WindowStyle','docked');imshow(bone,[]);
  
%blob=imadjust(grf,[0.0 1],[0.1,1]);
  
%   blob = imbinarize(blob,0.11);
 
 blob=imadjust(grf,[0.5 1],[0.0,1]);
  
   blob = imbinarize(blob,0.4);
   
    %%figure('WindowStyle','docked');imshow(blob);    title('first blob');

   %*************************************
 %   oriim=imadjust(oriim,[0.6 0.0 0.0; 1 1 1],[]);
   grf=rgb2gray(oriim);
   %%figure('WindowStyle','docked');imshow(blob);
  % imhist(grf);title('orig hist');
%   grf=oriim(:,:,2);
%grf = roughclustering(grf,256,x,y)
 %  grf = histeq(grf,256);
     grf = rescale(log(im2double(grf)+0.001));
     %%figure('WindowStyle','docked')
  %    imhist(grf);title('adjusted hist');
 %   blob=imadjust(grf,[0.5 1],[0.0,1]);
 % grf=adapthisteq(grf,'Distribution','uniform', 'NBins',256,'ClipLimit',0.0);    %was 256 & 80
   %%figure('WindowStyle','docked');imshow(grf);       title('adjusted gray');
      
%      haha
  blob=grf;
blob=imadjust(blob,[]);
 histy = imhist(blob);
 
  %%figure('WindowStyle','docked');
 % imhist(blob);
 
% histy(find(histy==max(histy)))=0;
  
 %histy(find(histy<max(histy)*0.001))=0;
% histy(find(histy==max(histy)))=0;

for i=2:130
   histy(i)=0; 
    
end


for i=170:256 %161
   histy(i)=0;     
end

histy(256)=histy(1);
 
 % bar(histy);

 maxi = islocalmax(histy);
 
% histy(find(histy<max(histy)*0.1))=0;
%  %%figure('WindowStyle','docked');
% bar(histy);
 
 % indy = find(histy,10,'first')/256;
  [locmin prom] = islocalmin(histy(find(histy)));
  
 % indo = find(prom==max(prom),1,'first') %earlier
 
  indo = find(prom,1,'first');
  if(isempty(indo))
      indo=round(length(prom)/2);
  end
  
  indy = find(histy,indo,'first');
  indy(indo);
  
 %blob=imadjust(blob,[indy(2) 1],[0 1]);
%  hold on
  
 %  plot(128,histy(128),'g*:')
  
 % plot(indy(indo),histy(indy(indo)),'r*:')
%legend('hist','mid','first local min');
  
%  hold off
  
  
  
  indexy = (indy(indo)/256);
  if indexy < 0.5  %was 0.7
     indexy = 0.5;
  end
 % indexy=0.625
   blob = imbinarize(blob,indexy);
   
 %   %%figure('WindowStyle','docked');imshow(blob);

   
   se2 = strel('disk', 0);
   blob=imdilate(blob, se2);
   blob=imerode(blob, se2);
   
    se2 = strel('disk', 0);
    blob=imerode(blob, se2);
  % blob=imdilate(blob, se2);
   se2 = strel('disk', 0);
    blob=imdilate(blob, se2);
    
 %   blob = medfilt2(blob);
    
  % blob = bwmorph(blob, 'spur', Inf);
  
     blob = bwmorph(blob, 'spur', 4);
   blob = bwmorph(blob, 'clean');
 
  % blob = bwmorph(blob, 'thicken');

   
    %%figure('WindowStyle','docked');imshow(blob);
%    
 % haha

   %**********************************
   


final = im2uint8(blob);


end