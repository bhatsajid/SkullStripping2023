function [mask skull] = test5(fimg,neeta,neetab)
%filling(imfill)


[x,y,z] = size(fimg);
if(z == 1)
fimg(:,:,2)=fimg(:,:,1);
fimg(:,:,3)=fimg(:,:,1);
end

originalimage=fimg;
    
  
   
    smoriginal=fimg;
    
   
   
   global histeqimg;
   
   histeqimg=histeq(originalimage,256);
   
   origimage = fimg;
   
   histeqimage = histeqimg;
   
     hsvimg = rgb2hsv(fimg);  
% 
  hsvimg = im2uint8(rescale(hsvimg));
  
  labimg = rgb2lab(fimg);  
% 
  labimg = im2uint8(rescale(labimg));
   

    
	
	[mask, skull]=binning(fimg, neeta, neetab); %Binning
	
end