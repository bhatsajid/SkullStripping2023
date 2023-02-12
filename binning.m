function [mask, skull] = binning(image, neeta, neetab)


image(:,:,2)=image(:,:,1);
image(:,:,3)=image(:,:,1);

prefillblob = getblob(image);
image = padarray(image,[55 55], 'both');  
[x,y,z] = size(image);	


prefillblob = padarray(prefillblob,[55 55], 'both');

blob = imbinarize(prefillblob);






oriim = image;
grfo = rgb2gray(oriim);
[Xa, Ya, ch] = size(grfo);
xer = round(Xa.*0.05);  

yer = round(Ya.*0.05);

if(rem(xer,2)==0)
    xer = xer + 1;
end


if(rem(yer,2)==0)
    yer = yer + 1;
end

h=ones(xer,yer);
stdfiltimg=stdfilt(grfo,h); 
stdvect=stdfiltimg(:);
meanstd=mean(stdvect);

meanstd=round(meanstd/double((max(max(grfo)))),2);


 
  CC = bwconncomp(blob);
 
 numPixels = cellfun(@numel,CC.PixelIdxList);
 [biggest,idx] = max(numPixels);
 blobby=zeros(x,y);
 if(~isempty(idx))  
    blobby(CC.PixelIdxList{idx}) = 1;
   
 end
 
 

   
  
   blobby = bwmorph(blobby,'spur', Inf);

      blobby = bwmorph(blobby,'close');

      blobby=imfill(blobby,'holes');
      
       blobby=bwmorph(blobby,'thin',4);
     blobby = bwmorph(blobby,'spur', Inf);
       
       blobby = bwmorph(blobby,'thicken',4);
       
 
 blobby =   bwmorph(blobby,'bothat');
  se0 = strel('disk',2);
 blobby = imdilate(blobby,se0);
 
 
 
 
 prefillblob = prefillblob - im2uint8(blobby);
 

 
 
 
 
 
 se3 = strel('disk',1);
 blob = imdilate(blob,se3);
 blob = bwmorph(blob,'close',3);
 

 [x,y,z] = size(blob);

se3 = strel('disk',round(min(x,y)*0.01));
blob = imdilate(blob,se3,'full');


 blob=imfill(blob,8,'holes');


se3 = strel('disk',round(min(x,y)*0.01));
blob = imerode(blob,se3,'full');



 se3 = strel('disk',round(min(x,y)*0.03));
blob = imerode(blob,se3,'full');
blob = imdilate(blob,se3,'full'); 
 
 
[x,y,z] = size(blob);
 CC = bwconncomp(blob);
 
 numPixels = cellfun(@numel,CC.PixelIdxList);
 [biggest,idx] = max(numPixels);
 
 blob=zeros(x,y,1);
 blob(CC.PixelIdxList{idx}) = 1;

  

 blob = im2uint8(blob);

   




[x y z]=size(image); 



if(z == 1)
image(:,:,2)=image(:,:,1);
image(:,:,3)=image(:,:,1);
end



[x1,y1,z1] = size(blob);
 
xlow=round((x1-x)/2);
xhigh=round(x1-(x1-x)/2)-1;
ylow=round((y1-y)/2);
yhigh=round(y1-(y1-y)/2)-1;


blobnew(1:x,1:y) = blob(xlow:xhigh,ylow:yhigh);



outline=imbinarize(blobnew,'adaptive','sensitivity', 0.3);

thinoutline=imbinarize(blobnew,'adaptive','sensitivity', 0.24);

thickthresh = round(min(size(blobnew))*0.05);

se11 = strel('disk',thickthresh); 

thinoutline = imdilate(thinoutline,se11);
thinoutline = imerode(thinoutline,se11);



 

thinthinoutline=imbinarize(im2uint8(bwmorph(blobnew,'thicken',1)),'adaptive','sensitivity', 0.0); %was 0.1  
 



 prefillblob = and(blobnew,prefillblob);
 
 
 





CC = bwconncomp(outline);

numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);

size(outline);
outline=zeros(x,y);

if(~isempty(idx))
    outline(CC.PixelIdxList{idx}) = 1;
end



[blobx, bloby,ch] = size(blob);
[imx, imy, ch2] = size(oriim);

difx = round((blobx-imx)/2);
dify = round((bloby-imy)/2);

crop_blob = blob(difx:imx+difx-1, dify:imy+dify-1);

oriim = bitand(crop_blob, oriim);


image=oriim;
  
    
    grf=rgb2gray(image);
    
    
hister=imhist(grf);
minovalu = find(hister,1,'first');
hister(minovalu) = 0;
newovalu = find(hister,1,'first');
zerrs = find(grf==(minovalu-1));

if(~isempty(newovalu))
grf(zerrs)=newovalu;
end
   
    
origrf = grf; 



grf=rescale(log(im2double(grf)+meanstd/5)); 
 
grf=adapthisteq(grf,'NumTiles',[round(imx/70) round(imy/70)],'Distribution','uniform', 'NBins',256,'ClipLimit',0.004);    
grf=imadjust(grf,[],[]); 
 
     
     grfchange=grf;

[Xa, Ya] = size(grf);
xer = round(Xa.*0.5); 

yer = round(Ya.*0.5);

if(rem(xer,2)==0)
    xer = xer + 1;
end


if(rem(yer,2)==0)
    yer = yer + 1;
end

h=ones(xer,yer);
stdfiltimg=stdfilt(grf,h);
stdvect=stdfiltimg(:);
meanstd=mean(stdvect);

meanstd=round(meanstd/double((max(max(grf)))),2);






trygrf = grf;



histo=imhist(trygrf);


firstidx =  find(histo,1,'first');
secndidx =  find(histo,1,'last');

 
histo(firstidx) = 0;
histo(secndidx) = 0;

firstidx =  find(histo,1,'first');
secndidx =  find(histo,1,'last');
 

 
 intensities = find(histo>0 & histo<round(0.02*max(histo))); 
 
for i = 1:length(intensities)
    
grfidx = find(trygrf==intensities(i)-1);
trygrf(grfidx) = round((firstidx+secndidx)/2);  

end


mingrf = find(histo,1,'first');

maxgrf = find(histo,1,'last');




medgrf = round(round(double(((maxgrf+mingrf)/2)))/255,1); 
if(isempty(medgrf))
    disp('empty medgrf')
    medgrf = 0.6;
    
end

grf_tes=imbinarize(im2double(trygrf),medgrf);

test_seeds=bwmorph(blobnew, 'remove').*grf_tes;  

test_seeds = bwmorph(test_seeds, 'shrink',Inf);





test_outline = imreconstruct(test_seeds, grf_tes);
thickthresh = round(min(size(blobnew))*0.05);

se11 = strel('disk',thickthresh); 

test_outline = imdilate(test_outline,se11);
test_outline = imerode(test_outline,se11);
testflag = find(test_outline);


base_outline = test_outline;

test_outline = imfill(test_outline,'holes');

thickthresh = round(min(size(blobnew))*0.1);
se11 = strel('disk',thickthresh);






if(length(testflag)>0)
   medgrf=medgrf * neeta;%0.87 ;    1.5 for T1, 0.015 for GADO 0.56 else for flair    -0.15 for 3DFLAIR  stenc +0.15
else
    medgrf = 0.75 ; %0.56
end




stenc =imbinarize(im2double(trygrf),medgrf+0.1); 




grf=imbinarize(im2double(trygrf),medgrf);




centerblob = imfill(thinoutline,'holes');
centerblob = xor(centerblob,thinoutline);
GG = bwconncomp(centerblob);


numPixels = cellfun(@numel,GG.PixelIdxList);
[biggest,idx] = max(numPixels);

largestcc = zeros(x,y);
centerline=zeros(x,y);

if(~isempty(idx))
largestcc(GG.PixelIdxList{idx}) = 1;

blobcenter = regionprops(largestcc,'Centroid','BoundingBox');
cent=ceil(blobcenter.Centroid);
htwd=ceil(blobcenter.BoundingBox);
else
    mask = zeros(x,y);
    skull = zeros(x,y);
   return; 
end
centerline(round(cent(2)-((htwd(4))*.7)/2):round(cent(2)+((htwd(4))*.7)/2),cent(1)) = 1; 
centerline(cent(2),round(cent(1)-((htwd(3))*.7)/2):round(cent(1)+((htwd(3))*.7)/2)) = 1; 



se6 = strel('disk',4);
centerline = imdilate(centerline,se6);


centerline=insertShape(centerline,'FilledCircle',[cent(1) cent(2) abs(min(htwd(3), htwd(4))/2 * 0.4)]);
centerline=imbinarize(rgb2gray(centerline),0);


fixed = stenc; 
   
   se7 = strel('disk',0); 
  fixed = imdilate(fixed, se7);
  
 



tesimg(:,:,1)=im2uint8(fixed);
tesimg(:,:,1)=im2uint8(centerline);




tesimg(:,:,3) = im2uint8(fixed);
 

seeds=centerline.*fixed;
seeds = bwmorph(seeds, 'shrink',Inf);





removal1 = imreconstruct(seeds, fixed);

 se7 = strel('disk',1);  
 removal1 = imdilate(removal1, se7);

 removal1 = imerode(removal1, se7);






 se7 = strel('disk',1); 
removal1 = imdilate(removal1, se7);

removed = grf - removal1;
idx = find(removed<0);
removed(idx)=0;



tesimg(:,:,1)=im2uint8(removed);
tesimg(:,:,1)=im2uint8(centerline);

removed = bwmorph(removed, 'close',0); 


seeds=centerline.*removed;
seeds = bwmorph(seeds, 'shrink',Inf);


removal2 = imreconstruct(seeds, logical(removed));
 se7 = strel('disk',0); 
removal2 = imdilate(removal2, se7);


removed = removed - removal2;
 

base_filled = imfill(base_outline,'holes');
test_img = bitxor(base_outline,base_filled);

flagtest= find(test_img);

idx = find(removed<0);
removed(idx)=0;





removal = or(removal1,removal2);


removed = bwmorph(removed, 'clean',0);





seeds=thinthinoutline.*removed;    
seeds = bwmorph(seeds, 'shrink',Inf);
less = imreconstruct(seeds, removed);



se7 = strel('disk',2);  
fixed = imerode(fixed, se7);
fixed = imdilate(fixed, se7);




seeds=outline.*fixed;  
seeds = bwmorph(seeds, 'shrink',Inf);
less2 = imreconstruct(seeds, fixed);    

less =or(less,less2);
removed=or(less,removed);





 se7 = strel('disk',1);  
removal = imdilate(removal, se7);

seeds=removed.*removal;
seeds = bwmorph(seeds, 'shrink',Inf);


removal = imreconstruct(seeds, logical(removed));


removal=removal-less;
idx = find(removal<0);


se7 = strel('disk',1); 
removed = imerode(removed, se7);
removed = imdilate(removed, se7);
se7 = strel('disk',1);  
removed = imdilate(removed, se7);
removed = imerode(removed, se7);

seeds=outline.*removed;
seeds = bwmorph(seeds, 'shrink',Inf);
finalstensil = imreconstruct(seeds, removed);






 se7 = strel('disk',0); 
finalstensil = imerode(finalstensil, se7);



se11 = strel('disk',0);  
finalstensil = imdilate(finalstensil,se11);





seeds=thinthinoutline.*prefillblob;  



seeds = bwmorph(seeds, 'shrink',Inf);

seedo=seeds;






finalstensilextra = imreconstruct(seeds, prefillblob);

diffstenc = imfill(xor(prefillblob, finalstensilextra),'holes');
diffstenc = bwmorph(diffstenc,'clean');


thickthresh = round(min(size(prefillblob))*0.05);

se11 = strel('disk',thickthresh); 
diffstenc = imerode(diffstenc,se11);
diffstenc = imdilate(diffstenc,se11);

testflag = find(diffstenc);



 idx = find(removed<0);
 removed(idx)=0;
 removed = logical(removed);

if(false)
   finalstensil=finalstensilextra;
   
   
   CH = bwconvhull(finalstensil);


thinthinoutline=imbinarize(im2uint8(CH),'adaptive','sensitivity', 0.0); 
thinthinoutline=bwmorph(thinthinoutline,'thin',Inf);
finalstensil = bwmorph(finalstensil,'thicken',1);
   
else
 
 
 
 
 
 
 
   
seeds=bwmorph(thinoutline,'thicken',3).*removed;
seeds = bwmorph(seeds, 'shrink',Inf);




se10 = strel('disk',round(min(x,y)*0.015));
closedstenc = imdilate(removed,se10);
closedstenc = imfill(closedstenc,'holes');
se10 = strel('disk',round(min(x,y)*0.016)); 
closedstenc = imerode(closedstenc,se10);

closedstenc = bwmorph(closedstenc,'remove');

 



seeds=or(seeds,closedstenc);



se11 = strel('disk',thickthresh);
finalstensil = imdilate(finalstensil,se11);
finalstensil = imerode(finalstensil,se11);


seedsdisp=thinoutline.*finalstensil;    

finsil = finalstensil;
outlinefinsil2 = bwmorph(finsil,'remove');

finsil = bwmorph(finsil,'bridge',0);
finsil = bwmorph(finsil,'thicken',4); %was 4
finsil = bwmorph(finsil,'bridge',1);
finsil = bwmorph(finsil,'fill');

outlinefinsil1 = bwmorph(finsil,'remove');





refine_img = origrf.*uint8(finsil);   



tester = refine_img;
refine_img = rescale(log(im2double(refine_img)+meanstd));  
refine_img = imadjust(refine_img);

hister = imhist(refine_img);

hister(find(hister,1,'first'))=0;
hister(find(hister,1,'last'))=0;


firstthresho = find(hister==max(hister),1,'first');

secondthresho = find(hister==max(hister),1,'first');




thresho = round(round((find(hister==max(hister),1,'first')/256),2),1); 



refine_img = imbinarize(refine_img,thresho*1.8);   


refine_img = bwmorph(refine_img,'thicken',0);
refine_img = bwmorph(refine_img,'bridge',0);
refine_img = bwmorph(refine_img,'fill');
refine_img = bwmorph(refine_img,'thin',Inf);

refine_img = bwmorph(refine_img,'spur',0);
refine_img = bwmorph(refine_img,'clean');


finalstensil = or(finalstensil,refine_img);

finalstensil = bwmorph(finalstensil,'clean');
finalstensil = bwmorph(finalstensil,'thicken',2);




se10 = strel('disk',round(min(x,y)*0.020)); 
closedstenc = imdilate(finalstensil,se10);
closedstenc = imfill(closedstenc,'holes');
se10 = strel('disk',round(min(x,y)*0.021)); 
closedstenc = imerode(closedstenc,se10);

closedstenc = bwmorph(closedstenc,'remove');





teststenc=finalstensil;
CH = bwconvhull(teststenc);

thinthinoutline=imbinarize(im2uint8(blobnew),'adaptive','sensitivity', 0.0); 

thinthinoutline= bwmorph(thinthinoutline,'thin',Inf);

thinthinoutline2=imbinarize(im2uint8(CH),'adaptive','sensitivity', 0.0);
thinthinoutline2= bwmorph(thinthinoutline2,'thin',Inf);

thinthinoutline= or(or(thinthinoutline,thinthinoutline2),closedstenc);


end


finalmask = finalstensil;



finalmask = or(seedo,finalmask);










 se11 = strel('disk',0); 
 finalmask = imdilate(finalmask,se11);
  finalmask = imerode(finalmask,se11);

  

 
seeds=thinthinoutline.*finalmask;    
seeds = bwmorph(seeds, 'shrink',Inf);
finalmask = imreconstruct(seeds, finalmask);
 
se11 = strel('disk',thickthresh);
finalmask = imdilate(finalmask,se11);
finalmask = imerode(finalmask,se11);




finalmask = or(finalmask,thinthinoutline);

 


se11 = strel('disk',3);
dilatedmask = imdilate(finalmask,se11);







se11 = strel('disk',3);
erodedmask = imerode(dilatedmask,se11);



testmask(:,:,1)=im2uint8(finalmask);
testmask(:,:,2)=im2uint8(grfchange);
testmask(:,:,3)=im2uint8(thinoutline);





testmask(:,:,1)=im2uint8(centerline);
testmask(:,:,2)=im2uint8(thinoutline);
testmask(:,:,3)=im2uint8(grf);

 

 

final(:,:,1)=im2uint8(finalmask);
final(:,:,2)=im2uint8(imcomplement(imfill(prefillblob,'holes')));
final(:,:,3)=im2uint8(grfchange);

[msx,msy]=size(grfchange);
se11 = strel('disk',round(min(msx,msy)*0.0025));
negation = imerode(final(:,:,1),se11);



negation = or(negation, final(:,:,2));
se11 = strel('disk',round(min(msx,msy)*0.02));
negation = imdilate(negation,se11);
negation = imerode(negation,se11);



final(:,:,1)=im2uint8(finalmask);
final(:,:,2)=negation;
final(:,:,3)=im2uint8(grfchange);



brain = bitand(imcomplement(im2uint8(negation)), im2uint8(grfchange));

[Xa, Ya] = size(grf);
xer = round(Xa.*0.2);  

yer = round(Ya.*0.2);

if(rem(xer,2)==0)
    xer = xer + 1;
end


if(rem(yer,2)==0)
    yer = yer + 1;
end


brain=adapthisteq(brain,'NumTiles',[round(x/15) round(y/15)],'Distribution','rayleigh', 'NBins',18,'ClipLimit',0.00,'Alpha', 0.9); 
brain=imadjust(brain,[0.0,0.99],[]);


idxs=find(brain);



std_brain = std(double(brain(idxs)))/double((max(max(brain)))) *(neeta+neetab);  





%======================find a thresh for brain=========================
histo=imhist(brain);


firstidx =  find(histo,1,'first');
secndidx =  find(histo,1,'last');
 
histo(firstidx) = 0;


firstidx =  find(histo,1,'first');
secndidx =  find(histo,1,'last');
 

 intensities = find(histo>0 & histo<round(0.05*max(histo))); %was 0.15
 
for i = 1:length(intensities)
   
grfidx = find(trygrf==intensities(i)-1);
trygrf(grfidx) = round((firstidx+secndidx)/2);  
histo(intensities(i))=0;
end

 

mingrf = find(histo,1,'first');

maxgrf = find(histo,1,'last');



brainmedgrf = round(round(double(((maxgrf-mingrf)/2)+mingrf))/255,2)*(1/std_brain); 
if(isempty(brainmedgrf))
   brainmedgrf = 0.5; 
end


drain= imbinarize(brain,0.79);
brain= imbinarize(brain,brainmedgrf);
brain=or(brain,drain);

grain = bwconncomp(brain);


numPixels = cellfun(@numel,grain.PixelIdxList);
[biggest,idx] = max(numPixels);

largestcc = zeros(x,y);

if(~isempty(find(brain)))
if(~isempty(idx))
    largestcc(grain.PixelIdxList{idx}) = 1;
end

    regprop = regionprops(largestcc,'BoundingBox');
    brain_dim=ceil(regprop.BoundingBox);
  
end


se11 = strel('disk',round(min(msx,msy)*0.005)); 

brain = imerode(brain,se11);
brain = imdilate(brain,se11);




brain=imfill(brain,'holes');

CC = bwconncomp(brain);
 numPixels = cellfun(@numel,CC.PixelIdxList);
 [biggest] = max(numPixels);
 if(isempty(biggest))
    biggest = 0; 
 end
 

se11 = strel('disk',round(min(msx,msy)*0.01)); 
brain = imerode(brain,se11);
brain = imdilate(brain,se11);


se11 = strel('disk',round(min(msx,msy)*0.01)); 

brain = imdilate(brain,se11);
brain = imerode(brain,se11);





final(:,:,1)=im2uint8(finalmask);
final(:,:,2)=im2uint8(bwmorph(logical((im2uint8(brain))),'remove'));
final(:,:,3)=im2uint8(grfchange);

mask = brain;

skull = finalmask; 



end