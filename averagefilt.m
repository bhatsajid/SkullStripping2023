function img = averagefilt(img)
%Performs guassian smoothing (uses averaging 3x3 filter)

 h = fspecial('average',3);
 img=imfilter(img,h);


end