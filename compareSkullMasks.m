path = 'C:\Users\Raptor\Documents\MATLAB\MRI Experiment\NFBSexp\';
caseFile= 'A00033747';


groundTruth = [path,caseFile,'\','sub-',caseFile,'_ses-NFB3_T1w_brainmask.nii.gz']
%testFile = [path, caseFile,'\','FSL_BET-T1_mask.nii.gz']
%testFile = [path, caseFile,'\','BEAST_T1.nii.gz']
%testFile = [path, caseFile,'\','BrainSuite_T1.nii.gz']
%testFile = [path, caseFile,'\','HD_BET_T1_mask.nii.gz']
%testFile = [path, caseFile,'\','Robex_T1.nii']
%testFile = [path, caseFile,'\','SKST_T1.nii.gz']
%testFile = [path, caseFile,'\','SwissSS_T1.nii.gz']
testFile = [path, caseFile,'\','SKST_T1Adjusted.nii.gz']
%testFile = [path, caseFile,'\','sub-',caseFile,'_ses-NFB3_T1w_MBSS.nii.nii']



mask_vol= niftiread(groundTruth);
test_vol=niftiread(testFile);  

[x,y,z1,t]=size(mask_vol)

[x,y,z2,t]=size(test_vol)
     
    %  test_vol=permute(test_vol,[2,3,1]);  %for Robex and SwissSS and FSLBET and HDBET
    %   test_vol=permute(test_vol,[2,3,1]);  %for BEAST, SKST
   % mask_vol=permute(mask_vol,[3,1,2]);  %for BEAST, SKST
    
[x,y,z1,t]=size(mask_vol)

[x,y,z2,t]=size(test_vol)

Results=zeros(1,8);
count =0;
for slice = 1:z1
%slice = 50;
%round(z1/2)+1
maskImg =  imbinarize(im2uint8(rescale(mask_vol(:,:,slice))),0);

testImg =  imrotate(im2uint8(rescale(test_vol(:,:,slice))),0); %for BEAST


%round(z2/2)+1
%testImg =  imrotate(im2uint8(rescale(test_vol(:,:,slice))),0);  %for MBSS, FSL, HDBET, Robex
testImg = imbinarize(testImg, 0.0);  % for Robex use 0 for MBSS use 0.8

%haha
[Accuracy, Sensitivity, FMeasure, Precision, MCC, Dice, Jaccard, Specificity]=EvaluateImageSegmentationScores(maskImg,testImg);
if(~isnan(Precision)&&~isnan(FMeasure)&&~isnan(Sensitivity)&&~isnan(MCC)&&~isnan(Dice))
Results(1) = Results(1)+Accuracy;
Results(2) = Results(2)+Sensitivity;
Results(3) = Results(3)+FMeasure;
Results(4) = Results(4)+Precision;
Results(5) = Results(5)+MCC;
Results(6) = Results(6)+Dice;
Results(7) = Results(7)+Jaccard;
Results(8) = Results(8)+Specificity;
count=count+1;
end
end
disp('Accuracy, Sensitivity, FMeasure, Precision, MCC, Dice, Jaccard, Specificity')
Results./count
figure; subplot(1,2,1);
size(maskImg)
imshow(mask_vol(:,:,110),[]); title('Masks orient');

 subplot(1,2,2);
 imshow(imrotate(test_vol(:,:,110),0),[]); 

%no roation for FSLBET