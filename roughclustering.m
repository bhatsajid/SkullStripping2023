function roughclusters = roughclustering(fimg,clusters,x,y)

	%creating bins
	bins=ceil(255/clusters);
   
	

	for i=1:x
		for j=1:y
			for k=1:bins
				if fimg(i,j) < k*bins
					fimg(i,j)=k*bins;
					break;
				end
			end
		end
	end
	
	roughclusters=fimg;
    
%    roughclusters=histeq(fimg, clusters);

end