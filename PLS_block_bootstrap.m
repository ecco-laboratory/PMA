parfor ib = 1:3000
bs_inds=zeros(size(masked_dat.dat,2),1); %initialize

for i=1:max(study) %for each study
    
    study_inds=find(study==i); %find images that belong to this study
    bs_inds(study_inds)=study_inds(randi([1,length(study_inds)],1,length(study_inds))); %randomly resample with replacement

end


[~,~,~,~,bs_b(ib,:,:)]=plsregress(masked_dat.dat(:,bs_inds)',zscore(Y(~wh_removed,:)),16);
end


mean_bs = mean(bs_b);
std_bs = std(bs_b);
bs_Z=mean_bs./std_bs;
bs_Z=reshape(bs_Z,[],47);
bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
