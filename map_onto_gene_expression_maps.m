%% Step 0: load the signature
load([resultsdir filesep 'pos_affect_pattern.mat'])
write_maps=false; %don't rewrite data
%% Step 1: load the gene expression maps
files = dir([datadir filesep 'ReceptorABA' filesep '*.gz']);
P={};
for f=1:length(files)
    P{f}=[files(f).folder filesep files(f).name];
end

nt_maps = fmri_data(P);
nt_maps = preprocess(nt_maps,'smooth',6);


leg = {'DRD1','DRD2','DRD3','OPRD1','OPRK1','OPRM1'};
if write_maps

    for ii=1:length(leg)
        tv = nt_maps;
        tv.dat=tv.dat(:,ii);
        tv.fullpath = [datadir filesep 'neurotransmitter_maps_' leg{ii} '.nii'];
        write(tv,'overwrite');
    end

end

%% load bootstrap distribution of beta coefficients from PLS regression and resize it
load([tempdir filesep 'bs_b.mat'])
resized_nt_maps = nt_maps; %create temporary variable to resize nt maps
resized_nt_maps = apply_mask(resized_nt_maps,pos_affect_pattern);
resized_nt_maps = resample_space(resized_nt_maps,pos_affect_pattern);
resized_nt_maps = replace_empty(resized_nt_maps);
%% for all bootstrap samples, compute correlation and regression coefficients relating gene expression and signature weights
for k=1:size(bs_b,1)

    load([resultsdir filesep 'pos_affect_pattern.mat']); %initialize with proper voxels included

    pos_affect_pattern.dat = squeeze(bs_b(k,2:end,2))'; %grab bootstrap betas from a given iteration
    pos_affect_pattern=replace_empty(pos_affect_pattern); %replace empty voxels to match size with receptor maps


    %spatial regression
    [b(k,:), dev, st]=glmfit(resized_nt_maps.dat(resized_nt_maps.dat(:,1)~=0 & pos_affect_pattern.dat~=0,:),pos_affect_pattern.dat(resized_nt_maps.dat(:,1)~=0 & pos_affect_pattern.dat~=0));

    %simple correlation
    R = corr((b(k,2:end)*(resized_nt_maps.dat(resized_nt_maps.dat(:,1)~=0 & pos_affect_pattern.dat~=0,:))')'+b(k,1),pos_affect_pattern.dat(resized_nt_maps.dat(:,1)~=0 & pos_affect_pattern.dat~=0));

    for i=1:6 %for each map
        rs(k,i)= corr((b(k,2:end)*(resized_nt_maps.dat(resized_nt_maps.dat(:,1)~=0 & pos_affect_pattern.dat~=0,:))')'+b(k,1),(resized_nt_maps.dat(resized_nt_maps.dat(:,i)~=0 & pos_affect_pattern.dat~=0,i)));
    end

    stats = image_similarity_plot(pos_affect_pattern, 'correlation', 'mapset',resized_nt_maps,'notable','nofigure','noplot');
    r_nt_bootstrap(k,:)=stats.r;
end

%compute stats for pairwise differences 
for i=1:6
    for j=1:6
        mean_bs = mean(b(:,1+i)-b(:,1+j));
        std_bs = std(b(:,1+i)-b(:,1+j));
        bs_Z_paired(i,j)=mean_bs./std_bs;
        bs_P_paired(i,j) = 2*normcdf(-1*abs(bs_Z_paired(i,j)),0,1);
    end
end

%compute stats for spatial regression coefficients
mean_bs = mean(b(:,2:end));
std_bs = std(b(:,2:end));
bs_Z=mean_bs./std_bs;
bs_P = 2*normcdf(-1*abs(bs_Z),0,1);


%compute stats for differences in spatial regression coefficients
mean_bs_diff = mean(b(:,end)-mean(b(:,2:end-1),2));
std_bs_diff = std(b(:,end)-mean(b(:,2:end-1),2));
bs_Z_diff=mean_bs_diff./std_bs_diff;
bs_P_diff = 2*normcdf(-1*abs(bs_Z_diff),0,1);

%compute stats for differences in spatial regression coefficients for
%figure panels
mean_bs_diff_pair = mean(b(:,end)-b(:,3));
std_bs_diff_pair = std(b(:,end)-b(:,3));
bs_Z_diff_pair=mean_bs_diff_pair./std_bs_diff_pair;
bs_P_diff_pair = 2*normcdf(-1*abs(bs_Z_diff_pair),0,1);

C = seaborn_colors(10);
figure;
distributionPlot(b(:,2:end),'showMM',0,'color',C([4 3 6 1 5 2],:)); legend(leg);
set(gca,'XTick','')
ylabel 'Beta Coefficient'
drawnow;snapnow; 
%% Create figure and estimate stats using correlation
figure;
distributionPlot(r_nt_bootstrap,'showMM',0,'color',C([4 3 6 1 5 2],:)); legend(leg);
set(gca,'XTick','')
ylabel 'Correlation Coefficient'
drawnow;snapnow; 

%average effects
mean_bs_corr = mean(r_nt_bootstrap);
std_bs_corr = std(r_nt_bootstrap);
bs_Z_corr=mean_bs./std_bs_corr;
bs_P_corr = 2*normcdf(-1*abs(bs_Z_corr),0,1);

%differences
mean_bs_diff_corr = mean(r_nt_bootstrap(:,end)-mean(r_nt_bootstrap(:,2:end-1),2));
std_bs_diff_corr = std(r_nt_bootstrap(:,end)-mean(r_nt_bootstrap(:,2:end-1),2));
bs_Z_diff_corr=mean_bs_diff_corr./std_bs_diff_corr;
bs_P_diff_corr = 2*normcdf(-1*abs(bs_Z_diff_corr),0,1);