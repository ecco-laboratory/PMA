%% Step 0: Configure Paths and Flag options
set_up_dirs_paths; %configure folders to load data, masks, code, and for writing output
write_out = false; %flag to write results to disk
%% Step 1: Load Training Data
load([datadir filesep 'alldata.mat']); %load processed training data


%% Step 2: clean up csf, include voxels that are either .3 probability in SPM TPM or are included in Pauli subcortical atlas
tpm = fmri_data(which('TPM.nii')); % tissue probability map from SPM
tpm.dat = double(tpm.dat(:,1)>.3); % include GM voxels with liberal threshold
subcort=(load_atlas('CIT168')); %Pauli atlas
subcort= replace_empty(resample_space(subcort,tpm));
tpm=replace_empty(tpm);
tpm.dat=tpm.dat+double(subcort.dat>0);
alldat = apply_mask(alldat,tpm); %mask out voxels 


%% Step 3: create a combined mask
%load individual regions and create mask of all regions
Amygdala=select_atlas_subset(load_atlas('Canlab2018'),{'Amy'});
subcort=merge_atlases(load_atlas('CIT168'),Amygdala);
insula = select_atlas_subset(load_atlas('glasser'),{'Ctx_PoI1_L','Ctx_Ig_L','Ctx_PI_L','Ctx_POI2_L','Ctx_MI_L','Ctx_AI_L','Ctx_AVI_L','Ctx_AAIC_L','Ctx_PoI1_R','Ctx_Ig_R','Ctx_PI_R','Ctx_POI2_R','Ctx_MI_R','Ctx_AI_R','Ctx_AVI_R','Ctx_AAIC_R'});
insula = resample_space(insula,subcort);
MFC=fmri_data(which('AllRegions.nii'));
all_regs=MFC;
all_regs=resample_space(all_regs,subcort);
all_regs.dat(subcort.dat>0)=1;
all_regs.dat(insula.dat>0)=1; %mask with all regions

% create mask with all cortical areas
cortex = MFC;
cortex=resample_space(cortex,subcort);
cortex.dat(insula.dat>0)=1;

%create masked data object using all regions
masked_dat=apply_mask(alldat,all_regs);
wh_removed=(masked_dat.removed_images);

%% Step 4: Perform univariate t-tests to examine voxelwise stats

% model average response
masked_dat.X=ones(size(masked_dat.dat,2),1);
st=regress(masked_dat); %  GLM

% create stats object and write
t=st.t;
t.dat=t.dat(:,end);
t.p=t.p(:,end);
t.sig=t.sig(:,end);

orthviews(threshold(t,.05,'FDR','k',10)); drawnow;snapnow;

t.fullpath = [resultsdir '\mean_response.nii'];

if write_out
    write(t,'overwrite');
end

% model positive vs negative affect 
masked_dat.X=[zscore(study<11)];
st=regress(masked_dat); %pos vs neg

% create stats object and write
t=st.t;
t.dat=t.dat(:,1);
t.p=t.p(:,1);
t.sig=t.sig(:,1);
orthviews(threshold(t,.01,'FDR','k',10)); drawnow;snapnow;
t.fullpath = [resultsdir '\positive_vs_negative_affect.nii'];

if write_out
    write(t,'overwrite')
end

%% Step 5: Specify Model for PLS regression, following Kragel et al. 2018

studyInds=condf2indic(study); %matrix of 0/1 based on study membership
subdomainInds=[[condf2indic(ceil(study(1:224)/2)); zeros(270,5)] [zeros(224,9); condf2indic(ceil((study(225:end)-10)/2))]]; %matrix of 0/1 based on subdomain membership
domInds=[[ones(224,1); zeros(270,1)] [zeros(224,3); condf2indic(ceil((study(225:end)-10)/6))]];%matrix of 0/1 based on domain membership

Y = [double(study<29) domInds subdomainInds studyInds]; % combine into single matrix

figure; imagesc(squareform(pdist(Y)));colorbar;colormap(gray); drawnow;snapnow; %show model
title('Distances between images based on functional ontology')
drawnow;snapnow;
%% Step 6: Run regressions separately for each functional domain

names={'' 'pleasure' 'pain' 'cognitive' 'negative'};

for r=2:5
    masked_dat=apply_mask(alldat,all_regs);
    wh_removed=(masked_dat.removed_images); % find images with missing coverage in an roi
    tv = condf2indic(study(Y(:,r)>0));
    masked_dat.dat = masked_dat.dat(:,Y(:,r)>0);
    masked_dat.X = zscore(tv);
    st=regress(masked_dat); % simple GLM

    t=st.t;
    t.dat=t.dat(:,end);
    t.p=t.p(:,end);
    t.sig=t.sig(:,end);
    orthviews(threshold(t,.05,'FDR','k',10)); drawnow; snapnow;%threshold

    if write_out
        t.fullpath = [resultsdir '\mean_response_' names{r} '.nii'];
        write(t,'overwrite');
    end
end




%% Step 7: Perform RSA using data from all regions, following methods in Kragel et al. 2018, see https://doi.org/10.1038/s41593-017-0051-7
masked_dat=apply_mask(alldat,all_regs);
    wh_removed=(masked_dat.removed_images); % find images with missing coverage in an roi
stats_allreg = rsa_regression(masked_dat,Y(~wh_removed,2:end),study(~wh_removed)); 

%% Step 8: Perform RSA using data from individual ROIS

Amygdala=select_atlas_subset(load_atlas('Canlab2018'),{'Amy'});
masked_dat=apply_mask(alldat,Amygdala);
wh_removed=(masked_dat.removed_images);
stats_Amygdala = rsa_regression(masked_dat,Y(~wh_removed,2:end),study(~wh_removed));

masked_dat=apply_mask(alldat,insula);
wh_removed=(masked_dat.removed_images);
stats_insula = rsa_regression(masked_dat,Y(~wh_removed,2:end),study(~wh_removed));

NAc=select_atlas_subset(load_atlas('CIT168'),{'NAC'});
masked_dat=apply_mask(alldat,NAc);
wh_removed=(masked_dat.removed_images);
stats_NAc = rsa_regression(masked_dat,Y(~wh_removed,2:end),study(~wh_removed));

VeP= select_atlas_subset(load_atlas('CIT168'),{'VeP'});
masked_dat=apply_mask(alldat,VeP);
wh_removed=(masked_dat.removed_images);
stats_VeP = rsa_regression(masked_dat,Y(~wh_removed,2:end),study(~wh_removed));

vmPFC=fmri_data(which('vmPFC.nii'));
masked_dat=apply_mask(alldat,vmPFC);
wh_removed=(masked_dat.removed_images);
stats_vmPFC = rsa_regression(masked_dat,Y(~wh_removed,2:end),study(~wh_removed));


%% Step 9: Plot results from RSA 
names = {'Pleasure' 'Pain' 'Cog' 'Neg Affect' 'Music' 'Food' 'Sexual' 'Money' 'Social'};

C=seaborn_colors(10);
C=C([2 6:10 3:5],:);


figure;
subplot(1,2,1)
distributionPlot(stats_insula.bs_gen_index(:,[2 6:10 3:5]),'showMM',0,'color',C)
set(gca,'XTickLabel',names([1 5:9 2:4]))
title 'Generalizability of Insula Representations'
ylabel 'Generalization Index'
subplot(1,2,2)
imagesc(stats_insula.RDM(~isnan(stats_insula.RDM(:,1)),~isnan(stats_insula.RDM(:,1))))
set(gca,'Visible','off')
colorbar;
colormap(gray)
drawnow;snapnow;


figure; subplot(1,2,1)
distributionPlot(stats_NAc.bs_gen_index(:,[2 6:10 3:5]),'showMM',0,'color',C)
set(gca,'XTickLabel',names([1 5:9 2:4]))
title 'Generalizability of NAc Representations'
ylabel 'Generalization Index'
subplot(1,2,2)
imagesc(stats_NAc.RDM(~isnan(stats_NAc.RDM(:,1)),~isnan(stats_NAc.RDM(:,1))))
set(gca,'Visible','off')
colorbar;
colormap(gray)
drawnow;snapnow;

figure; subplot(1,2,1)
distributionPlot(stats_VeP.bs_gen_index(:,[2 6:10 3:5]),'showMM',0,'color',C)
set(gca,'XTickLabel',names([1 5:9 2:4]))
title 'Generalizability of VeP Representations'
ylabel 'Generalization Index'
subplot(1,2,2)
imagesc(stats_VeP.RDM(~isnan(stats_VeP.RDM(:,1)),~isnan(stats_VeP.RDM(:,1))))
set(gca,'Visible','off')
colorbar;
colormap(gray)
drawnow;snapnow;


figure; subplot(1,2,1)
distributionPlot(stats_vmPFC.bs_gen_index(:,[2 6:10 3:5]),'showMM',0,'color',C)
set(gca,'XTickLabel',names([1 5:9 2:4]))
title 'Generalizability of vmPFC Representations'
ylabel 'Generalization Index'
subplot(1,2,2)
imagesc(stats_vmPFC.RDM(~isnan(stats_vmPFC.RDM(:,1)),~isnan(stats_vmPFC.RDM(:,1))))
set(gca,'Visible','off')
colorbar;
colormap(gray)
drawnow;snapnow;



%% Step 10: Single PLS model using all regions, save weights for generalization tests
masked_dat=apply_mask(alldat,all_regs);
wh_removed=(masked_dat.removed_images);

[xl,yl,xs,ys,bpls]=plsregress(masked_dat.dat',zscore(Y(~wh_removed,:)),16);

save([resultsdir filesep 'bpls.mat'], 'bpls', 'masked_dat')
pos_affect_pattern = masked_dat;
pos_affect_pattern.dat = bpls(2:end,2);
pos_affect_pattern.removed_images =0;
save([resultsdir filesep 'pos_affect_pattern.mat'], 'pos_affect_pattern');
PLS_block_bootstrap; %for inference

%2018 paper and wold2001 (chemometrics)


%% Step 11:  Perform 10-fold cross-validation to estimate performance for all domains
masked_dat=apply_mask(alldat,all_regs);
wh_removed=(masked_dat.removed_images);
kinds = crossvalind('Kfold',study,10);
clear yhat
for k=1:max(kinds)
    train = kinds ~=k;
    test = ~train;
    [xl,yl,xs,ys,bpls]=plsregress(masked_dat.dat(:,train)',zscore(Y(train,:)),16);
    yhat(test,:)=[ones(length(find(test)),1) masked_dat.dat(:,test)']*bpls(:,2:5);
end

figure;
subplot(1,2,1);
[~,inds]=max(yhat,[],2);
cm = confusionmat(double(domInds*(1:4)'),inds);
cm = bsxfun(@rdivide,cm,sum(cm,2));
imagesc(cm);colormap(flipud(gray)); colorbar;
set(gca,'XTickLabel',{'Pleasure' 'Pain' 'Cognitive Control' 'Negative Emotion'})
set(gca,'YTickLabel',{'Pleasure' 'Pain' 'Cognitive Control' 'Negative Emotion'})
set(gca,'XTick',1:4);
set(gca,'YTick',1:4);
axis square
subplot(1,2,2); hold all; C ={[.98 .5 .2],[1 0 0],[0 1 0],[0 0 1]};

for c=1:4
    ROC_domain(c)=roc_plot(yhat(:,c),Y(:,1+c)>0,'color',C{c},'plotmethod','observed');
end
plot([0 1],[0 1],'k-')
set(gca,'FontSize',10)
axis square
drawnow;snapnow;

%% Step 12: Create surface plots for regions of interest
volume_renders
PlotPatternsByROI
quantify_gradients

%% Step 13: Test in 4 independent studies
generalization_tests

%% Step 14: Spatial regressions (and correlations) with gene expression maps
map_onto_gene_expression_maps

%% Step 15: Evaluate expression in data from naloxone challenge
naloxone_test