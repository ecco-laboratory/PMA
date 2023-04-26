make_plot_do_stats = true;

load([resultsdir filesep 'pos_affect_pattern'])

subcort=select_atlas_subset(load_atlas('CIT168'),{'NAC','VeP'});
insula = select_atlas_subset(load_atlas('glasser'),{'Ctx_PoI1_L','Ctx_Ig_L','Ctx_PI_L','Ctx_POI2_L','Ctx_MI_L','Ctx_AI_L','Ctx_AVI_L','Ctx_AAIC_L','Ctx_PoI1_R','Ctx_Ig_R','Ctx_PI_R','Ctx_POI2_R','Ctx_MI_R','Ctx_AI_R','Ctx_AVI_R','Ctx_AAIC_R'});
insula = resample_space(insula,subcort);
vmpfc=fmri_data(which('vmpfc.nii'));
all_regs=vmpfc;
all_regs=resample_space(all_regs,subcort);
all_regs.dat(subcort.dat>0)=1;
all_regs.dat(insula.dat>0)=1;

pos_affect_pattern_hotspots = apply_mask(pos_affect_pattern,all_regs);

brs = fmri_data('C:\Users\phili\Downloads\unthresholded_Z_map_1.nii.gz');

clear P
files = dir('C:\Users\phili\ds000219-download\**\con_0001.nii');
for f=1:length(files)
    P{f}=[files(f).folder filesep files(f).name];
end

pleasant_taste = preprocess(fmri_data(P),'smooth',6);

files = dir('C:\Users\phili\ds000219-download\**\con_0002.nii');
for f=1:length(files)
    P{f}=[files(f).folder filesep files(f).name];
end

not_pleasant_taste = preprocess(fmri_data(P),'smooth',6);

pexp_pleasant_taste = apply_mask(pleasant_taste,pos_affect_pattern_hotspots,'pattern_expression','cosine_similarity');
pexp_not_pleasant_taste = apply_mask(not_pleasant_taste,pos_affect_pattern_hotspots,'pattern_expression','cosine_similarity');

EBDM_Cue2_hrew = fmri_data([datadir filesep 'EBDM_Cue2_low_effort_high_reward.nii']);
EBDM_Cue2_lrew = fmri_data([datadir filesep 'EBDM_Cue2_low_effort_low_reward.nii']);
EBDM_Cue1_hef = fmri_data([datadir filesep 'EBDM_cue1_high_effort.nii']);
EBDM_Cue1_lef = fmri_data([datadir filesep 'EBDM_cue1_low_effort.nii']);

pexp_ebdm_high_rew = apply_mask(EBDM_Cue2_hrew,pos_affect_pattern_hotspots,'pattern_expression','cosine_similarity');
pexp_ebdm_low_rew = apply_mask(EBDM_Cue2_lrew,pos_affect_pattern_hotspots,'pattern_expression','cosine_similarity');

pexp_ebdm_high_eff = apply_mask(EBDM_Cue1_hef,pos_affect_pattern_hotspots,'pattern_expression','cosine_similarity');
pexp_ebdm_low_eff = apply_mask(EBDM_Cue1_lef,pos_affect_pattern_hotspots,'pattern_expression','cosine_similarity');

load('C:\Users\phili\OneDrive - Emory University\Documents\DeepRL_Valence\SavedCode\outcomes.mat')
pexp_outcome_loss_incorrect = apply_mask(outcome_loss_incorrect,pos_affect_pattern_hotspots,'pattern_expression','cosine_similarity');
pexp_outcome_gain_correct = apply_mask(outcome_gain_correct,pos_affect_pattern_hotspots,'pattern_expression','cosine_similarity');

files = dir([datadir filesep 'Humor' filesep 'ISR2' filesep 'Caption_Reward-Control' filesep 'HC' filesep '*.img']);
HC_files = strrep({files(:).name},'con_0006','con_0002');
cd ([datadir filesep 'Humor' filesep 'ISR2' filesep 'Caption_Reward' ]);
dat_HC_caption = fmri_data(HC_files);
cd ([datadir filesep 'Humor' filesep 'ISR2' filesep 'Caption_Reward-Control' filesep 'HC' ]);
dat_HC_caption_vs_control = fmri_data({files(:).name});

dat_control = dat_HC_caption_vs_control;
dat_control.dat=dat_control.dat*-1;
dat_control.dat = dat_control.dat + dat_HC_caption.dat;

pexp_HC_caption = apply_mask(dat_HC_caption,pos_affect_pattern_hotspots,'pattern_expression','cosine_similarity');
pexp_HC_control = apply_mask(dat_control,pos_affect_pattern_hotspots,'pattern_expression','cosine_similarity');

load([datadir filesep 'groupvars.mat']);clear ROC_hotspots;

condition = categorical([ones(length(pexp_ebdm_high_rew),1);-1*ones(length(pexp_ebdm_high_rew),1);ones(length(pexp_outcome_gain_correct(group==2)),1);-1*ones(length(pexp_outcome_gain_correct(group==2)),1); ones(length(pexp_not_pleasant_taste),1); -1*ones(length(pexp_not_pleasant_taste),1); ones(length(pexp_HC_caption),1); -1*ones(length(pexp_HC_caption),1)]);
isPleasure = categorical([-1*ones(length(pexp_ebdm_high_rew),1);-1*ones(length(pexp_ebdm_high_rew),1);-1*ones(length(pexp_outcome_gain_correct(group==2)),1);-1*ones(length(pexp_outcome_gain_correct(group==2)),1); ones(length(pexp_not_pleasant_taste),1); ones(length(pexp_not_pleasant_taste),1); ones(length(pexp_HC_caption),1); ones(length(pexp_HC_caption),1)]);
studyNumber = [ones(2*length(pexp_ebdm_high_rew),1);2*ones(2*length(pexp_outcome_gain_correct(group==2)),1); 3*ones(2*length(pexp_not_pleasant_taste),1); 4*ones(2*length(pexp_HC_caption),1)];
pexp = [pexp_ebdm_high_rew;pexp_ebdm_low_rew; pexp_outcome_gain_correct(group==2);pexp_outcome_loss_incorrect(group==2);pexp_pleasant_taste;pexp_not_pleasant_taste;pexp_HC_caption;pexp_HC_control];
subject=categorical([[1:length(pexp_ebdm_high_rew) 1:length(pexp_ebdm_high_rew)]  26+[1:26 1:26]  52+[1:13 1:13]  65+[1:24 1:24]] )';
studyNumber = categorical(studyNumber);
t=table(pexp,condition,isPleasure,studyNumber,subject);
t_hotspots = t;

if make_plot_do_stats

%% plot results
figure;
subplot(1,3,1:2)
UnivarScatter([  [pexp_ebdm_high_rew]-[pexp_ebdm_low_rew] [pexp_outcome_gain_correct(group==2)]-[pexp_outcome_loss_incorrect(group==2)]    [pexp_pleasant_taste; nan(13,1)]-[pexp_not_pleasant_taste; nan(13,1)] [pexp_HC_caption-pexp_HC_control; nan(2,1)]])
set(gca,'XTickLabel',{'Large vs. Small Reward Cue','Gain vs. Loss Feedback','Pleasant vs. Neutral Taste','Caption Task vs. Control'})
hSc=findobj(gcf,'Type','scatter');
set(hSc,'MarkerEdgeAlpha',.25)
set(hSc,'MarkerFaceAlpha',.5)
for i=1:length(hSc);set(hSc(i),'XData',get(hSc(i),'XData')+.2);end

ylabel '\Delta Signature Response (Cosine Sim)'
xlabel Condition
subplot(1,3,3)
hold all;



ROC_hotspots(1) = roc_plot([pexp_ebdm_high_rew;pexp_ebdm_low_rew],[ones(26,1);zeros(26,1)],'twochoice','nofig','color',hSc(4).CData);
ROC_hotspots(2) = roc_plot([pexp_outcome_gain_correct(group==2);pexp_outcome_loss_incorrect(group==2)],[ones(length(pexp_outcome_gain_correct(group==2)),1); zeros(length(pexp_outcome_loss_incorrect(group==2)),1)],'twochoice','color',hSc(3).CData);
ROC_hotspots(3) = roc_plot([pexp_pleasant_taste;pexp_not_pleasant_taste],[ones(length(pexp_pleasant_taste),1); zeros(length(pexp_not_pleasant_taste),1)],'twochoice','color',hSc(2).CData);
ROC_hotspots(4) = roc_plot([pexp_HC_caption; pexp_HC_control],[ones(length(pexp_HC_caption),1);zeros(length(pexp_HC_caption),1)],'twochoice','color',hSc(1).CData);
% save ROC_hotspots_hotspots_only ROC_hotspots

hSc=findobj(gca,'Type','line');

for i=1:length(hSc)
set(hSc(i),'Color',[get(hSc(i),'Color') .5]);
end


axis square
xlabel '1 - Specificity'
ylabel 'Sensitivity'

hL=findobj(gca,'Type','line');
set(hL(2:2:9),'Marker','.','MarkerSize',20)
legend(hL(1:2:8),fliplr({'Large vs. Small Reward','Gain vs. Loss Feedback','Pleasant vs. Neutral Taste','Captioning vs. Control'}))

for l=1:length(hL)
    try
        drawnow;
        hMarkers = hL(l).MarkerHandle;

        hMarkers.FaceColorType = 'truecoloralpha';

        hMarkers.EdgeColorData = uint8(double(hMarkers.EdgeColorData).*[1 1 1 .5]');
        hMarkers.FaceColorData = hMarkers.EdgeColorData;

    end
end


set(gca,'FontSize',12)


%% ANOVA w/VF
condition = categorical([ones(length(pexp_ebdm_high_rew),1);-1*ones(length(pexp_ebdm_high_rew),1);ones(length(pexp_outcome_gain_correct(group==2)),1);-1*ones(length(pexp_outcome_gain_correct(group==2)),1); ones(length(pexp_not_pleasant_taste),1); -1*ones(length(pexp_not_pleasant_taste),1); ones(length(pexp_HC_caption),1); -1*ones(length(pexp_HC_caption),1)]);
isPleasure = categorical([-1*ones(length(pexp_ebdm_high_rew),1);-1*ones(length(pexp_ebdm_high_rew),1);-1*ones(length(pexp_outcome_gain_correct(group==2)),1);-1*ones(length(pexp_outcome_gain_correct(group==2)),1); ones(length(pexp_not_pleasant_taste),1); ones(length(pexp_not_pleasant_taste),1); ones(length(pexp_HC_caption),1); ones(length(pexp_HC_caption),1)]);
studyNumber = [ones(2*length(pexp_ebdm_high_rew),1);2*ones(2*length(pexp_outcome_gain_correct(group==2)),1); 3*ones(2*length(pexp_not_pleasant_taste),1); 4*ones(2*length(pexp_HC_caption),1)];
pexp = [pexp_ebdm_high_rew;pexp_ebdm_low_rew; pexp_outcome_gain_correct(group==2);pexp_outcome_loss_incorrect(group==2);pexp_pleasant_taste;pexp_not_pleasant_taste;pexp_HC_caption;pexp_HC_control];
subject=categorical([[1:length(pexp_ebdm_high_rew) 1:length(pexp_ebdm_high_rew)]  26+[1:26 1:26]  52+[1:13 1:13]  65+[1:24 1:24]] )';
studyNumber = categorical(studyNumber);
t=table(pexp,condition,isPleasure,studyNumber,subject);

model = fitlme(t,'pexp~condition:isPleasure+studyNumber+(1|subject)');
[p,F,r]=coefTest(model)

%% perm tests
rng(0,'twister')

for nit = 1:3000

    gt = [ones(26,1);zeros(26,1)];
    gt=gt(randperm(2*26));
    tv = [pexp_ebdm_high_rew;pexp_ebdm_low_rew];
    null_ROC_hotspots(nit,1) = roc_plot(tv,gt,'twochoice','noplot','nooutput','noboot');
    null_es(nit,1) = compute_d_a(gt==1,tv);

    gt = [ones(26,1);zeros(26,1)];
    gt=gt(randperm(2*26));
   
    tv=[pexp_outcome_gain_correct(group==2);pexp_outcome_loss_incorrect(group==2)];
    null_ROC_hotspots(nit,2) = roc_plot(tv,gt,'twochoice','noplot','nooutput','noboot');
    null_es(nit,2) = compute_d_a(gt==1,tv);

    tv = [pexp_pleasant_taste;pexp_not_pleasant_taste];
    gt = [ones(13,1);zeros(13,1)];
    gt=gt(randperm(26));

    null_ROC_hotspots(nit,3) = roc_plot(tv,gt,'twochoice','noplot','nooutput','noboot');
    null_es(nit,3) = compute_d_a(gt==1,tv);

    tv = [pexp_HC_caption; pexp_HC_control];
    gt = [ones(24,1);zeros(24,1)];
    gt=gt(randperm(2*24));

    null_ROC_hotspots(nit,4) = roc_plot(tv,gt,'twochoice','noplot','nooutput','noboot');
    null_es(nit,4) = compute_d_a(gt==1,tv);

    
    nit/3000
end

% save null_ROC_hotspots_hotspots null_ROC_hotspots
phat(1)=1-(1+sum([null_ROC_hotspots(:,1).AUC]<ROC_hotspots(1).AUC))/(nit+1);
phat(2)=1-(1+sum([null_ROC_hotspots(:,2).AUC]<ROC_hotspots(2).AUC))/(nit+1);
phat(3)=1-(1+sum([null_ROC_hotspots(:,3).AUC]<ROC_hotspots(3).AUC))/(nit+1);
phat(4)=1-(1+sum([null_ROC_hotspots(:,4).AUC]<ROC_hotspots(4).AUC))/(nit+1);


%% ANOVA w/VF


%parametric bootstrap
model_orig = fitglme(t,'pexp ~ condition:isPleasure+studyNumber+(1|subject:studyNumber)+(1|subject)');
rng(0,'twister')
for it=1:10000
ynew = random(model_orig);
model = refit(model_orig,ynew);
bs_c(it,:) = model.Coefficients.Estimate;
it/10000
end

bs_P = 2*normcdf(-1*abs(mean(bs_c)./std(bs_c)),0,1)

for i=1:5
bs_ci(i,:)=[2*model_orig.Coefficients.Estimate(i)-prctile(bs_c(:,i),[97.5 2.5])];
end
end