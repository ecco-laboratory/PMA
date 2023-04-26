% average signature weights within each region and test generalizability
load_avg_map = true;

if load_avg_map
avg_pos_affect_pattern = fmri_data(which('avg_pos_affect_pattern.nii'));
else
Amygdala = select_atlas_subset(load_atlas('Canlab2018'),{'Amy'});
subcort=merge_atlases(load_atlas('CIT168'),Amygdala);
insula = select_atlas_subset(load_atlas('glasser'),{'Ctx_PoI1_L','Ctx_Ig_L','Ctx_PI_L','Ctx_POI2_L','Ctx_MI_L','Ctx_AI_L','Ctx_AVI_L','Ctx_AAIC_L','Ctx_PoI1_R','Ctx_Ig_R','Ctx_PI_R','Ctx_POI2_R','Ctx_MI_R','Ctx_AI_R','Ctx_AVI_R','Ctx_AAIC_R'});
insula = resample_space(insula,subcort);

MFC=fmri_data(which('allregions.nii'));
vmpfc = resample_space(fmri_data(which('vmpfc.nii')),subcort);
dmpfc = resample_space(fmri_data(which('dmpfc.nii')),subcort);
sgACC = resample_space(fmri_data(which('sgACC.nii')),subcort);
pMCC = resample_space(fmri_data(which('pMCC.nii')),subcort);
pgACC =  resample_space(fmri_data(which('pgACC.nii')),subcort);
aMCC =  resample_space(fmri_data(which('aMCC.nii')),subcort);

all_regs_by_roi=vmpfc;
all_regs_by_roi=resample_space(all_regs_by_roi,subcort);
all_regs_by_roi.dat(all_regs_by_roi.dat>0)=1;
all_regs_by_roi.dat(dmpfc.dat>0)=max(all_regs_by_roi.dat)+ceil(dmpfc.dat(dmpfc.dat>0));
all_regs_by_roi.dat(sgACC.dat>0)=max(all_regs_by_roi.dat)+ceil(sgACC.dat(sgACC.dat>0));
all_regs_by_roi.dat(pMCC.dat>0)=max(all_regs_by_roi.dat)+ceil(pMCC.dat(pMCC.dat>0));
all_regs_by_roi.dat(pgACC.dat>0)=max(all_regs_by_roi.dat)+ceil(pgACC.dat(pgACC.dat>0));
all_regs_by_roi.dat(aMCC.dat>0)=max(all_regs_by_roi.dat)+ceil(aMCC.dat(aMCC.dat>0));

all_regs_by_roi.dat(subcort.dat>0)=max(all_regs_by_roi.dat)+subcort.dat(subcort.dat>0);
all_regs_by_roi.dat(insula.dat>0)=max(all_regs_by_roi.dat)+insula.dat(insula.dat>0);




load([resultsdir filesep 'pos_affect_pattern'])
all_regs_by_roi = resample_space(all_regs_by_roi,pos_affect_pattern,'nearest');
all_regs_by_roi = replace_empty(all_regs_by_roi);
pos_affect_pattern = replace_empty(pos_affect_pattern);

avg_pos_affect_pattern = replace_empty(pos_affect_pattern);
avg_pos_affect_pattern.dat = zeros(size(avg_pos_affect_pattern.dat));
for i = 1:max(all_regs_by_roi.dat)
    avg_pos_affect_pattern.dat(all_regs_by_roi.dat==i) = mean(pos_affect_pattern.dat(all_regs_by_roi.dat==i));
end

avg_pos_affect_pattern.fullpath = 'avg_pos_affect_pattern.nii';
% avg_pos_affect_pattern.write;
end

%%
clear P 
files = dir('C:\Users\phili\ds000219-download\**\con_0001.nii');
for f=1:length(files); P{f}=[files(f).folder filesep files(f).name];end
pleasant_taste = preprocess(fmri_data(P),'smooth',6);
files = dir('C:\Users\phili\ds000219-download\**\con_0002.nii');
clear P 
for f=1:length(files); P{f}=[files(f).folder filesep files(f).name];end
not_pleasant_taste = preprocess(fmri_data(P),'smooth',6);

pexp_pleasant_taste = apply_mask(pleasant_taste,avg_pos_affect_pattern,'pattern_expression','cosine_similarity');
pexp_not_pleasant_taste = apply_mask(not_pleasant_taste,avg_pos_affect_pattern,'pattern_expression','cosine_similarity');

EBDM_Cue1_hrew = fmri_data([datadir filesep 'EBDM_Cue2_low_effort_high_reward.nii']);
EBDM_Cue1_lrew = fmri_data([datadir filesep 'EBDM_Cue2_low_effort_low_reward.nii']);
EBDM_Cue1_hef = fmri_data([datadir filesep 'EBDM_cue1_high_effort.nii']);
EBDM_Cue1_lef = fmri_data([datadir filesep 'EBDM_cue1_low_effort.nii']);

pexp_ebdm_high_rew = apply_mask(EBDM_Cue1_hrew,avg_pos_affect_pattern,'pattern_expression','cosine_similarity');
pexp_ebdm_low_rew = apply_mask(EBDM_Cue1_lrew,avg_pos_affect_pattern,'pattern_expression','cosine_similarity');

pexp_ebdm_high_eff = apply_mask(EBDM_Cue1_hef,avg_pos_affect_pattern,'pattern_expression','cosine_similarity');
pexp_ebdm_low_eff = apply_mask(EBDM_Cue1_lef,avg_pos_affect_pattern,'pattern_expression','cosine_similarity');

load('C:\Users\phili\OneDrive - Emory University\Documents\DeepRL_Valence\SavedCode\outcomes.mat')
pexp_outcome_loss_incorrect = apply_mask(outcome_loss_incorrect,avg_pos_affect_pattern,'pattern_expression','cosine_similarity');
pexp_outcome_gain_correct = apply_mask(outcome_gain_correct,avg_pos_affect_pattern,'pattern_expression','cosine_similarity');

    files = dir([datadir filesep 'Humor' filesep 'ISR2' filesep 'Caption_Reward-Control' filesep 'HC' filesep '*.img']);
HC_files = strrep({files(:).name},'con_0006','con_0002');
cd ([datadir filesep 'Humor' filesep 'ISR2' filesep 'Caption_Reward' ]);
dat_HC_caption = fmri_data(HC_files);
cd ([datadir filesep 'Humor' filesep 'ISR2' filesep 'Caption_Reward-Control' filesep 'HC' ]);
dat_HC_caption_vs_control = fmri_data({files(:).name});

dat_control = dat_HC_caption_vs_control;
dat_control.dat=dat_control.dat*-1;
dat_control.dat = dat_control.dat + dat_HC_caption.dat;

pexp_HC_caption = apply_mask(dat_HC_caption,avg_pos_affect_pattern,'pattern_expression','cosine_similarity');
pexp_HC_control = apply_mask(dat_control,avg_pos_affect_pattern,'pattern_expression','cosine_similarity');

%% plot results
load([datadir filesep 'groupvars.mat']);clear ROC;
figure;
subplot(1,3,1:2)
UnivarScatter([  [pexp_ebdm_high_rew]-[pexp_ebdm_low_rew] [pexp_outcome_gain_correct(group==2)]-[pexp_outcome_loss_incorrect(group==2)]    [pexp_pleasant_taste; nan(13,1)]-[pexp_not_pleasant_taste; nan(13,1)] [pexp_HC_caption-pexp_HC_control; nan(2,1)]])
set(gca,'XTickLabel',{'Large vs. Small Reward Cue','Gain vs. Loss Feedback','Pleasant vs. Neutral Taste','Caption Task vs. Control'})
hSc=findobj(gcf,'Type','scatter');
set(hSc,'MarkerEdgeAlpha',.25)
set(hSc,'MarkerFaceAlpha',.5)
for i=1:length(hSc);set(hSc(i),'XData',get(hSc(i),'XData')+.2);end

ylabel '\Delta Pattern Expression (Cosine Sim)'
xlabel Condition
subplot(1,3,3)
hold all;



ROC(1) = roc_plot([pexp_ebdm_high_rew;pexp_ebdm_low_rew],[ones(26,1);zeros(26,1)],'twochoice','nofig','color',hSc(4).CData);
ROC(2) = roc_plot([pexp_outcome_gain_correct(group==2);pexp_outcome_loss_incorrect(group==2)],[ones(length(pexp_outcome_gain_correct(group==2)),1); zeros(length(pexp_outcome_loss_incorrect(group==2)),1)],'twochoice','color',hSc(3).CData);
ROC(3) = roc_plot([pexp_pleasant_taste;pexp_not_pleasant_taste],[ones(length(pexp_pleasant_taste),1); zeros(length(pexp_not_pleasant_taste),1)],'twochoice','color',hSc(2).CData);
ROC(4) = roc_plot([pexp_HC_caption; pexp_HC_control],[ones(length(pexp_HC_caption),1);zeros(length(pexp_HC_caption),1)],'twochoice','color',hSc(1).CData);


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
isPleasure = categorical([-1*ones(length(pexp_ebdm_high_rew),1);-1*ones(length(pexp_ebdm_high_rew),1);-1*ones(length(pexp_outcome_gain_correct(group==2)),1);-1*ones(length(pexp_outcome_gain_correct(group==2)),1); ones(length(pexp_not_pleasant_taste),1); -1*ones(length(pexp_not_pleasant_taste),1); ones(length(pexp_HC_caption),1); -1*ones(length(pexp_HC_caption),1)]);
studyNumber = [ones(2*length(pexp_ebdm_high_rew),1);2*ones(2*length(pexp_outcome_gain_correct(group==2)),1); 3*ones(2*length(pexp_not_pleasant_taste),1); 4*ones(2*length(pexp_HC_caption),1)];
pexp = [pexp_ebdm_high_rew;pexp_ebdm_low_rew; pexp_outcome_gain_correct(group==2);pexp_outcome_loss_incorrect(group==2);pexp_pleasant_taste;pexp_not_pleasant_taste;pexp_HC_caption;pexp_HC_control];
subject=categorical([[1:length(pexp_ebdm_high_rew) 1:length(pexp_ebdm_high_rew)]  26+[1:26 1:26]  52+[1:13 1:13]  65+[1:24 1:24]] )';
studyNumber = categorical(studyNumber);
t=table(pexp,isPleasure,studyNumber,subject);

model = fitlme(t,'pexp~isPleasure+studyNumber+(1|subject)');

BF_no_pleasure = bf.anova(t,'pexp~studyNumber+(1|subject)');
BF_pleasure = bf.anova(t,'pexp~isPleasure+studyNumber+(1|subject)');

BF_pleasure/BF_no_pleasure

%% perm tests

for nit = 1:3000

    gt = [ones(26,1);zeros(26,1)];
    gt=gt(randperm(2*26));
    null_ROC(nit,1) = roc_plot([pexp_ebdm_high_rew;pexp_ebdm_low_rew],gt,'twochoice','noplot','nooutput','noboot');
    null_ROC(nit,2) = roc_plot([pexp_outcome_gain_correct(group==2);pexp_outcome_loss_incorrect(group==2)],gt,'twochoice','noplot','nooutput','noboot');

    gt = [ones(13,1);zeros(13,1)];
    gt=gt(randperm(26));
    null_ROC(nit,3) = roc_plot([pexp_pleasant_taste;pexp_not_pleasant_taste],gt,'twochoice','noplot','nooutput','noboot');


    gt = [ones(24,1);zeros(24,1)];
    gt=gt(randperm(2*24));
    null_ROC(nit,4) = roc_plot([pexp_HC_caption; pexp_HC_control],gt,'twochoice','noplot','nooutput','noboot');

    nit/3000
end

phat(1)=1-(1+sum([null_ROC(:,1).AUC]<ROC(1).AUC))/(nit+1);
phat(2)=1-(1+sum([null_ROC(:,2).AUC]<ROC(2).AUC))/(nit+1);
phat(3)=1-(1+sum([null_ROC(:,3).AUC]<ROC(3).AUC))/(nit+1);
phat(4)=1-(1+sum([null_ROC(:,4).AUC]<ROC(4).AUC))/(nit+1);


