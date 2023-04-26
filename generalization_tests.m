make_plot_do_stats = true; %default: true
exclude_outliers = false; %default: true
do_perms = false; %default: true
test_reward = false; % default: false

load([datadir filesep 'groupvars.mat']); %find control

load([tempdir filesep 'pos_affect_pattern.mat']); %pleasure signature from this study

brs = fmri_data([datadir filesep 'unthresholded_Z_map_1.nii.gz']); %reward signature from https://neurovault.org/images/775976/; see https://doi.org/10.1016/j.neuroimage.2023.119990

load([tempdir  filesep 'bpls.mat']) %betas for full model, including other domains

if test_reward
    pos_affect_pattern = brs;
end
clear ROC
for mdl = 2 %[3:5 2] %2 is the pleasure signature; 3-5 are other domains

    if ~test_reward
        pos_affect_pattern.dat = bpls(2:end,mdl);
    end

    % load data from Dalenberg et al. 2017 (validation study 1)
    P = {};
    files = dir([datadir filesep 'taste' filesep '**' filesep 'con_0001.nii']);
    for f=1:length(files)
        P{f}=[files(f).folder filesep files(f).name];
    end
    pleasant_taste = preprocess(fmri_data(P),'smooth',6); %add light smoothing
    outliers_pleasant_taste = plot(pleasant_taste,'norunmontages' );

    P = {};
    files = dir([datadir filesep 'taste' filesep '**' filesep 'con_0002.nii']);
    for f=1:length(files)
        P{f}=[files(f).folder filesep files(f).name];
    end
    not_pleasant_taste = preprocess(fmri_data(P),'smooth',6);  %add light smoothing

    %identify potential outliers using mahalanobis distance
    outliers_not_pleasant_taste = plot(not_pleasant_taste,'norunmontages' );

    pexp_pleasant_taste = apply_mask(pleasant_taste,pos_affect_pattern,'pattern_expression','cosine_similarity');
    pexp_not_pleasant_taste = apply_mask(not_pleasant_taste,pos_affect_pattern,'pattern_expression','cosine_similarity');


    % load data from Admon et al. 2015 (validation study 2)

    files = dir([datadir filesep 'Humor' filesep 'ISR2' filesep 'Caption_Reward-Control' filesep 'HC' filesep '*.img']);
    HC_files = strrep({files(:).name},'con_0006','con_0002');
    cd([datadir filesep 'Humor' filesep 'ISR2' filesep 'Caption_Reward' ]);
    dat_HC_caption = fmri_data(HC_files);
    cd([datadir filesep 'Humor' filesep 'ISR2' filesep 'Caption_Reward-Control' filesep 'HC' ]);
    dat_HC_caption_vs_control = fmri_data({files(:).name});
    cd(codedir)

    dat_control = dat_HC_caption_vs_control;
    dat_control.dat=dat_control.dat*-1;
    dat_control.dat = dat_control.dat + dat_HC_caption.dat;

    %identify potential outliers using mahalanobis distance
    outliers_caption = plot(dat_HC_caption,'norunmontages' );
    outliers_control = plot(dat_control,'norunmontages' );

    pexp_HC_caption = apply_mask(dat_HC_caption,pos_affect_pattern,'pattern_expression','cosine_similarity');
    pexp_HC_control = apply_mask(dat_control,pos_affect_pattern,'pattern_expression','cosine_similarity');

    % load data from Arulpragasm et al. 2018 (validation study 3)

    EBDM_Cue2_hrew = fmri_data([datadir filesep 'EBDM_cue2_low_effort_high_reward.nii']);
    EBDM_Cue2_lrew = fmri_data([datadir filesep 'EBDM_cue2_low_effort_low_reward.nii']);

    %identify potential outliers using mahalanobis distance
    outliers_ebdm_high_rew = plot(EBDM_Cue2_hrew,'norunmontages' );
    outliers_ebdm_low_rew = plot(EBDM_Cue2_lrew,'norunmontages' );

    %compute pattern expression
    pexp_ebdm_high_rew = apply_mask(EBDM_Cue2_hrew,pos_affect_pattern,'pattern_expression','cosine_similarity');
    pexp_ebdm_low_rew = apply_mask(EBDM_Cue2_lrew,pos_affect_pattern,'pattern_expression','cosine_similarity');

    % load data from  validation study 4

    load([datadir filesep 'outcomes.mat'])

    %identify potential outliers using mahalanobis distance
    outliers_loss = plot(outcome_loss_incorrect,'norunmontages' );
    outliers_gain = plot(outcome_gain_correct,'norunmontages' );

    %compute pattern expression
    pexp_outcome_loss_incorrect = apply_mask(outcome_loss_incorrect,pos_affect_pattern,'pattern_expression','cosine_similarity');
    pexp_outcome_gain_correct = apply_mask(outcome_gain_correct,pos_affect_pattern,'pattern_expression','cosine_similarity');



    %create variables for linear mixed effects model
    condition = categorical([ones(length(pexp_ebdm_high_rew),1);-1*ones(length(pexp_ebdm_high_rew),1);ones(length(pexp_outcome_gain_correct(group==2)),1);-1*ones(length(pexp_outcome_gain_correct(group==2)),1); ones(length(pexp_not_pleasant_taste),1); -1*ones(length(pexp_not_pleasant_taste),1); ones(length(pexp_HC_caption),1); -1*ones(length(pexp_HC_caption),1)]);
    isPleasure = categorical([-1*ones(length(pexp_ebdm_high_rew),1);-1*ones(length(pexp_ebdm_high_rew),1);-1*ones(length(pexp_outcome_gain_correct(group==2)),1);-1*ones(length(pexp_outcome_gain_correct(group==2)),1); ones(length(pexp_not_pleasant_taste),1); ones(length(pexp_not_pleasant_taste),1); ones(length(pexp_HC_caption),1); ones(length(pexp_HC_caption),1)]);
    studyNumber = [ones(2*length(pexp_ebdm_high_rew),1);2*ones(2*length(pexp_outcome_gain_correct(group==2)),1); 3*ones(2*length(pexp_not_pleasant_taste),1); 4*ones(2*length(pexp_HC_caption),1)];
    pexp = [pexp_ebdm_high_rew;pexp_ebdm_low_rew; pexp_outcome_gain_correct(group==2);pexp_outcome_loss_incorrect(group==2);pexp_pleasant_taste;pexp_not_pleasant_taste;pexp_HC_caption;pexp_HC_control];
    outliers = [outliers_ebdm_high_rew;outliers_ebdm_low_rew;outliers_gain(group==2);outliers_loss(group==2);outliers_pleasant_taste;outliers_not_pleasant_taste;outliers_caption;outliers_control];
    subject=categorical([[1:length(pexp_ebdm_high_rew) 1:length(pexp_ebdm_high_rew)]  26+[1:26 1:26]  52+[1:13 1:13]  65+[1:24 1:24]] )';
    studyNumber = categorical(studyNumber);

    %put data into a table for analysis
    t=table(pexp,condition,isPleasure,studyNumber,subject);

    if exclude_outliers
        pexp(outliers ==1) =NaN;
    end

    t_signature = t; %rename table to save output later

    if make_plot_do_stats
        %% plot results
        figure;
        subplot(1,3,1:2)
        UnivarScatter([  pexp(studyNumber==categorical(1) & condition==categorical(1)) - pexp(studyNumber==categorical(1) & condition==categorical(-1)) pexp(studyNumber==categorical(2)& condition==categorical(1)) - pexp(studyNumber==categorical(2)& condition==categorical(-1))    [pexp(studyNumber==categorical(3)& condition==categorical(1)) - pexp(studyNumber==categorical(3)& condition==categorical(-1)); nan(13,1)] [pexp(studyNumber==categorical(4)& condition==categorical(1)) - pexp(studyNumber==categorical(4)& condition==categorical(-1)); nan(2,1)]])
        set(gca,'XTickLabel',{'Large vs. Small Reward Cue','Gain vs. Loss Feedback','Pleasant vs. Neutral Taste','Caption Task vs. Control'})
        hSc=findobj(gcf,'Type','scatter');
        set(hSc,'MarkerEdgeAlpha',.25)
        set(hSc,'MarkerFaceAlpha',.5)
        for i=1:length(hSc);set(hSc(i),'XData',get(hSc(i),'XData')+.2);end

        ylabel '\Delta Signature Response (Cosine Sim)'
        xlabel Condition
        subplot(1,3,3)
        hold all;

        for ss=1:4
            tv = [pexp(studyNumber==categorical(ss) & condition==categorical(1)); pexp(studyNumber==categorical(ss) & condition==categorical(-1))];
            tv = reshape(tv,length(tv)/2,2);
            tv = tv(~any(isnan(tv')),:);
            tv = tv(:);
            ROC(ss) = roc_plot(tv,[ones(length(tv)/2,1);zeros(length(tv)/2,1)],'twochoice','nofig','color',hSc(5-ss).CData);
        end

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
        drawnow;snapnow;

        condition = categorical([ones(length(pexp_ebdm_high_rew),1);-1*ones(length(pexp_ebdm_high_rew),1);ones(length(pexp_outcome_gain_correct(group==2)),1);-1*ones(length(pexp_outcome_gain_correct(group==2)),1); ones(length(pexp_not_pleasant_taste),1); -1*ones(length(pexp_not_pleasant_taste),1); ones(length(pexp_HC_caption),1); -1*ones(length(pexp_HC_caption),1)]);
        isPleasure = categorical([-1*ones(length(pexp_ebdm_high_rew),1);-1*ones(length(pexp_ebdm_high_rew),1);-1*ones(length(pexp_outcome_gain_correct(group==2)),1);-1*ones(length(pexp_outcome_gain_correct(group==2)),1); ones(length(pexp_not_pleasant_taste),1); ones(length(pexp_not_pleasant_taste),1); ones(length(pexp_HC_caption),1); ones(length(pexp_HC_caption),1)]);
        studyNumber = [ones(2*length(pexp_ebdm_high_rew),1);2*ones(2*length(pexp_outcome_gain_correct(group==2)),1); 3*ones(2*length(pexp_not_pleasant_taste),1); 4*ones(2*length(pexp_HC_caption),1)];
        pexp = [pexp_ebdm_high_rew;pexp_ebdm_low_rew; pexp_outcome_gain_correct(group==2);pexp_outcome_loss_incorrect(group==2);pexp_pleasant_taste;pexp_not_pleasant_taste;pexp_HC_caption;pexp_HC_control];
        subject=categorical([[1:length(pexp_ebdm_high_rew) 1:length(pexp_ebdm_high_rew)]  26+[1:26 1:26]  52+[1:13 1:13]  65+[1:24 1:24]] )';
        studyNumber = categorical(studyNumber);
        t=table(pexp,condition,isPleasure,studyNumber,subject);

        model = fitlme(t,'pexp~condition:isPleasure+studyNumber+(1|subject)');
        [p,F,r]=coefTest(model);
    end
    AUC_by_model(mdl-1,:)=[ROC(:).AUC];

    %% perm tests

    if do_perms
        rng(0,'twister'); %for reproducibility

        for nit = 1:10000
            try
                gt = [ones(26,1);zeros(26,1)];
                gt=gt(randperm(2*26));
                tv = [pexp_ebdm_high_rew;pexp_ebdm_low_rew];
                null_ROC(nit,1) = roc_plot(tv,gt,'twochoice','noplot','nooutput','noboot');
                null_es(nit,1) = compute_d_a(gt==1,tv);

                gt = [ones(26,1);zeros(26,1)];
                gt=gt(randperm(2*26));

                tv=[pexp_outcome_gain_correct(group==2);pexp_outcome_loss_incorrect(group==2)];
                null_ROC(nit,2) = roc_plot(tv,gt,'twochoice','noplot','nooutput','noboot');
                null_es(nit,2) = compute_d_a(gt==1,tv);

                tv = [pexp_pleasant_taste;pexp_not_pleasant_taste];
                gt = [ones(13,1);zeros(13,1)];
                gt=gt(randperm(26));

                null_ROC(nit,3) = roc_plot(tv,gt,'twochoice','noplot','nooutput','noboot');
                null_es(nit,3) = compute_d_a(gt==1,tv);

                tv = [pexp_HC_caption; pexp_HC_control];
                gt = [ones(24,1);zeros(24,1)];
                gt=gt(randperm(2*24));

                null_ROC(nit,4) = roc_plot(tv,gt,'twochoice','noplot','nooutput','noboot');
                null_es(nit,4) = compute_d_a(gt==1,tv);

            catch
                nit=nit-1;
            end
            nit/10000
        end

        phat(mdl-1,1)=1-(1+sum([null_ROC(:,1).AUC]<ROC(1).AUC))/(nit+1);
        phat(mdl-1,2)=1-(1+sum([null_ROC(:,2).AUC]<ROC(2).AUC))/(nit+1);
        phat(mdl-1,3)=1-(1+sum([null_ROC(:,3).AUC]<ROC(3).AUC))/(nit+1);
        phat(mdl-1,4)=1-(1+sum([null_ROC(:,4).AUC]<ROC(4).AUC))/(nit+1);

        pleas_coef = model.Coefficients.Estimate(5);
    end

end

%parametric bootstrap
model_orig = fitglme(t,'pexp ~ condition:isPleasure+studyNumber+(1|subject:studyNumber)+(1|subject)');
rng(0,'twister'); %for reproducibility
for it=1:10000
    ynew = random(model_orig);
    model = refit(model_orig,ynew);
    bs_c(it,:) = model.Coefficients.Estimate;
end
[p,F,r]=coefTest(model);

bs_P = 2*normcdf(-1*abs(mean(bs_c)./std(bs_c)),0,1);

for i=1:5
    bs_ci(i,:)=[2*model_orig.Coefficients.Estimate(i)-prctile(bs_c(:,i),[97.5 2.5])];
end

