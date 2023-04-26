PMA = spm_read_vols(spm_vol([resultsdir '\pos_affect_pattern.nii']));

%% INS (sagittal) R
% figure;imagesc(flipud(squeeze(PMA(55,:,:))'))
% PMA(PMA==0)=nan;
figure;imagesc(flipud(squeeze(PMA(55,35:73,19:50))'))
colormap(canlabCmap)
caxis([-.0003 .0003])
colorbar;
set(gca,'visible','off')
drawnow;snapnow;
%% INS (sagittal) L

figure;imagesc(flipud(squeeze(PMA(20,35:73,19:50))'))
colormap(canlabCmap)
caxis([-.0003 .0003])
colorbar;
set(gca,'visible','off')
drawnow;snapnow;

%% cingulate/striatum (coronal)
% figure;imagesc(flipud(squeeze(PMA(:,60,:))')) %MCC VS
figure;imagesc(flipud(squeeze(PMA(:,70,:))')) %insula

colormap(canlabCmap)
caxis([-.0003 .0003])
colorbar;
set(gca,'visible','off')
drawnow;snapnow;


%% midbrain transverse

figure;imagesc(flipud(squeeze(PMA(:,:,27))'))
colormap(canlabCmap)
caxis([-.0003 .0003])
colorbar;
set(gca,'visible','off')
drawnow;snapnow;

%% aMCC

figure;imagesc(flipud(squeeze(PMA(40,:,:))'))
colormap(canlabCmap)
caxis([-.0003 .0003])
colorbar;
set(gca,'visible','off')
drawnow;snapnow;
