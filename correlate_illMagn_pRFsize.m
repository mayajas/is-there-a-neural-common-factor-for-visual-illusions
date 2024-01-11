%% Correlate illusion magnitudes and pRF sizes at relevant eccentricities
clear all;clc; close all;

mainpath=fileparts(matlab.desktop.editor.getActiveFilename);

illMagn     = readtable(fullfile(mainpath,'ill_fMRI_illusionMagnitudes.csv'));
illMagn_num = illMagn.Variables;
illMagn     = illMagn.Properties.VariableNames;

load(fullfile(mainpath,'pRF_measures.mat'))

olt = 3.5;
modified_z = true; % if true, use modified zscore (Iglewicz & Hoaglin); if false, then use standard zscore
%% clean illusions
for ill = 1:13
    if modified_z
        M_z = modified_zscore(illMagn_num(:,ill));
        illMagn_num(~(abs(M_z) < olt),ill) = NaN;
    else
        z = zscore(illMagn_num(:,ill));
        illMagn_num(~(abs(z) < olt),ill)=NaN;
    end
end

%% Delboeuf illusion
% Eccentricities to consider:
% Top left 1.96°/3.92° - right ventral
% Bottom right 1.96°/3.92° - left dorsal
ill             = 3;
eoi             = [1.96 3.92]; % eccentricities of interest
subjects        = cellstr(unique(T_pRFs.subject,'rows'))';
ROIs            = {'V1v','V2v','V3v','V4','V1d','V2d','V3d'};
hems            = {'lh','rh'};
eccentricities  = round(cell2mat(T_pRFs.eccentricity(1,:)),1);
ecc_1_8        = find(eccentricities >= 1 & eccentricities <= 8);

clear Sigma 
for ROI = 1:length(ROIs)
    for sub = 1:length(subjects)
        for hem = 1:length(hems)
            roi_idx = strfind(T_pRFs.ROI,ROIs{ROI});
            roi_idx = find(~cellfun(@isempty,roi_idx));
            hem_idx = all(ismember(T_pRFs.hemisphere,hems{hem}),2);
            hem_idx = find(hem_idx);
            
            sub_idx = (T_pRFs.subject(:,1) == subjects{sub}(:,1)) & (T_pRFs.subject(:,2) == subjects{sub}(:,2));
            sub_idx = find(sub_idx);

            Sigma = cell2mat(T_pRFs.sigma(intersect(intersect(roi_idx,hem_idx),sub_idx),:));
            r2s   = cell2mat(T_pRFs.r2(intersect(intersect(roi_idx,hem_idx),sub_idx),:));
    
            % left dorsal rois
            if hem == 1 && (ROI == 5 || ROI == 6 || ROI == 7)
                Ld1(sub,ROI-4) = predict(fitlm(eccentricities(ecc_1_8),Sigma(ecc_1_8),'Weights',r2s(ecc_1_8)),eoi(1));
                Ld2(sub,ROI-4) = predict(fitlm(eccentricities(ecc_1_8),Sigma(ecc_1_8),'Weights',r2s(ecc_1_8)),eoi(2));
            end

            % right ventral rois
            if hem == 2 && (ROI == 1 || ROI == 2 || ROI == 3 || ROI == 4)
                Rv1(sub,ROI) = predict(fitlm(eccentricities(ecc_1_8),Sigma(ecc_1_8),'Weights',r2s(ecc_1_8)),eoi(1));
                Rv2(sub,ROI) = predict(fitlm(eccentricities(ecc_1_8),Sigma(ecc_1_8),'Weights',r2s(ecc_1_8)),eoi(2));
            end
        end
    end
end

b=figure('innerposition',[0 0 2000 500],'units','centimeters','name','Delboeuf');
LdRv1 = [nanmean(cat(3,Ld1,Rv1(:,1:3)),3) Rv1(:,4)];
LdRv2 = [nanmean(cat(3,Ld2,Rv2(:,1:3)),3) Rv2(:,4)];

% 1.96
for i=1:4
    ax1 = subplot(2,4,i);
    % remove outliers & NaNs
    subs = 1:30;
    if modified_z
        M_z = modified_zscore(LdRv1(:,i));
        LdRv1(~(abs(M_z) < olt),i) = NaN;
    else
        z = zscore(LdRv1(:,i));
        LdRv1(~(abs(z) < olt),i) = NaN;
    end
    subs = 1:30;
    subs = subs(~isnan(LdRv1(subs,i)) & ~isnan(illMagn_num(:,ill)));
    % save pRF sizes for table
    pRF_1_96(:,i)               = LdRv1(:,i);
    pRF_1_96(setdiff(1:30,subs),i)  = NaN;
    % scatter plot
    scatter(illMagn_num(subs,ill),LdRv1(subs,i),7,'filled','MarkerEdgeColor','k','MarkerFaceColor','k')
    [r,p]= corr(illMagn_num(subs,ill),LdRv1(subs,i));
    r_DB_1(i) = r; p_DB_1(i)=p;
    h1 = lsline(ax1);
    h1.Color = 'r';
    h1.LineWidth = 0.5;
    title(['r=' num2str(r,'%.3f') ', p=' num2str(p,'%.3f') ', N=' num2str(length(subs))])
    if i == 1
        ylabel('1.96 deg')
    end
    set(gca,'fontname','arial','fontsize',12)
    xlim([-0.13 0.33])

    if i==1
        ylim([0.25 1.05])
    elseif i==2
        ylim([0.4 1.6])
    elseif i==3
        ylim([0.75 1.9])
    elseif i==4
        ylim([1 3.15])
    end
end
% 3.92
for i=1:4
    ax1 = subplot(2,4,i+4);
    % remove outliers
    subs = 1:30;
    if modified_z
        M_z = modified_zscore(LdRv2(:,i));
        LdRv2(~(abs(M_z) < olt),i) = NaN;
    else
        z = zscore(LdRv2(:,i));
        LdRv2(~(abs(z) < olt),i) = NaN;
    end
    subs = 1:30;
    subs = subs(~isnan(LdRv2(subs,i)) & ~isnan(illMagn_num(:,ill)));
    % save pRF sizes for table
    pRF_3_92(:,i)               = LdRv2(:,i);
    pRF_3_92(setdiff(1:30,subs),i)  = NaN;
    % scatter plot
    scatter(illMagn_num(subs,ill),LdRv2(subs,i),7,'filled','MarkerEdgeColor','k','MarkerFaceColor','k')
    [r,p]= corr(illMagn_num(subs,ill),LdRv2(subs,i));
    r_DB_2(i) = r; p_DB_2(i)=p;
    h1 = lsline(ax1);
    h1.Color = 'r';
    h1.LineWidth = 0.5;
    title(['r=' num2str(r,'%.3f') ', p=' num2str(p,'%.3f') ', N=' num2str(length(subs))])
    if i == 1
        ylabel('3.92 deg')
    end
    set(gca,'fontname','arial','fontsize',12)
    xlim([-0.13 0.33])

    if i==1
        ylim([0.4 1.3])
    elseif i==2
        ylim([0.7 1.8])
    elseif i==3
        ylim([1.2 2.05])
    elseif i==4
        ylim([1.25 3.05])
    end
end
saveas(b,fullfile(mainpath,'figures','db.svg'))



%% Ebbinghaus 2 illusion
% Eccentricities to consider:
% Left 4.65°/9.3° 
% Right 4.65°/9.3° 
ill             = 5;
eoi             = [4.65 9.3]; % eccentricities of interest
ROIs            = {'V1','V2','V3','V4'};

clear Sigma
for ROI = 1:length(ROIs)
    for sub = 1:length(subjects)
        for hem = 1:length(hems)
            roi_idx = strfind(T_pRFs.ROI,ROIs{ROI});
            roi_idx = find(~cellfun(@isempty,roi_idx));
            hem_idx = all(ismember(T_pRFs.hemisphere,hems{hem}),2);
            hem_idx = find(hem_idx);
            
            sub_idx = (T_pRFs.subject(:,1) == subjects{sub}(:,1)) & (T_pRFs.subject(:,2) == subjects{sub}(:,2));
            sub_idx = find(sub_idx);

            Sigma = cell2mat(T_pRFs.sigma(intersect(intersect(roi_idx,hem_idx),sub_idx),:));
            r2s   = cell2mat(T_pRFs.r2(intersect(intersect(roi_idx,hem_idx),sub_idx),:));
    
            % left rois
            if hem == 1 
                L1(sub,ROI) = predict(fitlm(eccentricities(ecc_1_8),Sigma(ecc_1_8),'Weights',r2s(ecc_1_8)),eoi(1));
                L2(sub,ROI) = predict(fitlm(eccentricities(ecc_1_8),Sigma(ecc_1_8),'Weights',r2s(ecc_1_8)),eoi(2));
            end

            % right rois
            if hem == 2 
                R1(sub,ROI) = predict(fitlm(eccentricities(ecc_1_8),Sigma(ecc_1_8),'Weights',r2s(ecc_1_8)),eoi(1));
                R2(sub,ROI) = predict(fitlm(eccentricities(ecc_1_8),Sigma(ecc_1_8),'Weights',r2s(ecc_1_8)),eoi(2));
            end
        end
    end
end


b=figure('innerposition',[0 0 2000 500],'units','centimeters','name','Ebbinghaus 2');
LR1 = nanmean(cat(3,L1,R1),3);
LR2 = nanmean(cat(3,L2,R2),3);

% 4.65°
for i=1:4
    ax1 = subplot(2,4,i);
    % remove outliers & NaNs
    if modified_z
        M_z = modified_zscore(LR1(:,i));
        LR1(~(abs(M_z) < olt),i) = NaN;
    else
        z = zscore(LR1(:,i));
        LR1(~(abs(z) < olt),i) = NaN;
    end
    subs = 1:30;
    subs = subs(~isnan(LR1(subs,i)) & ~isnan(illMagn_num(:,ill)));
    % save pRF sizes for table
    pRF_4_65(:,i)               = LR1(:,i);
    pRF_4_65(setdiff(1:30,subs),i)  = NaN;
    % scatter plot
    scatter(illMagn_num(subs,ill),LR1(subs,i),7,'filled','MarkerEdgeColor','k','MarkerFaceColor','k')
    [r,p]= corr(illMagn_num(subs,ill),LR1(subs,i));
    r_EB2_1(i) = r; p_EB2_1(i)=p;
    h1 = lsline(ax1);
    h1.Color = 'r';
    h1.LineWidth = 0.5;
    title(['r=' num2str(r,'%.3f') ', p=' num2str(p,'%.3f') ', N=' num2str(length(subs))])
    if i == 1
        ylabel('4.65 deg')
    end
    set(gca,'fontname','arial','fontsize',12)
    xlim([-0.15 0.55])
    if i==1
        ylim([0.25 0.95])
    elseif i==2
        ylim([0.4 1.4])
    elseif i==3
        ylim([0.75 1.6])
    elseif i==4
        ylim([1.4 2.9])
    end
end
% 9.3°
for i=1:4
    ax1 = subplot(2,4,i+4);
    % remove outliers
    if modified_z
        M_z = modified_zscore(LR2(:,i));
        LR2(~(abs(M_z) < olt),i) = NaN;
    else
        z = zscore(LR2(:,i));
        LR2(~(abs(z) < olt),i) = NaN;
    end
    subs = 1:30;
    subs = subs(~isnan(LR2(subs,i)) & ~isnan(illMagn_num(:,ill)));
    % save pRF sizes for table
    pRF_9_3(:,i)               = LR2(:,i);
    pRF_9_3(setdiff(1:30,subs),i)  = NaN;
    % scatter plot
    scatter(illMagn_num(subs,ill),LR2(subs,i),7,'filled','MarkerEdgeColor','k','MarkerFaceColor','k')
    [r,p]= corr(illMagn_num(subs,ill),LR2(subs,i));
    r_EB2_2(i) = r; p_EB2_2(i)=p;
    h1 = lsline(ax1);
    h1.Color = 'r';
    h1.LineWidth = 0.5;
    title(['r=' num2str(r,'%.3f') ', p=' num2str(p,'%.3f') ', N=' num2str(length(subs))])
    if i == 1
        ylabel('9.3 deg')
    end
    set(gca,'fontname','arial','fontsize',12)
    xlim([-0.15 0.55])
    if i==1
        ylim([0.4 1.45])
    elseif i==2
        ylim([0.65 1.85])
    elseif i==3
        ylim([1.0 2.35])
    elseif i==4
        ylim([1.8 4])
    end
end
saveas(b,fullfile(mainpath,'figures','eb2.svg'))


%% write table
ROIs            = {'V1','V2','V3','V4'};
clear T
% slope and intercept
for i = 1:length(ROIs)
    % pRF small DB
    v = genvarname([ROIs{i} '_pRF_1_96']);
    eval([v '= squeeze(pRF_1_96(:,i));']);
    if i == 1
        eval(['T = table(' v ');'])
    else
        eval(['T=[T table(' v ')];'])
    end
    eval(['clear ' v])
        
    % pRF big DB
    v = genvarname([ROIs{i} '_pRF_3_92']);
    eval([v '= squeeze(pRF_3_92(:,i));']);
    eval(['T=[T table(' v ')];'])
    eval(['clear ' v])
    
    % pRF small EB2
    v = genvarname([ROIs{i} '_pRF_4_65']);
    eval([v '= squeeze(pRF_4_65(:,i));']);
    eval(['T=[T table(' v ')];'])
    eval(['clear ' v])
    
    % pRF big EB2
    v = genvarname([ROIs{i} '_pRF_9_3']);
    eval([v '= squeeze(pRF_9_3(:,i));']);
    eval(['T=[T table(' v ')];'])
    eval(['clear ' v])
end

% illusion magnitude
for i = [3 5] % db and eb2
    v = genvarname([illMagn{1,i}]);
    eval([v '= illMagn_num(:,i);']);
    eval(['T=[T table(' v ')];'])
    eval(['clear ' v])
end

writetable(T,'illMagn_pRF_size.csv')

