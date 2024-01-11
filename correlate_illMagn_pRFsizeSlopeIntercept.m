%% Correlate illusion magnitudes and slope/intercept of pRF size as a function of eccentricity
clear all; clc; close all;

mainpath=fileparts(matlab.desktop.editor.getActiveFilename);

illMagn     = readtable(fullfile(mainpath,'ill_fMRI_illusionMagnitudes.csv'));
illMagn_num = illMagn.Variables;
illMagn     = illMagn.Properties.VariableNames;

load(fullfile(mainpath,'pRF_measures.mat'))

pltFig = true;
olt = 3.5;
modified_z = true; % if true, use modified zscore (Iglewicz & Hoaglin); if false, then use standard zscore

%% Calculate slope and intercept - pRF size as a function of eccentricity
subjects        = cellstr(unique(T_pRFs.subject,'rows'))';
ROIs            = {'V1','V2','V3','V4'};
hems            = {'lh','rh'};
eccentricities  = round(cell2mat(T_pRFs.eccentricity(1,:)),1);
ecc_1_8        = find(eccentricities >= 1 & eccentricities <= 8);

for sub = 1:length(subjects)
    if pltFig
        b = figure('innerposition',[0 0 2000 500],'units','centimeters');
    end
    i = 1;
    for hem = 1:length(hems)
        for ROI = 1:length(ROIs)
            roi_idx = strcmp(T_pRFs.ROI,ROIs{ROI});
            roi_idx = find(roi_idx);
            hem_idx = all(ismember(T_pRFs.hemisphere,hems{hem}),2);
            hem_idx = find(hem_idx);
            
            sub_idx = (T_pRFs.subject(:,1) == subjects{sub}(:,1)) & (T_pRFs.subject(:,2) == subjects{sub}(:,2));
            sub_idx = find(sub_idx);
            
            Sigma = cell2mat(T_pRFs.sigma(intersect(intersect(roi_idx,hem_idx),sub_idx),:));
            r2s   = cell2mat(T_pRFs.r2(intersect(intersect(roi_idx,hem_idx),sub_idx),:));
            
            % plot weighted linear regression line
            mdl = fitlm(eccentricities(ecc_1_8),Sigma(ecc_1_8),'Weights',r2s(ecc_1_8));

            Slope(sub,hem,ROI) = mdl.Coefficients.Estimate(2);
            Intercept(sub,hem,ROI) = mdl.Coefficients.Estimate(1);
            r2_adj(sub,hem,ROI) = mdl.Rsquared.Adjusted;
            r2(sub,hem,ROI) = mdl.Rsquared.Ordinary;
            if pltFig
                ax1 = subplot(2,length(ROIs),i);
                h1 = plot(mdl);
                h1(1).Color = 'k';
                h1(1).Marker = 'o';
                h1(1).MarkerFaceColor = 'k';
                h1(1).MarkerSize = 2;
                h1(2).LineWidth = 1.5;
                h1(3).LineWidth = 1.2;
                h1(4).LineWidth = 1.2;
                
                legend('off')
                if ROI == 1
                    if hem == 1
                        ylabel(['left hem.'])
                    elseif hem == 2
                        ylabel(['right hem.'])
                    end
                end
                
                if hem == 1
                    title(ROIs{ROI})
                else
                    title('')
                end
                xlim([1 8])
                xlabel('')
                ylabel('')
                i = i+1;
                set(gca,'fontname','arial','fontsize',12)
            end
        end
    end
    if pltFig
        if ~exist(fullfile(mainpath,'figures'),'dir')
            mkdir(fullfile(mainpath,'figures'))
        end
        saveas(b,fullfile(mainpath,'figures',['sub' num2str(sub) '.png']))
    end
end

r2(Slope<=0 | Intercept==0)=NaN;
Slope(Slope==0)     = NaN;
Intercept(Slope<0)  = NaN;
Slope(Slope<0)      = NaN;
Slope               = squeeze(nanmean(Slope,2));
r2 = squeeze(nanmean(r2,2));
Intercept(Intercept==0) = NaN;
Intercept = squeeze(nanmean(Intercept,2));


% remove illusion magn outliers
for ill = 1:13
    subs = 1:30;
    if modified_z
        M_z = modified_zscore(illMagn_num(subs,ill));
        subs = subs((abs(M_z) < olt));
    else
        z = zscore(illMagn_num(subs,ill));
        subs = subs((abs(z) < olt));
    end
    Illusions(:,ill) = illMagn_num(:,ill);
    Illusions(setdiff(1:30,subs),ill)  = NaN;
end

%% write table
clear T
% slope and intercept
for i = 1:length(ROIs)
    % slope
    v = genvarname([ROIs{i} '_slope']);
    eval([v '= squeeze(Slope(:,i));']);
    if i == 1
        eval(['T = table(' v ');'])
    else
        eval(['T=[T table(' v ')];'])
    end
    eval(['clear ' v])
    
    % intercept
    v = genvarname([ROIs{i} '_intercept']);
    eval([v '= squeeze(Intercept(:,i));']);
    eval(['T=[T table(' v ')];'])
    eval(['clear ' v])
    % r2
    v = genvarname([ROIs{i} '_r2']);
    eval([v '= squeeze(r2(:,i));']);
    eval(['T=[T table(' v ')];'])
    eval(['clear ' v])
end

% illusion
for i = 1:size(illMagn,2)
    v = genvarname([illMagn{1,i}]);
    eval([v '= Illusions(:,i);']);
    eval(['T=[T table(' v ')];'])
    eval(['clear ' v])
end
writetable(T,'illMagn_pRFsizeSlopeIntercept.csv')