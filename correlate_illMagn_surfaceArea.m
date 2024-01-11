function correlate_illMagn_surfaceArea()
%% Correlate illusion magnitudes and V1-V4 surface area
clear all;clc; close all;
load('surface_areas_0_8.mat','Area_entire','ROIs_all')
T = table();
mainpath=fileparts(matlab.desktop.editor.getActiveFilename);

illMagn     = readtable(fullfile(mainpath,'ill_fMRI_illusionMagnitudes.csv'));
illMagn_num = illMagn.Variables;
illMagn     = illMagn.Properties.VariableNames;

pltFig = true; % whether to save figures

olt = 3.5;
modified_z = true; % if true, use modified zscore (Iglewicz & Hoaglin); if false, then use standard zscore

%% write table
Area_entire_bothhem = squeeze(sum(Area_entire,2));
for ii = 1:length(ROIs_all)
    eval(['T.' ROIs_all{ii} '= Area_entire_bothhem(:,ii);'])
    % left hem
    eval(['T.' ROIs_all{ii} '_SA_left = squeeze(Area_entire(:,1,ii));'])
    % right hem
    eval(['T.' ROIs_all{ii} '_SA_right = squeeze(Area_entire(:,2,ii));'])
end

for ii = 1:size(illMagn,2)
    eval(['T.' illMagn{1,ii} '= illMagn_num(:,ii);'])
end

%% Plots
roi_idx = [1 4 7 10]; % V1, V2, V3, V4 (other indices correspond to dorsal/ventral subregions)
EB_PZ = [3 5 11]; % illusions previously tested by Schwarzkopf et al.
other_size_illusions = [1 4 8 10]; % other size illusions
other_illusions = [2 6 7 9 12 13]; % other illusions of uniform texture, perceived orientation, and contrast

SA = surfaceAreaOutlierRemoval();
Illusions = illusionsOutlierRemoval();
rr = [];
pp = [];
surfaceArea_vs_illMagn(EB_PZ,'Schwarzkopf')
surfaceArea_vs_illMagn(other_size_illusions,'other_size')
surfaceArea_vs_illMagn(other_illusions,'other')
mean_explained_variance = mean(rr.^2);
for tt = 1:4
    disp([ROIs_all{tt} ': ' num2str(round(mean_explained_variance(tt)*100,1)) '% explained variance'])
end
rr = round(rr,3);
if pltFig
    corr_heatmap()
end

%% write table
clear T
T = table();
for ii = 1:length(ROIs_all)
    eval(['T.' ROIs_all{ii} '= SA(:,ii);'])    
end

for ii = 1:size(illMagn,2)
    eval(['T.' illMagn{1,ii} '= Illusions(:,ii);'])    
end
writetable(T,fullfile(mainpath,'illMagn_surfaceArea_0_8.csv')   )

allROIsModel_adjR2 = [];
for ii = 1:size(illMagn,2)
    lm = fitlm(T,[illMagn{ii} ' ~ V1 + V2 + V3 + V4']);
    allROIsModel_adjR2 = [allROIsModel_adjR2 lm.Rsquared.Adjusted];
end

%% functions
    function SA = surfaceAreaOutlierRemoval()
        for roi = 1:4
            subs = 1:30;
            if modified_z
                M_z = modified_zscore(table2array(T(subs,roi_idx(roi))));
                subs = subs((abs(M_z) < olt));
            else
                z = zscore(table2array(T(subs,roi_idx(roi))));
                subs = subs((abs(z) < olt));
            end
            % save SA for table
            SA(:,roi)                   = table2array(T(:,roi_idx(roi)));
            SA(setdiff(1:30,subs),roi)  = NaN;
        end
    end
    function Illusions = illusionsOutlierRemoval()
        % remove illusion magn outliers
        for ill = 1:13
            subs = 1:30;
            if modified_z
                M_z = modified_zscore(table2array(T(subs,ill+12)));
                subs = subs(abs(M_z) < olt);
            else
                z = zscore(table2array(T(subs,ill+12)));
                subs = subs((abs(z) < olt));
            end
            Illusions(:,ill) = table2array(T(:,ill+12));
            Illusions(setdiff(1:30,subs),ill)  = NaN;
        end
    end
    function surfaceArea_vs_illMagn(illVals,whichIlls)
        a = figure('innerposition',[0 0 2000 1000],'units','centimeters');
        i = 1;
        for roi = 1:4
            for ill = 1:length(illVals)
                subs = 1:30;
                subs = subs(~isnan(SA(:,roi)) & ~isnan(Illusions(:,illVals(ill))));
                
                [r,p]= corr(Illusions(subs,illVals(ill)),SA(subs,roi));
                rr(illVals(ill),roi) = r;
                pp(illVals(ill),roi) = p;
                
                ax1 = subplot(4,length(illVals),i);
                i = i+1;
                scatter(Illusions(subs,illVals(ill)),SA(subs,roi),7,'filled','MarkerEdgeColor','k','MarkerFaceColor','k')
                h1 = lsline(ax1);
                h1.Color = 'r';
                h1.LineWidth = 0.5;
                
                
                if ill == 1
                    ylabel([ROIs_all{roi}])
                else
                    yticklabels([])
                end
                if roi == 1
                    title({illMagn{illVals(ill)}, ['r=' num2str(r,'%.3f') ', p=' num2str(p,'%.3f') ', N=' num2str(length(subs))]})
                    xticklabels([])
                elseif roi < 4
                    xticklabels([])
                    title(['r=' num2str(r,'%.3f') ', p=' num2str(p,'%.3f') ', N=' num2str(length(subs))])
                else
                    title(['r=' num2str(r,'%.3f') ', p=' num2str(p,'%.3f') ', N=' num2str(length(subs))])
                end

                if strcmp(whichIlls,'other')
                    set(gca,'fontname','arial','fontsize',10)
                else
                    set(gca,'fontname','arial','fontsize',12)
                end
                
                subs = ~isnan(Illusions(:,illVals(ill)));
                xlim([h1.XData(1)-0.1*abs(h1.XData(2)) h1.XData(2)+0.1*abs(h1.XData(2))])

                if roi == 1
                    ylim([1200 3900])
                    yticks(1500:1000:4000) 
                elseif roi == 2
                    ylim([700 3500])
                    yticks(1000:1000:4000)
                elseif roi == 3
                    ylim([950 2500])
                    yticks(0:500:3250)
                elseif roi == 4
                    ylim([300 1300])
                    yticks(0:250:2000)
                end
                
                
            end
        end
        if pltFig
            if ~exist(fullfile(mainpath,'figures'),'dir')
                mkdir(fullfile(mainpath,'figures'))
            end
            if strcmp(whichIlls,'Schwarzkopf')
                saveas(a,fullfile(mainpath,'figures','_correlations_Schwarzkopf_et_al_withoutOutliers.svg'))
            elseif strcmp(whichIlls,'other_size')
                saveas(a,fullfile(mainpath,'figures','_correlations_size_other_withoutOutliers.svg'))
            elseif strcmp(whichIlls,'other')
                saveas(a,fullfile(mainpath,'figures','_correlations_other.svg'))
            end
        end
    end
    function corr_heatmap()
        a = figure;
        m = 64;
        mcc = 52;
        fntsz=16;
        if (mod(m,2) == 0)
            % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
            m1 = m*0.5;
            r = (0:m1-1)'/max(m1-1,1);
            g = r;
            r = [r; ones(m1,1)];
            g = [g; flipud(g)];
            b = flipud(r);
        else
            % From [0 0 1] to [1 1 1] to [1 0 0];
            m1 = floor(m*0.5);
            r = (0:m1-1)'/max(m1,1);
            g = r;
            r = [r; ones(m1+1,1)];
            g = [g; 1; flipud(g)];
            b = flipud(r);
        end
        c = [r g b];
        subplot(131)
        hh = heatmap(rr);
        hh.ColorLimits = [-1 1];
        hh.Colormap = c;
        hh.XDisplayLabels = ROIs_all;
        hh.YDisplayLabels = illMagn;
        title('After outlier removal')
        set(gca,'FontSize',fntsz)
        subplot(132)
        m = 64;
        if (mod(m,2) == 0)
            % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
            m1 = m*0.5;
            r = (0:m1-1)'/max(m1-1,1);
            g = r;
            r = [r; ones(m1,1)];
            g = [g; flipud(g)];
            b = flipud(r);
        else
            % From [0 0 1] to [1 1 1] to [1 0 0];
            m1 = floor(m*0.5);
            r = (0:m1-1)'/max(m1,1);
            g = r;
            r = [r; ones(m1+1,1)];
            g = [g; 1; flipud(g)];
            b = flipud(r);
        end
        c = [r g b];
        
        rr(pp>.05) = NaN;
        hh = heatmap(rr);
        hh.ColorLimits = [-1 1];
        hh.Colormap = c;
        hh.XDisplayLabels = ROIs_all;
        hh.YDisplayLabels = illMagn;
        title('Noncorrected significant')
        set(gca,'FontSize',fntsz)
        subplot(133)
        m = 64;
        if (mod(m,2) == 0)
            % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
            m1 = m*0.5;
            r = (0:m1-1)'/max(m1-1,1);
            g = r;
            r = [r; ones(m1,1)];
            g = [g; flipud(g)];
            b = flipud(r);
        else
            % From [0 0 1] to [1 1 1] to [1 0 0];
            m1 = floor(m*0.5);
            r = (0:m1-1)'/max(m1,1);
            g = r;
            r = [r; ones(m1+1,1)];
            g = [g; 1; flipud(g)];
            b = flipud(r);
        end
        c = [r g b];
        rr(pp>.05/mcc) = NaN;
        hh = heatmap(rr);
        hh.ColorLimits = [-1 1];
        hh.Colormap = c;
        hh.XDisplayLabels = ROIs_all;
        hh.YDisplayLabels = illMagn;
        title('Corrected significant')
        set(gca,'FontSize',fntsz)
    end
end


