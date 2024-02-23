function [block_output, trial_output] = apply_cutoff_D_performance(taskType, block_output, trial_output, plot_summary_fig)
% Applies cutoff based on performance towrad the end of discrimination phase
% Outputs erased data for subjects who did not pass cut-offs

select_option = 2;
    % Option1: From whole-block metrics, last-day performance > 54%
    % Option2: From trials metrics (running window average), t-test for 50% --- this is default

    Exc_mean = [];  % bar plots data bin
    Inc_mean = [];

    dFname = ["output/plot_data/Option"+select_option+"_cutoff_subjects.txt"];
    delete(dFname);
    diary(dFname);
    diary on;

    disp("Task selected: "+taskType);
    switch select_option
        case 1
            disp("Option1: From whole-block metrics, Discrim. last-day performance > 54%");
            allGroups = fieldnames(block_output.(taskType));
            for g = 1:numel(allGroups)
                allSubjects = fieldnames(block_output.(taskType).(allGroups{g}));
                disp(" >>> "+allGroups{g}+" (N = "+numel(allSubjects)+")");
                exclude_cnt = 0;
                include_cnt = 0;
                no_data = 0;
                for s = 1:numel(allSubjects)
                    SubData = block_output.(taskType).(allGroups{g}).(allSubjects{s});
                    if ~isfield(SubData,'pbetter')
                        disp("   "+s+">> "+allSubjects{s}+" : NO DATA");
                        no_data = no_data + 1;
                        continue;
                    end
                    DiscrimPerf = SubData.pbetter(SubData.D==1);
                    if DiscrimPerf(end)<.54
                        verdict = " (Exclude from analysis)";
                        exclude_cnt = exclude_cnt + 1;
                        Exc_mean = [Exc_mean; mean(DiscrimPerf)];
                        block_output.(taskType).(allGroups{g}).(allSubjects{s}) = struct;   % erase data
                        trial_output.(taskType).(allGroups{g}).(allSubjects{s}) = struct;
                    else
                        verdict = "";
                        include_cnt = include_cnt + 1;
                        Inc_mean = [Inc_mean; mean(DiscrimPerf)];
                    end
                    disp("   "+s+">> "+allSubjects{s}+" : "+DiscrimPerf(end)+verdict);
                end
                disp(" > Excluded: "+exclude_cnt+" / "+num2str(numel(allSubjects)-no_data)+", N = "+include_cnt+" left");
                fprintf('\n');
            end
            fprintf('\n\n');
        case 2
            disp("Option2: From trials metrics (running window average), t-test for >50%");
            allGroups = fieldnames(trial_output.(taskType));
            for g = 1:numel(allGroups)
                allSubjects = fieldnames(trial_output.(taskType).(allGroups{g}));
                disp(" >>> "+allGroups{g}+" (N = "+numel(allSubjects)+")");
                exclude_cnt = 0;
                include_cnt = 0;
                no_data = 0;
                for s = 1:numel(allSubjects)
                    SubData = trial_output.(taskType).(allGroups{g}).(allSubjects{s}).beforeRev; % last day of phase before reversal
                    if ~isfield(SubData.Run,'pbetter')
                        disp("   "+s+">> "+allSubjects{s}+" : NO DATA");
                        no_data = no_data + 1;
                        continue;
                    end
                    DiscrimPerf = SubData.Run.pbetter(SubData.D==1,:);      % choose last Discrim. day
                    if ~ttest(DiscrimPerf,0.5,'Alpha',0.05,'Tail','right')
                        verdict = " (Exclude from analysis)";
                        exclude_cnt = exclude_cnt + 1;
                        Exc_mean = [Exc_mean; mean(DiscrimPerf)];
                        block_output.(taskType).(allGroups{g}).(allSubjects{s}) = struct;   % erase data
                        trial_output.(taskType).(allGroups{g}).(allSubjects{s}) = struct;
                    else
                        verdict = "";
                        include_cnt = include_cnt + 1;
                        Inc_mean = [Inc_mean; mean(DiscrimPerf)];
                    end
                    disp("   "+s+">> "+allSubjects{s}+" : "+mean(DiscrimPerf,'omitnan')+verdict);
                end
                disp(" > Excluded: "+exclude_cnt+" / "+num2str(numel(allSubjects)-no_data)+", N = "+include_cnt+" left");
                fprintf('\n');
            end
            fprintf('\n\n');
    end

    diary off

    %% Compare excluded vs. included subjects performance (last day data)
    if plot_summary_fig    
        
        figure; clf
        bar([mean(Exc_mean,'omitnan'),mean(Inc_mean,'omitnan')],'LineStyle','none'); hold on
        errorbar(1:2,[mean(Exc_mean,'omitnan'),mean(Inc_mean,'omitnan')],[std(Exc_mean,'omitnan'),std(Inc_mean,'omitnan')],'LineStyle','none','LineWidth',1);
        xticklabels(["Excluded (n="+sum(~isnan(Exc_mean))+")","Included (n="+sum(~isnan(Inc_mean))+")"]);
        xlabel("All Subjects");
        ylabel("mean P(Correct) +/- SD");

        if strcmp(taskType,"SO")    
            title(["Stimulus-Outcome task","mean P(Correct) in last Disrcim. day"]);
        elseif strcmp(taskType,"SO")    
            title(["Action-Outcome task","mean P(Correct) in last Disrcim. day"]);
        end

        [h,pval,~,stats] = ttest2(Exc_mean,Inc_mean);
        if ~isnan(h)&&h
            text(1.5, 0.95, "∗",'HorizontalAlignment','center','FontSize',20);
        end
        text(0.05,0.95,["Two-sampel t-test: ","\it{t}\rm("+stats.df+") = "+num2str(stats.tstat,4)+", \it{p} = \rm"+num2str(pval,3)],'fontsize',12,'Units','normalized','VerticalAlignment','top');
        set(gca,'FontName','Helvetica','FontSize',14,'FontWeight','normal','LineWidth',1, 'tickdir','out','Box','off');
        
    end

end