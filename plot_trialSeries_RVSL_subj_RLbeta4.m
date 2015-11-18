function plot_trialSeries_RVSL_subj_RLbeta4(subj, plot_ppc)
% This founction plots 1st choice, switch not not, outcomes and missing data point

load '\projects\SocialInflu\sit_stan\_data\data3_129.mat'
load '\projects\SocialInflu\sit_stan\_outputs\param_beta4_alt6'

delete('trialSeries.ps');

for j = 1:length(subj)
    
    
    % --- load data ------------------------------------------
    data = squeeze(data3(:,:, subj(j) ));
    nt   = size(data,1);
    choice1  = data(:,3);
    missInd  = data(:,41);
    swch     = data(:,5);
    otcm2    = data(:,14);
    reversal = data(:,2);
    
    % --- get internal model variables -----------------------
    [~,~,~,model] = RevLearn_RLbeta4_alt6(parMat(subj(j),:), data);
    myValue = model.myValue;
    otherValue = model.otherValue;
    with    = model.with;
    agst    = model.agst;
    wgtWith = model.wgtWith;
    wgtAgst = model.wgtAgst;    
    
    %%% --- plot -----
    f(j) = figure;
    set(f(j),'color',[1 1 1], 'position', [20 200 1300 800]);
    
    % --- plot choice1's history
    plot(1:nt, data(:,3), 'k:', 'linewidth', 1)
    hold on
    
    nt_1 = find(data(:,3) == 1);
    nt_2 = find(data(:,3) == 2);
    plot(nt_1, ones(length(nt_1),1), 'ko', 'MarkerSize',5, 'MarkerFaceColor', 'm')
    plot(nt_2, 2*ones(length(nt_2),1), 'ko', 'MarkerSize',5, 'MarkerFaceColor', 'g')
        

    % --- plot missing choice
    missChoice = choice1(logical(missInd));
    missTrial  = find(missInd==1);
    plot(missTrial,missChoice,'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r')
    
    % --- plot 2nd outcome
    plot(find(otcm2==1),  otcm2(otcm2==1)*2.5, 'g.')
    plot(find(otcm2==-1), otcm2(otcm2==-1)*(-2.5), 'r.')
    
    % --- plot switch or not
    plot(find(swch==1), swch(swch==1)*3, 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b')
    
    % --- plot my Value of option1 and option2
    plot(1:nt, myValue(:,1)-0.5,'m.-')
    plot(1:nt, myValue(:,2)-0.5,'g.-')
    H1 = text(-16,-1, 'my value');
    set(H1, 'rotation',90, 'position', [-3 -1 0])
    
    % --- plot other Value of option1 and option2
    plot(1:nt, otherValue(:,1)-4.5,'m.-')
    plot(1:nt, otherValue(:,2)-4.5,'g.-')
    H2 = text(-16,-5, 'other value');
    set(H2, 'rotation',90, 'position', [-3 -5 0])
    
    % text(101, -1.8, 'option 1: magenta')
    % text(101, -2.3, 'option 2: green')
    
    
    % --- plot with /against
    plot(1:nt, with-9, 'g.-')
    plot(1:nt, agst-9, 'r.-')
    H3 = text(-15.5,-8.5, 'unweighted with /against');
    
    % --- plot weighted with / agsint
    plot(1:nt, wgtWith-10.5, 'g.-')
    plot(1:nt, wgtAgst-10.5, 'r.-')
    H4 = text(-15.5,-10, 'weighted with / against');
        
    line(xlim,[-7.5 -7.5], 'color','k','lineStyle','-')
    text(101, -9, 'with: green')
    text(101, -9.5, 'against: red')
    
    %%% plot ppc -------------------------------------------------
    if nargin > 1 && plot_ppc
        load '\projects\SocialInflu\sit_stan\_outputs\ppc_RLbeta4_a6_pointwise.mat'
        
        acc = ppc_RLbeta4_a6_pointwise(subj(j),:);
        y_pos = -13;
        plot(1:nt, acc+y_pos, 'b-.', 'linewidth', 2.5)
        
        line(xlim,[y_pos+1.5 y_pos+1.5], 'color','k','lineStyle','-')
        line(xlim,[y_pos+0.5, y_pos+0.5 ], 'color','k','lineStyle',':')        
        
        H5 = text(-16,-9, 'predictive accuracy');
        set(H5, 'rotation',90, 'position', [-5 y_pos-1.5 0])
        set(gca, 'YTick', [y_pos y_pos+1, -10.5 -9.5 -9 -8 -4.5 -0.5, 1 2 2.5 3], 'YTickLabel', ...
        {'0.00','1.00','0','1','0','1','0','0','1st choice: option1','1st choice: option2','outcome','switch'} )
        
        ylim([-14 4])
    
    elseif nargin == 1
        ylim([-8 , 3.5])
        set(gca, 'YTick', [1 2 2.5 3], 'YTickLabel', ...
        {'1st choice: option1','1st choice: option2','outcome','switch'} )
    end
    
    % --- plot rewardProb reversal
    rev_ind = find(reversal == 1);
    for r = 1:length(rev_ind)
        line([rev_ind(r), rev_ind(r)], ylim, 'color','r','lineStyle','--')                
    end      
    
    title(sprintf('subject No. %d; PARAM: lr = %3.2f, disc = %3.2f, beta1-6 = [%3.2f %3.2f %3.2f %3.2f %3.2f %3.2f]', ...
        subj(j), parMat(subj(j),:)))
    xlabel('trials')
    
    hold off
    
    % --- save plots into file
    print('-f', '-dpsc2','-append','-loose','-r150', 'trialSeries.ps')
    
end
