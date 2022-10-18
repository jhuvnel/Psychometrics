clc; clear; close all

% Author: CFB. Note -- code for plotting figures was adapted from AA's
% scripts

path1 = 'H:\MVI\DATA SUMMARY\IN PROGRESS\Psychometrics\Stimulation Pause Detection\JVR-DMD-2023\';
patients = {'MVI001R019','MVI002R004','MVI003R140','MVI004R201','MVI006R296','MVI009R908','MVI010R141'};
path2 = 'NoDistractionsANALYZED.mat';
%patients = {'MVI001','MVI002','MVI003','MVI004','MVI006','MVI009','MVI010'};

sub_mark = 'xdo^ps+hv<'; %MVI001-MVI010

means = NaN(1,length(patients));
stds = means;

MRL = NaN(2,length(patients)); % first row is for the max average response latency, second row is for the std of that max
mRL = MRL; % same as above but for the minimum RL
DRL = MRL;


% First extract data and plot tiled layout of psychometric curves
stdcolor = [192 0 0]./250;
meancolor = [192 0 0]./250;
fig1 = figure(1);
clf;
set(fig1,'Color',[1,1,1],'Position',[115   172   978   1000])
ha = gobjects(1,6);
for i = 1:6
    ha(i) = subplot(2,3,i);
end
xmin = 0.06;
xmax = 0.99;
xspac = 0.01;
ymin = 0.08;
ymax = 0.9;
yspac = 0.03;
ywid = (ymax-ymin-yspac)/2;
xwid = (xmax-xmin-2*xspac)/3;
xpos = repmat(xmin:(xwid+xspac):xmax,1,2);
ypos = reshape(repmat([(ymax-ywid),ymin],3,1),[],1)';
for i = 1:6
    ha(i).Position = [xpos(i),ypos(i),xwid,ywid];
end
set(ha,'box','on')
annotation('textbox',[0.417 .02 1 .02],'String','Stimulation Pause Duration (ms)','FontSize',14,...
    'HorizontalAlignment','left','EdgeColor','none');
annotation('textarrow',[0.02,0.02],[0.5,1],'String','Fraction Detected',...
    'TextRotation',90,'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle',...
    'LineStyle','none','HeadStyle','none');

for i = 1:length(patients)
    %load([patients{i},'.mat'])
    load([path1 patients{i}(1:6) '\' patients{i} path2])
    fun = @(x,xdata)cdf('Normal',xdata,x(1),x(2));
    x0 = [1,20,10];
    x = lsqcurvefit(fun,x0,xdata,ydata);
    times = linspace(xdata(1),xdata(end));
    [MRL(1,i),Midx] = max(detectionsummary(:,5));
    [mRL(1,i),midx] = min(detectionsummary(:,5));
    MRL(2,i) = detectionsummary(Midx,6);
    mRL(2,i) = detectionsummary(midx,6);
    RL(:,i) = detectionsummary(:,5); % RL column removing initial NaN values
    DRL(1,i) = RL(end,i)-RL(1,i);
    if i == 1
        means(i) = x(1);
        stds(i) = x(2);
    else
        axes(ha(i-1))
        means(i) = x(1);
        stds(i) = x(2);
        plot(xdata,ydata,'ko','MarkerFaceColor','k')
        hold on
        plot(times,fun(x,times),'k-','LineWidth',2)
%         plot(x(1),0.5,'x','MarkerSize',10,'LineWidth',3,'Color',[.5 .5 .5])
%         plot(x(1) + [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',[.5 .5 .5])
%         plot(x(1) - [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',[.5 .5 .5])
        plot(x(1),0.5,'x','MarkerSize',10,'LineWidth',3,'Color',meancolor)
        plot(x(1) + [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',stdcolor)
        plot(x(1) - [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',stdcolor)
        hl(i-1,1) = plot(NaN,NaN,'x','MarkerSize',10,'LineWidth',3,'Color',meancolor);
        hl(i-1,2) = plot(NaN,NaN,':','LineWidth',2,'Color',stdcolor);
        grid on
        grid minor
        %legend('Data',sprintf('Cum Gaussian Dist Fit (mean %.3g.2,sd %.3g)',x(1),x(2)));
        legend(ha(i-1),hl(i-1,:),['\mu ' sprintf('= %.3g ms',x(1))],['\sigma ' sprintf('= %.3g ms',x(2))],'Location','southeast','FontSize',12)
        % title(sprintf('Detection Rate: %s analyzed %s\nusing %s',filename,char(datetime),nameofthismfile),'FontSize',8);
        title(patients{i},'FontSize',14)
        set(gcf,'color','w');
        box on
        axis([0 1.1*max(xdata) 0 1.05]);
    end
end

set(ha(1:3),'XTickLabel',[])
set(ha([2,3,5,6]),'YTickLabel',[])
% xlabel(ha(4:6),'Stimulation Pause Duration (ms)','FontSize',14)
% ylabel(ha([1 4]),{'Fraction Detected'},'FontSize',14)

%% FIGURE 3 - Single psychometric plot with response latency 
path1 = 'H:\MVI\DATA SUMMARY\IN PROGRESS\Psychometrics\Stimulation Pause Detection\JVR-DMD-2023\';
patients = {'MVI001R019','MVI002R004','MVI003R140','MVI004R201','MVI005R107','MVI006R296','MVI007R765','MVI008R021','MVI009R908','MVI010R141'};
path2 = {'NoDistractionsANALYZED.mat', 'WithDistractionsANALYZED.mat'};

for j = 9 % which patient to plot
    figc = figure;
    set(figc,'Color',[1,1,1],'Position',[440  378  500   650])
    ha = gobjects(1,2);
    for i = 1:2
        ha(i) = subplot(1,2,i);
    end
    xmin = 0.13;
    xmax = 0.99;
    xspac = 0.02;
    ymin = 0.1;
    ymax = 0.95;
    yspac = 0.03;
    ywid = repmat((ymax-ymin-yspac)/2,1,2);
    xwid = xmax-xmin;
    ypos = [ymin ymin+ywid(1)+yspac];
    xpos = flip(reshape(repmat(xmin:(xwid+xspac):xmax,1,2),[],1)',2);
    xscale = 'log';
    load([path1 patients{j}(1:6) '\' patients{j} path2{i}])
    for i = 1:2
        ha(i).Position = [xpos(i),ypos(i),xwid,ywid(i)];
        axes(ha(i))
        grid on
        grid minor
        hold on
        ha(i).XAxis.FontSize = 16;
        ha(i).YAxis.FontSize = 16;
        if i == 2
            %legend('Data',sprintf('Cum Gaussian Dist Fit (mean %.3g.2,sd %.3g)',x(1),x(2)));
            fun = @(x,xdata)cdf('Normal',xdata,x(1),x(2));
            x0 = [1,20,10];
            x = lsqcurvefit(fun,x0,xdata,ydata);
            times = linspace(xdata(1),xdata(end));
            plot(xdata,ydata,'ko','MarkerFaceColor','none','Marker',sub_mark(j),'MarkerSize',7)
            plot(times,fun(x,times),'k-','LineWidth',2)
            hs(1) = plot(x(1),0.5,'*','MarkerSize',10,'LineWidth',1.5,'Color',meancolor);
            hs(2) = plot(x(1) + [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',stdcolor);
            plot(x(1) - [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',stdcolor)
            legend(ha(i),hs,['\mu ' sprintf('= %.3g ms',x(1))],['\sigma ' sprintf('= %.3g ms',x(2))],'Location','northwest','FontSize',12)
            title(patients{j},'FontSize',16)
            ylabel('Fraction Detected','FontSize',16);
            axis([3 1.1*300 0 1.05]);
        else
            errorbar(xdata,detectionsummary(:,5),detectionsummary(:,6),'k','Marker',sub_mark(j),'MarkerFaceColor','none','LineWidth',1)
            legend(ha(i),'\mu \pm \sigma','Location','southwest','FontSize',12)
            ylabel('Response Latency (s)','FontSize',16);
            axis([3 1.1*300 0.3 1.3]);
        end
        %     ylabel('Fraction Detected','FontSize',14);
        %     xlabel('Stimulation Pause Duration (ms)','FontSize',14);
        % title(sprintf('Detection Rate: %s analyzed %s\nusing %s',filename,char(datetime),nameofthismfile),'FontSize',8);
        set(gcf,'color','w');
        box on
        set(gca,'Xscale',xscale)
        % set(gca,'Xscale','log')
    end
end
annotation('textbox',[0.217 .035 1 .02],'String','Stimulation Pause Duration (ms)','FontSize',16,...
   'HorizontalAlignment','left','EdgeColor','none');
% annotation('textarrow',[0.015,1],[0.5,1],'String','Fraction Detected',...
%     'TextRotation',90,'FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle',...
%     'LineStyle','none','HeadStyle','none');
set(ha(2),'XTickLabels','')


%% FIGURE 4/6 -  Population cdf
path1 = 'H:\MVI\DATA SUMMARY\IN PROGRESS\Psychometrics\Stimulation Pause Detection\JVR-DMD-2023\';
path2 = {'NoDistractionsANALYZED.mat'};
patients = {'MVI001R019','MVI003R140','MVI004R201','MVI005R107','MVI006R296','MVI008R021','MVI009R908','MVI010R141'}; %MVI007 fell asleep
path2 = {'WithDistractionsANALYZED.mat'};
patients = {'MVI001R019','MVI003R140','MVI005R107','MVI006R296','MVI007R765','MVI008R021','MVI009R908','MVI010R141'}; %MVI004 button didn't work


xscale = 'log';
sub_mark = 'xdo^ps+hv<'; %MVI001-MVI010
stdcolor = [192 0 0]./250;
meancolor = [192 0 0]./250;
figc = figure;
set(figc,'Color',[1,1,1],'Position',[440  378  500   650])
ha = gobjects(1,1);
for i = 1:1
    ha(i) = subplot(1,1,i);
end
xmin = 0.13;
xmax = 0.99;
xspac = 0.02;
ymin = 0.1;
ymax = 0.95;
yspac = 0.03;
ywid = repmat((ymax-ymin-yspac)/1,1,1);
xwid = xmax-xmin;
ypos = [ymin ymin+ywid(1)+yspac];
xpos = flip(reshape(repmat(xmin:(xwid+xspac):xmax,1,1),[],1)',1);
xscale = 'log';
for i = 1:1
    ha(i).Position = [xpos(i),ypos(i),xwid,ywid(i)];
end

yspacBIG = 0.03;
xspacBIG = 0.02;

set(ha,'box','on')
xdataALL = [];
xdataPOOLED = [0 4 9 19 29 39 59 79 99 198 298]; % same as xdataALL but binning similar pauses like 297 ms and 298 ms for easier calculations
markALL = [];
ydataALL = [];
ydataPOOLED = zeros(length(patients),length(xdataPOOLED));
repdataPOOLED = ydataPOOLED; % number of reps per duration and per subject
if strcmp(xscale,'linear')
    jitf = 0.05;
%     jitf = 0;
elseif strcmp(xscale,'log')
    jitf = [0 0.4 1.2 1.5 1.3 3 3 3 5 12 15]/50;
%     jitf = 0;
end
cALL = [];
order = [5:8 4:-1:2];
RL = NaN(length(patients),length(xdataPOOLED));
for i = 1:length(patients)+1
    if i <= length(patients)
        sub_id = str2num(patients{i}(4:6));
        load([path1 patients{i}(1:6) '\' patients{i} path2{1}])
        xdataALL = [xdataALL; xdata];
        ydataALL = [ydataALL; ydata];
        ydataPOOLED(i,:) = detectionsummary(:,3)';
        repdataPOOLED(i,:) = detectionsummary(:,2)';
        markALL = [markALL repmat(sub_mark(sub_id),1,length(xdata))];
        fun = @(x,xdata)cdf('Normal',xdata,x(1),x(2));
        x0 = [1,20,10];
        x = lsqcurvefit(fun,x0,xdata,ydata);
        times = linspace(xdata(1),xdata(end));
%         [MRL(1,i),Midx] = max(detectionsummary(:,5));
%         [mRL(1,i),midx] = min(detectionsummary(:,5));
%         MRL(2,i) = detectionsummary(Midx,6);
%         mRL(2,i) = detectionsummary(midx,6);
        RL(i,:) = detectionsummary(:,5)'; % RL column removing initial NaN values
        omitnanRL = RL(i,~isnan(RL(i,:)));
        MRL(1,i) = omitnanRL(end); % RL associated with longest duration stim pause detected
        mRL(1,i) = omitnanRL(1); % RL associated with shortest duration stim pause detected
        DRL(1,i) = omitnanRL(end) - omitnanRL(1);
        means(i) = x(1);
        stds(i) = x(2);
        cALL(i,:) = fun(x,times);
    else
        fun = @(x,xdataPOOLED)cdf('Normal',xdataPOOLED,x(1),x(2));
        x0 = [1,20,10];
        x = lsqcurvefit(fun,x0,xdataPOOLED,sum(ydataPOOLED,1)./sum(repdataPOOLED,1)); %%%%%%%%%%%%%%%% this way of computing mean? and std?
        xp = lsqcurvefit(fun,x0,xdataPOOLED,sum(ydataPOOLED,1)./sum(repdataPOOLED,1) + std(ydataPOOLED./repdataPOOLED,1));
        xm = lsqcurvefit(fun,x0,xdataPOOLED,sum(ydataPOOLED,1)./sum(repdataPOOLED,1) - std(ydataPOOLED./repdataPOOLED,1));
%         times = linspace(xdataPOOLED(1),xdataPOOLED(end));
        times = xdataPOOLED(1):0.01:xdataPOOLED(end);
        meancdf = fun(x,times);
        meancdf_pts = meancdf(100*xdataPOOLED+1);
        prms = sum((ydataPOOLED./repdataPOOLED - meancdf_pts).^2,1)/length(xdataPOOLED);
%         xp =meancdf_pts+sum(prms);
%         xm = meancdf_pts-sum(prms);
        axes(ha(1))
        means(i) = x(1);
        stds(i) = x(2);
        for j = 1:size(ydataPOOLED,1)
            sub_id = str2num(patients{j}(4:6));
            if j == 1
                plot(xdataPOOLED+(jitf)*(1 + 10*abs(j-3))*sign(j-3),ydataPOOLED(j,:)./repdataPOOLED(j,:),'k','Marker',sub_mark(sub_id),'MarkerFaceColor','k','LineWidth',1,'LineStyle','none')
            else
                plot(xdataPOOLED+(jitf)*(1 + 10*abs(j-3))*sign(j-3),ydataPOOLED(j,:)./repdataPOOLED(j,:),'k','Marker',sub_mark(sub_id),'MarkerFaceColor','none','LineStyle','none')
            end
            hold on
        end
        set(gca,'Xscale',xscale)
        plot(times,fun(x,times),'k-','LineWidth',2)
        times(1) = 0.001;
        fill([times fliplr(times)], [fun(xp,times) fliplr(fun(xm,times))],[.5 .5 .5],'facealpha',.25);
        xdataPOOLED(1) = xdataPOOLED(1) + 0.001;
%         fill([xdataPOOLED fliplr(xdataPOOLED) xdataPOOLED(1)], [xp fliplr(xm) xp(1)],[.5 .5 .5],'facealpha',.25);
        %         plot(x(1),0.5,'x','MarkerSize',10,'LineWidth',3,'Color',[.5 .5 .5])
        %         plot(x(1) + [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',[.5 .5 .5])
        %         plot(x(1) - [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',[.5 .5 .5])
        plot(x(1),0.5,'*','MarkerSize',10,'LineWidth',1.5,'Color',meancolor)
        plot(x(1) + [x(2) x(2)]/2,[-10 1.5],':','LineWidth',2,'Color',stdcolor)
        plot(x(1) - [x(2) x(2)]/2,[-10 1.5],':','LineWidth',2,'Color',stdcolor)
        hl(1,1) = plot(NaN,NaN,'*','MarkerSize',10,'LineWidth',1.5,'Color',meancolor);
        hl(1,2) = plot(NaN,NaN,':','LineWidth',2,'Color',stdcolor);
        hl2(1) = plot(NaN,NaN,'LineWidth',2,'Color','k');
        hl2(2) = fill(NaN,NaN,[.5 .5 .5],'facealpha',.25);
        grid on
        grid minor
        %legend('Data',sprintf('Cum Gaussian Dist Fit (mean %.3g.2,sd %.3g)',x(1),x(2)));
        lps = legend(ha(1),[hl(1,:) hl2],['\mu ' sprintf('= %.3g ms',x(1))],['\sigma ' sprintf('= %.3g ms',x(2))],'CDF fit to group \mu','CDF fit to group \sigma','Location','northwest','FontSize',12);
        % title(sprintf('Detection Rate: %s analyzed %s\nusing %s',filename,char(datetime),nameofthismfile),'FontSize',8);
        if strcmp(xscale,'log')
            lps.Location = 'northwest';
%             lps.Position = [0.1495 0.855 0.05 0.05];
            lps.FontSize = 12.05;
        elseif strcmp(xscale,'linear')
            lps.Location = 'southeast';
        end
        title('Pooled Data for All Participants','FontSize',14)
%         title('Linear Scale','FontSize',14)
%         xlabel('Stimulation Pause Duration (ms)','FontSize',14)
%         ylabel('Fraction Detected','FontSize',14)
        set(gcf,'color','w');
        box on
        axis([3 1.1*max(xdataALL) -0.05 1.05]);
    end
end

xlabel('Stimulation Pause Duration (ms)','FontSize',14)
ylabel({'Fraction Detected'},'FontSize',14)
% 
% annotation('textbox',[0.417 .02 1 .02],'String','Stimulation Pause Duration (ms)','FontSize',14,...
%     'HorizontalAlignment','left','EdgeColor','none');
% annotation('textarrow',[0.02,0.02],[0.5,1],'String','Fraction Detected',...
%     'TextRotation',90,'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle',...
%     'LineStyle','none','HeadStyle','none');
% For insert figure
% axes(ha(1))
% title('Linear Scale','FontSize',14)
% xlabel('Stimulation Pause Duration (ms)','FontSize',14)
% ylabel('Fraction Detected','FontSize',14)
% axes(ha(5))
% title('')
% axes(ha(6))
% title('')
% axes(ha(7))
% title('')


%% FIGURE 5/7 - Population RL's
path1 = 'H:\MVI\DATA SUMMARY\IN PROGRESS\Psychometrics\Stimulation Pause Detection\JVR-DMD-2023\';
path2 = {'NoDistractionsANALYZED.mat'};
patients = {'MVI001R019','MVI003R140','MVI004R201','MVI005R107','MVI006R296','MVI008R021','MVI009R908','MVI010R141'}; %MVI007 fell asleep
path2 = {'WithDistractionsANALYZED.mat'};
patients = {'MVI001R019','MVI003R140','MVI005R107','MVI006R296','MVI007R765','MVI008R021','MVI009R908','MVI010R141'}; %MVI004 button didn't work

xscale = 'log';
sub_mark = 'xdo^ps+hv<'; %MVI001-MVI010
stdcolor = [192 0 0]./250;
meancolor = [192 0 0]./250;
figc = figure;
set(figc,'Color',[1,1,1],'Position',[440  378  500   650])
ha = gobjects(1,1);
for i = 1:1
    ha(i) = subplot(1,1,i);
end
xmin = 0.13;
xmax = 0.99;
xspac = 0.02;
ymin = 0.1;
ymax = 0.95;
yspac = 0.03;
ywid = repmat((ymax-ymin-yspac)/1,1,1);
xwid = xmax-xmin;
ypos = [ymin ymin+ywid(1)+yspac];
xpos = flip(reshape(repmat(xmin:(xwid+xspac):xmax,1,1),[],1)',1);
xscale = 'log';
for i = 1:1
    ha(i).Position = [xpos(i),ypos(i),xwid,ywid(i)];
end

yspacBIG = 0.03;
xspacBIG = 0.02;

set(ha,'box','on')
xdataALL = [];
xdataPOOLED = [0 4 9 19 29 39 59 79 99 198 298]; % same as xdataALL but binning similar pauses like 297 ms and 298 ms for easier calculations
rlPOOLED = zeros(length(patients),length(xdataPOOLED)); % response latencies
if strcmp(xscale,'linear')
    jitf = 0.05;
%     jitf = 0;
elseif strcmp(xscale,'log')
    jitf = [0 0.4 1.2 1.5 1.3 3 3 3 5 12 15]/50;
%     jitf = 0;
end
cALL = [];
order = [5:8 4:-1:2];
RL = NaN(length(patients),length(xdataPOOLED));
for i = 1:length(patients)+1
    if i <= length(patients)
        sub_id = str2num(patients{i}(4:6));
        load([path1 patients{i}(1:6) '\' patients{i} path2{1}])
        xdataALL = [xdataALL; xdata];
        rlPOOLED(i,:) = detectionsummary(:,5)';
    else
        times = linspace(xdataPOOLED(1),xdataPOOLED(end));
        axes(ha(1))
        means(i) = x(1);
        stds(i) = x(2);
        for j = 1:size(rlPOOLED,1)
            sub_id = str2num(patients{j}(4:6));
            if j == 1
                plot(xdataPOOLED + jitf,rlPOOLED(j,:),'k','Marker',sub_mark(sub_id),'MarkerFaceColor','k','LineWidth',0.7,'LineStyle','none')
            else
                plot(xdataPOOLED + jitf,rlPOOLED(j,:),'k','Marker',sub_mark(sub_id),'MarkerFaceColor','none','LineStyle','none')
            end
            hold on
        end
%         errorbar(xdataPOOLED,mean(rlPOOLED,1,'omitnan'),std(rlPOOLED,1,'omitnan'),'k','Marker','.','MarkerFaceColor','k','LineWidth',2);
        plot(xdataPOOLED,mean(rlPOOLED,1,'omitnan'),'k','Marker','.','MarkerFaceColor','k','LineWidth',2);
        xdataPOOLED(1) = 0.01;
        upp = mean(rlPOOLED,1,'omitnan')+std(rlPOOLED,1,'omitnan')./sum(~isnan(rlPOOLED));
        loww = fliplr(mean(rlPOOLED,1,'omitnan')-std(rlPOOLED,1,'omitnan')./sum(~isnan(rlPOOLED)));
%         fill([xdataPOOLED(1:5) fliplr(xdataPOOLED(1:5))], [upp(1:5) loww(end-4:end)],[.5 .5 .5],'facealpha',.1,'LineStyle','none'); % for quiet condition
        fill([xdataPOOLED(2:5) fliplr(xdataPOOLED(2:5))], [upp(2:5) loww(end-4:end-1)],[.5 .5 .5],'facealpha',.1,'LineStyle','none'); % fror distraction condition, exclude NaN
        fill([xdataPOOLED(5:end) fliplr(xdataPOOLED(5:end))], [upp(5:end) loww(1:end-4)],[.5 .5 .5],'facealpha',.4,'LineStyle','none');
        fill([xdataPOOLED(2:end) fliplr(xdataPOOLED(2:end))], [upp(2:end) loww(1:end-1)],[.5 .5 .5],'facealpha',0);
        plot([0.01 400],ones(1,2),'k:','LineWidth',1.5)
        %         plot(x(1),0.5,'x','MarkerSize',10,'LineWidth',3,'Color',[.5 .5 .5])
        %         plot(x(1) + [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',[.5 .5 .5])
        %         plot(x(1) - [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',[.5 .5 .5])
        hl(1,1) = errorbar(NaN,NaN,'k','Marker','.','MarkerFaceColor','k','LineWidth',2);
        hl(1,2) = plot(NaN,NaN,'k:','LineWidth',1.5);
        hl2(1) = fill(NaN,NaN,[.5 .5 .5],'facealpha',.25);
        grid on
        grid minor
        %legend('Data',sprintf('Cum Gaussian Dist Fit (mean %.3g.2,sd %.3g)',x(1),x(2)));
        legend(ha(1),[hl(1,1) hl2(1) hl(1,2)],'Group \mu','\pmSEM','Latency = 1 s','Location','southwest','FontSize',12)
        % title(sprintf('Detection Rate: %s analyzed %s\nusing %s',filename,char(datetime),nameofthismfile),'FontSize',8);
        title('Pooled Data for All Participants','FontSize',14)
        set(gcf,'color','w');
        box on
        set(gca,'Xscale',xscale)
        if strcmp(xscale,'linear')
            axis([0 1.1*max(xdataPOOLED) 0.35 1.1]);
        elseif strcmp(xscale,'log')
            axis([3.5 1.1*max(xdataPOOLED) 0.35 1.1]);
        end
    end
end

xlabel('Stimulation Pause Duration (ms)','FontSize',14)
ylabel({'Response Latency (s)'},'FontSize',14)

% annotation('textbox',[0.417 .02 1 .02],'String','Stimulation Pause Duration (ms)','FontSize',14,...
%     'HorizontalAlignment','left','EdgeColor','none');
% annotation('textarrow',[0.02,0.02],[0.5,1],'String','Fraction Detected',...
%     'TextRotation',90,'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle',...
%     'LineStyle','none','HeadStyle','none');
% For insert figure
% axes(ha(1))
% title('Linear Scale','FontSize',14)
% xlabel('Stimulation Pause Duration (ms)','FontSize',14)
% ylabel('Fraction Detected','FontSize',14)
% axes(ha(5))
% title('')
% axes(ha(6))
% title('')
% axes(ha(7))
% title('')




%% Side-by-side comparison of 2 psychometric plots
path1 = 'H:\MVI\DATA SUMMARY\IN PROGRESS\Psychometrics\Stimulation Pause Detection\JVR-DMD-2023\';
patients = {'MVI001R019','MVI002R004','MVI003R140','MVI004R201','MVI005R107','MVI006R296','MVI007R765','MVI008R021','MVI009R908','MVI010R141'};
path2 = {'NoDistractionsANALYZED.mat', 'WithDistractionsANALYZED.mat'};

for j = 1:10 % which patient to plot
    figc = figure;
    set(figc,'Color',[1,1,1],'Position',[440  378  1200   650])
    ha = gobjects(1,2);
    for i = 1:2
        ha(i) = subplot(1,2,i);
    end
    xmin = 0.07;
    xmax = 0.99;
    xspac = 0.02;
    ymin = 0.1;
    ymax = 0.95;
    yspac = 0.03;
    ywid = (ymax-ymin);
    xwid = repmat((xmax-xmin-xspac)/2,1,2);
    xpos = [xmin xmin+xwid(1)+xspac];
    ypos = flip(reshape(repmat(ymin:(ywid+yspac):ymax,1,2),[],1)',2);
    xscale = 'log';
    for i = 1:2
        ha(i).Position = [xpos(i),ypos(i),xwid(i),ywid];
        axes(ha(i))
        if i == 1
            load([path1 patients{j}(1:6) '\' patients{j} path2{i}])
            title('Quiet','FontSize',16)
        else
            load([path1 patients{j}(1:6) '\' patients{j} path2{i}])
            title('Distractions','FontSize',16)
        end
        hold on
        fun = @(x,xdata)cdf('Normal',xdata,x(1),x(2));
        x0 = [1,20,10];
        x = lsqcurvefit(fun,x0,xdata,ydata);
        times = linspace(xdata(1),xdata(end));
        plot(xdata,ydata,'ko','MarkerFaceColor','none','Marker',sub_mark(j),'MarkerSize',7)
        plot(times,fun(x,times),'k-','LineWidth',2)
        hs(1) = plot(x(1),0.5,'*','MarkerSize',10,'LineWidth',1.5,'Color',meancolor);
        hs(2) = plot(x(1) + [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',stdcolor);
        plot(x(1) - [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',stdcolor)
        grid on
        grid minor
        ha(i).XAxis.FontSize = 16;
        ha(i).YAxis.FontSize = 16;
        if i == 1
        %legend('Data',sprintf('Cum Gaussian Dist Fit (mean %.3g.2,sd %.3g)',x(1),x(2)));
            legend(ha(i),hs,['\mu ' sprintf('= %.3g ms',x(1))],['\sigma ' sprintf('= %.3g ms',x(2))],'Location','northwest','FontSize',12)
        else
            legend(ha(i),'Data','Fit',['\mu ' sprintf('= %.3g ms',x(1))],['\sigma ' sprintf('= %.3g ms',x(2))],'Location','northwest','FontSize',12)
        end
        %     ylabel('Fraction Detected','FontSize',14);
        %     xlabel('Stimulation Pause Duration (ms)','FontSize',14);
        % title(sprintf('Detection Rate: %s analyzed %s\nusing %s',filename,char(datetime),nameofthismfile),'FontSize',8);
        set(gcf,'color','w');
        box on
        set(gca,'Xscale',xscale)
        axis([0 1.1*300 0 1.05]);
        % set(gca,'Xscale','log')
    end
end
annotation('textbox',[0.417 .035 1 .02],'String','Stimulation Pause Duration (ms)','FontSize',16,...
    'HorizontalAlignment','left','EdgeColor','none');
annotation('textarrow',[0.015,0.02],[0.5,1],'String','Fraction Detected',...
    'TextRotation',90,'FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle',...
    'LineStyle','none','HeadStyle','none');
set(ha(2),'YTickLabels','')


%% Side-by-side comparison of 2 response latency plots
path1 = 'H:\MVI\DATA SUMMARY\IN PROGRESS\Psychometrics\Stimulation Pause Detection\JVR-DMD-2023\';
patients = {'MVI001R019','MVI002R004','MVI003R140','MVI004R201','MVI005R107','MVI006R296','MVI007R765','MVI008R021','MVI009R908','MVI010R141'};
path2 = {'NoDistractionsANALYZED.mat', 'WithDistractionsANALYZED.mat'};
figc = figure;
set(figc,'Color',[1,1,1],'Position',[440  378  1200   650])
ha = gobjects(1,2);
for i = 1:2
    ha(i) = subplot(1,2,i);
end
xmin = 0.07;
xmax = 0.99;
xspac = 0.02;
ymin = 0.1;
ymax = 0.95;
yspac = 0.03;
ywid = (ymax-ymin);
xwid = repmat((xmax-xmin-xspac)/2,1,2);
xpos = [xmin xmin+xwid(1)+xspac];
ypos = flip(reshape(repmat(ymin:(ywid+yspac):ymax,1,2),[],1)',2);
xscale = 'log';

for j = 9 % which patient to plot
    for i = 1:2
        ha(i).Position = [xpos(i),ypos(i),xwid(i),ywid];
        axes(ha(i))
        if i == 1
            load([path1 patients{j}(1:6) '\' patients{j} path2{i}])
            title('Quiet','FontSize',16)
        else
            load([path1 patients{j}(1:6) '\' patients{j} path2{i}])
            title('Distractions','FontSize',16)
        end
        hold on
        plot(xdata,detectionsummary(:,5),'k','Marker',sub_mark(j),'MarkerFaceColor','none','LineWidth',1.5)   
        plot(xdata,detectionsummary(:,6),'k--','Marker',sub_mark(j),'MarkerFaceColor','none','LineWidth',1.5)   
        grid on
        grid minor
        ha(i).XAxis.FontSize = 16;
        ha(i).YAxis.FontSize = 16;
        if i == 2
            leg = legend(ha(i),'Mean','Std','Location','northeast','FontSize',12);
            leg.ItemTokenSize(1) = 45;
        end
        %     ylabel('Fraction Detected','FontSize',14);
        %     xlabel('Stimulation Pause Duration (ms)','FontSize',14);
        % title(sprintf('Detection Rate: %s analyzed %s\nusing %s',filename,char(datetime),nameofthismfile),'FontSize',8);
        set(gcf,'color','w');
        box on
        set(gca,'Xscale',xscale)
        axis([0 1.1*300 0 1.3]);
        % set(gca,'Xscale','log')
    end
end
annotation('textbox',[0.417 .035 1 .02],'String','Stimulation Pause Duration (ms)','FontSize',16,...
    'HorizontalAlignment','left','EdgeColor','none');
annotation('textarrow',[0.015,0.02],[0.5,1],'String','Response Latency Parameters (s)',...
    'TextRotation',90,'FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle',...
    'LineStyle','none','HeadStyle','none');
set(ha(2),'YTickLabels','')

%% Boxplots comparing thresholds quiet and distractions psychometric plots


% First extract data
for j = 1:2
    path1 = 'H:\MVI\DATA SUMMARY\IN PROGRESS\Psychometrics\Stimulation Pause Detection\JVR-DMD-2023\';
    patients = {'MVI001R019','MVI003R140','MVI005R107','MVI006R296','MVI008R021','MVI009R908','MVI010R141'};
    path2 = {'NoDistractionsANALYZED.mat','WithDistractionsANALYZED.mat'};
    for i = 1:length(patients)
        %load([patients{i},'.mat'])
        load([path1 patients{i}(1:6) '\' patients{i} path2{j}])
        fun = @(x,xdata)cdf('Normal',xdata,x(1),x(2));
        x0 = [1,20,10];
        x = lsqcurvefit(fun,x0,xdata,ydata);
        times = linspace(xdata(1),xdata(end));
        if j == 1
            means_NODis(i) = x(1);
            stds_NODis(i) = x(2);
        else
            means_WITHDis(i) = x(1);
            stds_WITHDis(i) = x(2);
        end
    end
end

constcolor = [0.768 0 0];
inccolor = [0 0.768 0];
sub_col = [constcolor constcolor inccolor constcolor inccolor inccolor inccolor];

for i = 1:length(patients)
    leg1_plist{i} = num2str(str2num(patients{i}(4:6)));
    leg1_mlist(i) = sub_mark(str2num(patients{i}(4:6)));
end

line_norm = 0.5;
line_bold = 2;
line_med = 1.5;
mark_size_big = 16;
mark_size_med = 5;
fig4 = figure(4);
set(fig4,'Color',[1,1,1],'Position',[115   172   678   1000])
ha = gobjects(1,3);
for i = 1:2
   ha(i) = subplot(1,2,i); 
end
xmin = 0.1;
xmax = 0.95;
xspac = 0.01;
xspacBIG = 0.01;%0.07;
ymin = [0.08 0.08];
ymax = [0.95 0.83];
yspac = 0.03;
ywid = (ymax-ymin);
xwid = repmat((xmax-xmin-xspacBIG)/1.6,1,2);
xwid(1) = xwid(1)*1;
xwid(end) = xwid(end)*0.5;
xpos = [xmin xmin+xwid(1)+xspacBIG];
ypos = ymin;

axes(ha(1))
offs2 = 0.23*rand(length(patients),1)+1.15;
offs3 = 0.5+ offs2;
plot(NaN,NaN)
hold on
boxp = boxplot([means_NODis(1,:)' means_WITHDis(1,:)'],'Widths',0.1,'Colors','k','Symbol','','whisker',1000);
set(boxp,'LineWidth',line_med)
for i = 1:length(patients)
    plot(offs2(i),means_NODis(1,i),'Marker',sub_mark(str2num(patients{i}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor',sub_col(1+3*(i-1):3+3*(i-1)),'LineWidth',line_med)
    plot(offs3(i),means_WITHDis(1,i),'Marker',sub_mark(str2num(patients{i}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor',sub_col(1+3*(i-1):3+3*(i-1)),'LineWidth',line_med)
    plot([offs2(i), offs3(i)],[means_NODis(1,i) means_WITHDis(1,i)],'k','LineWidth',1,'LineStyle',':')
    %     errorbar(offs2(i),means(i),stds(i),'Marker','o','MarkerSize',7,'MarkerFaceColor',meancolor,'MarkerEdgeColor',meancolor,'Color',stdcolor,'CapSize',15,'LineWidth',2);
%     plot(4+offs2(i),diff_params(i,1),'k:','LineWidth',line_norm);
    h1(i) = plot(NaN,NaN,'Marker',sub_mark(str2num(patients{i}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',line_med,'LineStyle','none');
end
plot([0 10], ones(1,2),'k:','LineWidth',1.5)
h2 = plot(NaN,NaN,'k:','LineWidth',1.5);
% %%%plot([offs2(1), 1+offs2(1)],[mRL(1,1) MRL(1,1)],'Marker',sub_mark(str2num(patients{1}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',2,'LineStyle','none')
% plot([offs2(1), 1+offs2(1)],[means_NODis(1,1) means_WITHDis(1,1)],'k','LineWidth',1,'LineStyle',':')
%plot((4:5)+offs2(i),params(:,4:5),'k-','LineWidth',line_norm);
% plot mean and std
hold off


axes(ha(2))
offs2 = 0.35*rand(length(patients),1)+1.1;
plot(NaN,NaN)
hold on
boxp = boxplot(means_WITHDis(1,:) - means_NODis(1,:),'Widths',0.1,'Colors','k','Symbol','','whisker',1000);
set(boxp,'LineWidth',line_med)
for i = 1:length(patients)
    plot(offs2(i),means_WITHDis(1,i) - means_NODis(1,i),'Marker',sub_mark(str2num(patients{i}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor',sub_col(1+3*(i-1):3+3*(i-1)),'LineWidth',line_med)
    %     errorbar(offs2(i),means(i),stds(i),'Marker','o','MarkerSize',7,'MarkerFaceColor',meancolor,'MarkerEdgeColor',meancolor,'Color',stdcolor,'CapSize',15,'LineWidth',2);
%     plot(4+offs2(i),diff_params(i,1),'k:','LineWidth',line_norm);
end
plot([0 10], zeros(1,2),'k--','LineWidth',2)
h3 = plot(NaN,NaN,'k--','LineWidth',2);
% %%%plot(offs2(1),MRL(1,1) - mRL(1,1),'Marker',sub_mark(str2num(patients{1}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',2,'LineStyle','none')
%plot((4:5)+offs2(i),params(:,4:5),'k-','LineWidth',line_norm);
% plot mean and std
hold off

annotation('textarrow',[0.02,0.02],[0.5,1],'String','Time (ms)',...
    'TextRotation',90,'FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle',...
    'LineStyle','none','HeadStyle','none');
annotation('textarrow',[xmin+xwid(1)+xwid(2)+0.08,0.02],[0.5,1],'String','\DeltaTime (ms)',...
    'TextRotation',90,'FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle',...
    'LineStyle','none','HeadStyle','none');

annotation('textbox',[0.087 .9613 1 .02],'String','Detection Thresholds for Stimulation Pause Durations','FontSize',16,...
    'HorizontalAlignment','left','EdgeColor','none','FontWeight','bold');

set(ha(1),'XTickLabel','Threshold','FontSize',16)
% xlabel(ha(2),'Walking Speed (mph)')
set(ha(1:2),'Xgrid','on','YGrid','on')
set(ha(1),'XLim',[0.8 2.2],'YLim',[10 90])
set(ha(2),'XLim',[0.8 1.5],'YLim',[-5 35])
set(ha(1),'XTick',1:2,'XTickLabel',{'Quiet','Distraction'},'XTickLabelRotation',0,'box','on','FontSize',16)
set(ha(2),'XTick',1:5,'XTickLabel',{'Change'},'XTickLabelRotation',0,'box','on','FontSize',16)
% set(ha(3),'YAxisLocation','right')
% title(ha(2),'Response Latencies corresponding to','FontSize',14)
% title(ha(3),'Shortest and Longest Detected Stimulation Pauses','FontSize',14)
% leg3_reord = [1,6,2,7,3,8,4,9,5,10];
leg1 = legend(ha(1),h1,leg1_plist,'NumColumns',4,'Location',[0.635201584432111,0.8595,0.262386430678466,0.0905],'FontSize',16);
leg1.ItemTokenSize(1) = 13;
title(leg1,'Subjects')
%leg2 = legend(ha(2),h2,'Latency = 1s','NumColumns',1,'Location','north');
%leg2.ItemTokenSize(1) = 45;
leg2 = legend(ha(2),h2,' No difference','NumColumns',1,'Location',[0.635201584432111,0.828333333333334,0.262386430678466,0.031333333333333],'FontSize',16);
leg2.ItemTokenSize(1) = 22;
% set(ha(3),'XLim',[3.3 4.5],'YLim',YLim3)
set(ha(2),'XTick',1,'XTickLabel',{'Difference'},'XTickLabelRotation',0,'box','on')
set(ha(2),'YAxisLocation','right')
% title(ha(3),'D. \DeltaDVAS Relative to Tonic')
%ylabel(ha(2),'\DeltaTime (s)')
% leg5.Position = [0.6966    0.8737    0.2040    0.0820];
ha(1).XAxis.FontSize = 16;
ha(2).XAxis.FontSize = 16;

for i = 1:2
    ha(i).Position = [xpos(i),ypos(i),xwid(i),ywid(i)];
end


% set(ha(3),'Position',[xpos(2),ypos(1),xwid2,ymax-ymin])
% set(ha(3),'Position',[xpos(3),ypos(1),xwid2,ymax-ymin])
% set(ha(1:2),'XTick',0:0.5:3,'XLim',XLim,'YLim',YLim,'XGrid','on','YGrid','on','box','on')

%% Boxplots comparing thresholds quiet and distractions response latencies


% First extract data
for j = 1:2
    path1 = 'H:\MVI\DATA SUMMARY\IN PROGRESS\Psychometrics\Stimulation Pause Detection\JVR-DMD-2023\';
    patients = {'MVI001R019','MVI003R140','MVI005R107','MVI006R296','MVI008R021','MVI009R908','MVI010R141'};
    path2 = {'NoDistractionsANALYZED.mat','WithDistractionsANALYZED.mat'};
    for i = 1:length(patients)
        %load([patients{i},'.mat'])
        load([path1 patients{i}(1:6) '\' patients{i} path2{j}])
        fun = @(x,xdata)cdf('Normal',xdata,x(1),x(2));
        x0 = [1,20,10];
        x = lsqcurvefit(fun,x0,xdata,ydata);
        times = linspace(xdata(1),xdata(end));
        if j == 1
            means_NODis(i) = x(1);
            stds_NODis(i) = x(2);
            RLnan = detectionsummary(~isnan(detectionsummary(:,5)),5:6);
            MRL_NODis(1,i) = RLnan(end,1); %corresponding to longest SP
            mRL_NODis(1,i) = RLnan(1,1); %corresponding to shortest SP
            MRL_NODis(2,i) = RLnan(end,2);
            mRL_NODis(2,i) = RLnan(1,2);
            DRL_NODis(1,i) = MRL_NODis(1,i)-mRL_NODis(1,i);
        else
            means_WITHDis(i) = x(1);
            stds_WITHDis(i) = x(2);
            RLnan = detectionsummary(~isnan(detectionsummary(:,5)),5:6);
            MRL_WITHDis(1,i) = RLnan(end,1); %corresponding to longest SP
            mRL_WITHDis(1,i) = RLnan(1,1); %corresponding to shortest SP
            MRL_WITHDis(2,i) = RLnan(end,2);
            mRL_WITHDis(2,i) = RLnan(1,2);
            DRL_WITHDis(1,i) = MRL_WITHDis(1,i)-mRL_WITHDis(1,i);
        end
    end
end

constcolor = [0.768 0 0];
inccolor = [0 0.768 0];
sub_col = [constcolor constcolor inccolor constcolor inccolor inccolor inccolor];

for i = 1:length(patients)
    leg1_plist{i} = num2str(str2num(patients{i}(4:6)));
    leg1_mlist(i) = sub_mark(str2num(patients{i}(4:6)));
end

line_norm = 0.5;
line_bold = 2;
line_med = 1.5;
mark_size_big = 16;
mark_size_med = 5;
fig4 = figure(4);
set(fig4,'Color',[1,1,1],'Position',[115   172   978   1000])
ha = gobjects(1,3);
ha = gobjects(1,3);
for i = 1:3
   ha(i) = subplot(1,3,i); 
end
xmin = 0.07;
xmax = 0.99;
xspac = 0.01;
xspacBIG = 0.09;%0.07;
ymin = 0.08;
ymax = 0.95;
yspac = 0.03;
ywid = (ymax-ymin);
xwid = repmat((xmax-xmin-xspacBIG -xspacBIG)/3,1,3);
xwid(end) = xwid(end)*1;
xpos = [xmin xmin+xwid(1)+xspacBIG xmin+2*xwid(1)+xspacBIG+xspacBIG];
ypos = flip(reshape(repmat(ymin:(ywid+yspac):ymax,1,3),[],1)',2);

axes(ha(1))
offs2 = 0.4*rand(length(patients),1)+1.2;
plot(NaN,NaN)
hold on
boxp = boxplot([mRL_NODis(1,:)' MRL_NODis(1,:)'],'Colors','k','Symbol','','whisker',1000);
set(boxp,'LineWidth',line_med)
for i = 1:length(patients)
    plot(offs2(i),mRL_NODis(1,i),'Marker',sub_mark(str2num(patients{i}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',line_med)
    plot(1+offs2(i),MRL_NODis(1,i),'Marker',sub_mark(str2num(patients{i}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',line_med)
    plot([offs2(i), 1+offs2(i)],[mRL_NODis(1,i) MRL_NODis(1,i)],'k','LineWidth',1,'LineStyle',':')
    %     errorbar(offs2(i),means(i),stds(i),'Marker','o','MarkerSize',7,'MarkerFaceColor',meancolor,'MarkerEdgeColor',meancolor,'Color',stdcolor,'CapSize',15,'LineWidth',2);
%     plot(4+offs2(i),diff_params(i,1),'k:','LineWidth',line_norm);
    h1(i) = plot(NaN,NaN,sub_mark(str2num(patients{i}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',line_med);
end
plot([0 10], ones(1,2),'k:','LineWidth',1.5)
% %%%plot([offs2(1), 1+offs2(1)],[mRL(1,1) MRL(1,1)],'Marker',sub_mark(str2num(patients{1}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',2,'LineStyle','none')
plot([offs2(1), 1+offs2(1)],[mRL_NODis(1,1) MRL_NODis(1,1)],'k','LineWidth',1,'LineStyle',':')
%plot((4:5)+offs2(i),params(:,4:5),'k-','LineWidth',line_norm);
% plot mean and std
title('Quiet','FontSize',14)
hold off


axes(ha(2))
offs2 = 0.35*rand(length(patients),1)+1.1;
plot(NaN,NaN)
hold on
boxp = boxplot([DRL_NODis(1,:)' DRL_WITHDis(1,:)'],'Colors','k','Symbol','','whisker',1000);
set(boxp,'LineWidth',line_med)
for i = 1:length(patients)
    plot(offs2(i),MRL_NODis(1,i) - mRL_NODis(1,i),'Marker',sub_mark(str2num(patients{i}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',line_med)
    plot(1+offs2(i),MRL_WITHDis(1,i) - mRL_WITHDis(1,i),'Marker',sub_mark(str2num(patients{i}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',line_med)
    plot([offs2(i), 1+offs2(i)],[MRL_NODis(1,i) - mRL_NODis(1,i), MRL_WITHDis(1,i) - mRL_WITHDis(1,i)],'k','LineWidth',1,'LineStyle',':')
    %     errorbar(offs2(i),means(i),stds(i),'Marker','o','MarkerSize',7,'MarkerFaceColor',meancolor,'MarkerEdgeColor',meancolor,'Color',stdcolor,'CapSize',15,'LineWidth',2);
%     plot(4+offs2(i),diff_params(i,1),'k:','LineWidth',line_norm);
end
plot([0 10], zeros(1,2),'k--','LineWidth',2)
h2 = plot(NaN,NaN,'k--','LineWidth',2);
% %%%plot(offs2(1),MRL(1,1) - mRL(1,1),'Marker',sub_mark(str2num(patients{1}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',2,'LineStyle','none')
%plot((4:5)+offs2(i),params(:,4:5),'k-','LineWidth',line_norm);
% plot mean and std
title('Difference','FontSize',14)
hold off

axes(ha(3))
offs2 = 0.4*rand(length(patients),1)+1.2;
plot(NaN,NaN)
hold on
boxp = boxplot([mRL_WITHDis(1,:)' MRL_WITHDis(1,:)'],'Colors','k','Symbol','','whisker',1000);
set(boxp,'LineWidth',line_med)
for i = 1:length(patients)
    plot(offs2(i),mRL_WITHDis(1,i),'Marker',sub_mark(str2num(patients{i}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',line_med)
    plot(1+offs2(i),MRL_WITHDis(1,i),'Marker',sub_mark(str2num(patients{i}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',line_med)
    plot([offs2(i), 1+offs2(i)],[mRL_WITHDis(1,i) MRL_WITHDis(1,i)],'k','LineWidth',1,'LineStyle',':')
    %     errorbar(offs2(i),means(i),stds(i),'Marker','o','MarkerSize',7,'MarkerFaceColor',meancolor,'MarkerEdgeColor',meancolor,'Color',stdcolor,'CapSize',15,'LineWidth',2);
%     plot(4+offs2(i),diff_params(i,1),'k:','LineWidth',line_norm);
end
plot([0 10], ones(1,2),'k:','LineWidth',1.5)
% %%%plot([offs2(1), 1+offs2(1)],[mRL(1,1) MRL(1,1)],'Marker',sub_mark(str2num(patients{1}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',2,'LineStyle','none')
%plot((4:5)+offs2(i),params(:,4:5),'k-','LineWidth',line_norm);
% plot mean and std
title('Distraction','FontSize',14)
hold off

annotation('textarrow',[0.02,0.02],[0.5,1],'String','Time (ms)',...
    'TextRotation',90,'FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle',...
    'LineStyle','none','HeadStyle','none');
annotation('textarrow',[xmin+xwid(1)+0.02,0.02],[0.5,1],'String','\DeltaTime (ms)',...
    'TextRotation',90,'FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle',...
    'LineStyle','none','HeadStyle','none');
annotation('textarrow',[xmin+2*xwid(1)+xspacBIG+0.013,0.02],[0.5,1],'String','Time (ms)',...
    'TextRotation',90,'FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle',...
    'LineStyle','none','HeadStyle','none');

annotation('textbox',[0.13 .99 1 .02],'String','Response Latencies for Shortest and Longest Detected Stimulation Pauses','FontSize',16,...
    'HorizontalAlignment','left','EdgeColor','none','FontWeight','bold');

set(ha(1),'XTickLabel','Threshold','FontSize',16)
% xlabel(ha(2),'Walking Speed (mph)')
set(ha(1:3),'Xgrid','on','YGrid','on')
set(ha(1),'XLim',[0.8 2.7],'YLim',[0 1.7])
set(ha(2),'XLim',[0.8 2.7],'YLim',[-1.2 0.6])
set(ha(3),'XLim',[0.8 2.7],'YLim',[0 1.7])
set(ha(1),'XTick',1:2,'XTickLabel',{'Shortest SP','Longest SP'},'XTickLabelRotation',0,'box','on','FontSize',16)
set(ha(2),'XTick',1:5,'XTickLabel',{'Change'},'XTickLabelRotation',0,'box','on','FontSize',16)
set(ha(3),'XTick',1:2,'XTickLabel',{'Shortest SP','Longest SP'},'XTickLabelRotation',0,'box','on','FontSize',16)

% set(ha(3),'YAxisLocation','right')
% title(ha(2),'Response Latencies corresponding to','FontSize',14)
% title(ha(3),'Shortest and Longest Detected Stimulation Pauses','FontSize',14)
% leg3_reord = [1,6,2,7,3,8,4,9,5,10];
leg1 = legend(ha(1),h1,leg1_plist,'NumColumns',4,'Location','north','FontSize',16);
leg1.ItemTokenSize(1) = 15;
title(leg1,'Subjects')
%leg2 = legend(ha(2),h2,'Latency = 1s','NumColumns',1,'Location','north');
%leg2.ItemTokenSize(1) = 45;
leg2 = legend(ha(2),h2,' No difference','NumColumns',1,'Location','north','FontSize',16);
leg2.ItemTokenSize(1) = 55;
% set(ha(3),'XLim',[3.3 4.5],'YLim',YLim3)
set(ha(2),'XTick',1:2,'XTickLabel',{'Quiet','Disctraction'},'XTickLabelRotation',0,'box','on')
set(ha(2),'YAxisLocation','left')
% title(ha(3),'D. \DeltaDVAS Relative to Tonic')
%ylabel(ha(2),'\DeltaTime (s)')
% leg5.Position = [0.6966    0.8737    0.2040    0.0820];
ha(1).XAxis.FontSize = 16;
ha(1).Title.FontSize = 14;
ha(2).XAxis.FontSize = 16;
ha(2).Title.FontSize = 14;
ha(3).XAxis.FontSize = 16;
ha(3).Title.FontSize = 14;


for i = 1:3
    ha(i).Position = [xpos(i),ypos(i),xwid(i),ywid];
end


% set(ha(3),'Position',[xpos(2),ypos(1),xwid2,ymax-ymin])
% set(ha(3),'Position',[xpos(3),ypos(1),xwid2,ymax-ymin])
% set(ha(1:2),'XTick',0:0.5:3,'XLim',XLim,'YLim',YLim,'XGrid','on','YGrid','on','box','on')

%% Different take on the tiled plot: population cdf, and tiles around
xscale = 'log';
sub_mark = 'xdo^ps+hv<'; %MVI001-MVI010
stdcolor = [192 0 0]./250;
meancolor = [192 0 0]./250;
fig2 = figure(2);
clf;
set(fig2,'Color',[1,1,1],'Position',[115   172   978   1000])
ha = gobjects(1,8);
ha_pos = NaN(8,4);
restplots = [4 8 12 13 14 15 16];
cont = -1;
for i = 1:8
    cont = cont + 1;
    if i == 1 
        ha(i) = subplot(4,4,[1 2 3 5 6 7 9 10 11]);
    else
        ha(i) = subplot(4,4,restplots(cont));
    end
    ha_pos(i,:)=get(ha(i),'Position'); 
end
xmin = 0.06;
xmax = 0.99;
xspac = 0.015;
ymin = 0.08;
ymax = 0.95;
yspac = 0.03;
ywid = (ymax-ymin-3*yspac)/4;
xwid = (xmax-xmin-3*xspac)/4;
xpos = repmat(xmin:(xwid+xspac):xmax,1,4);
ypos = flip(reshape(repmat(ymin:(ywid+yspac):ymax,4,1),[],1)',2);

yspacBIG = 0.03;
xspacBIG = 0.02;

for i = 1:8
    if i == 1
        ha(i).Position = [xpos(i),ypos(end-4)+yspacBIG,3*xwid + 2*xspac - xspacBIG ,3*ywid + 2*yspac - yspacBIG];
    elseif i>1 && i<5
        ha(i).Position = [xpos(4),ypos(4*(i-1)),xwid,ywid];
    else
        ha(i).Position = [xpos(8+i),ypos(8+i),xwid,ywid];
    end
end
set(ha,'box','on')
xdataALL = [];
xdataPOOLED = [0 4 9 19 29 39 59 79 99 198 298]; % same as xdataALL but binning similar pauses like 297 ms and 298 ms for easier calculations
markALL = [];
ydataALL = [];
ydataPOOLED = zeros(length(patients),length(xdataPOOLED));
repdataPOOLED = ydataPOOLED; % number of reps per duration and per subject
if strcmp(xscale,'linear')
    jitf = 0.05;
%     jitf = 0;
elseif strcmp(xscale,'log')
    jitf = [0 0.4 1.2 1.5 1.3 3 3 3 5 12 15]/50;
%     jitf = 0;
end
cALL = [];
order = [5:8 4:-1:2];
RL = NaN(length(patients),length(xdataPOOLED));
for i = 1:length(patients)+1
    if i <= length(patients)
        sub_id = str2num(patients{i}(4:6));
        load([patients{i},'.mat'])
        xdataALL = [xdataALL; xdata];
        ydataALL = [ydataALL; ydata];
        ydataPOOLED(i,:) = detectionsummary(:,3)';
        repdataPOOLED(i,:) = detectionsummary(:,2)';
        markALL = [markALL repmat(sub_mark(sub_id),1,length(xdata))];
        fun = @(x,xdata)cdf('Normal',xdata,x(1),x(2));
        x0 = [1,20,10];
        x = lsqcurvefit(fun,x0,xdata,ydata);
        times = linspace(xdata(1),xdata(end));
%         [MRL(1,i),Midx] = max(detectionsummary(:,5));
%         [mRL(1,i),midx] = min(detectionsummary(:,5));
%         MRL(2,i) = detectionsummary(Midx,6);
%         mRL(2,i) = detectionsummary(midx,6);
        RL(i,:) = detectionsummary(:,5)'; % RL column removing initial NaN values
        omitnanRL = RL(i,~isnan(RL(i,:)));
        MRL(1,i) = omitnanRL(end); % RL associated with longest duration stim pause detected
        mRL(1,i) = omitnanRL(1); % RL associated with shortest duration stim pause detected
        DRL(1,i) = omitnanRL(end) - omitnanRL(1);
        axes(ha(order(i)))
        means(i) = x(1);
        stds(i) = x(2);
        if i == 1
            plot(xdata,ydata,'k','Marker',sub_mark(sub_id),'MarkerFaceColor','k','LineWidth',1,'LineStyle','none')
            set(gca,'Xscale',xscale)
        else
            plot(xdata,ydata,'k','Marker',sub_mark(sub_id),'MarkerFaceColor','none','LineStyle','none')
            set(gca,'Xscale',xscale)
        end            
        hold on
        plot(times,fun(x,times),'k-','LineWidth',2)
        %         plot(x(1),0.5,'x','MarkerSize',10,'LineWidth',3,'Color',[.5 .5 .5])
        %         plot(x(1) + [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',[.5 .5 .5])
        %         plot(x(1) - [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',[.5 .5 .5])
        plot(x(1),0.5,'*','MarkerSize',10,'LineWidth',1.5,'Color',meancolor)
        plot(x(1) + [x(2) x(2)]/2,[-10 1.5],':','LineWidth',2,'Color',stdcolor)
        plot(x(1) - [x(2) x(2)]/2,[-10 1.5],':','LineWidth',2,'Color',stdcolor)
        hl(order(i),1) = plot(NaN,NaN,'*','MarkerSize',10,'LineWidth',1.5,'Color',meancolor);
        hl(order(i),2) = plot(NaN,NaN,':','LineWidth',2,'Color',stdcolor);
        grid on
        grid minor
        %legend('Data',sprintf('Cum Gaussian Dist Fit (mean %.3g.2,sd %.3g)',x(1),x(2)));
        lss = legend(ha(order(i)),hl(order(i),:),['\mu ' sprintf('= %.3g ms',x(1))],['\sigma ' sprintf('= %.3g ms',x(2))],'FontSize',12);
        if strcmp(xscale,'log')
            lss.Location = 'northwest';
        elseif strcmp(xscale,'linear')
            lss.Location = 'southeast';
        end
        % title(sprintf('Detection Rate: %s analyzed %s\nusing %s',filename,char(datetime),nameofthismfile),'FontSize',8);
        title(patients{i},'FontSize',14)
        set(gcf,'color','w');
        box on
        axis([0 1.1*max(xdata) -0.05 1.05]);
        cALL(i,:) = fun(x,times);
    else
        fun = @(x,xdataPOOLED)cdf('Normal',xdataPOOLED,x(1),x(2));
        x0 = [1,20,10];
        x = lsqcurvefit(fun,x0,xdataPOOLED,sum(ydataPOOLED,1)./sum(repdataPOOLED,1)); %%%%%%%%%%%%%%%% this way of computing mean? and std?
        xp = lsqcurvefit(fun,x0,xdataPOOLED,sum(ydataPOOLED,1)./sum(repdataPOOLED,1) + std(ydataPOOLED./repdataPOOLED,1));
        xm = lsqcurvefit(fun,x0,xdataPOOLED,sum(ydataPOOLED,1)./sum(repdataPOOLED,1) - std(ydataPOOLED./repdataPOOLED,1));
%         times = linspace(xdataPOOLED(1),xdataPOOLED(end));
        times = xdataPOOLED(1):0.01:xdataPOOLED(end);
        meancdf = fun(x,times);
        meancdf_pts = meancdf(100*xdataPOOLED+1);
        prms = sum((ydataPOOLED./repdataPOOLED - meancdf_pts).^2,1)/length(xdataPOOLED);
%         xp =meancdf_pts+sum(prms);
%         xm = meancdf_pts-sum(prms);
        axes(ha(1))
        means(i) = x(1);
        stds(i) = x(2);
        for j = 1:size(ydataPOOLED,1)
            sub_id = str2num(patients{j}(4:6));
            if j == 1
                plot(xdataPOOLED+(jitf)*(1 + 10*abs(j-3))*sign(j-3),ydataPOOLED(j,:)./repdataPOOLED(j,:),'k','Marker',sub_mark(sub_id),'MarkerFaceColor','k','LineWidth',1,'LineStyle','none')
            else
                plot(xdataPOOLED+(jitf)*(1 + 10*abs(j-3))*sign(j-3),ydataPOOLED(j,:)./repdataPOOLED(j,:),'k','Marker',sub_mark(sub_id),'MarkerFaceColor','none','LineStyle','none')
            end
            hold on
        end
        set(gca,'Xscale',xscale)
        plot(times,fun(x,times),'k-','LineWidth',2)
        times(1) = 0.001;
        fill([times fliplr(times)], [fun(xp,times) fliplr(fun(xm,times))],[.5 .5 .5],'facealpha',.25);
        xdataPOOLED(1) = xdataPOOLED(1) + 0.001;
%         fill([xdataPOOLED fliplr(xdataPOOLED) xdataPOOLED(1)], [xp fliplr(xm) xp(1)],[.5 .5 .5],'facealpha',.25);
        %         plot(x(1),0.5,'x','MarkerSize',10,'LineWidth',3,'Color',[.5 .5 .5])
        %         plot(x(1) + [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',[.5 .5 .5])
        %         plot(x(1) - [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',[.5 .5 .5])
        plot(x(1),0.5,'*','MarkerSize',10,'LineWidth',1.5,'Color',meancolor)
        plot(x(1) + [x(2) x(2)]/2,[-10 1.5],':','LineWidth',2,'Color',stdcolor)
        plot(x(1) - [x(2) x(2)]/2,[-10 1.5],':','LineWidth',2,'Color',stdcolor)
        hl(1,1) = plot(NaN,NaN,'*','MarkerSize',10,'LineWidth',1.5,'Color',meancolor);
        hl(1,2) = plot(NaN,NaN,':','LineWidth',2,'Color',stdcolor);
        hl2(1) = plot(NaN,NaN,'LineWidth',2,'Color','k');
        hl2(2) = fill(NaN,NaN,[.5 .5 .5],'facealpha',.25);
        grid on
        grid minor
        %legend('Data',sprintf('Cum Gaussian Dist Fit (mean %.3g.2,sd %.3g)',x(1),x(2)));
        lps = legend(ha(1),[hl(1,:) hl2],['\mu ' sprintf('= %.3g ms',x(1))],['\sigma ' sprintf('= %.3g ms',x(2))],'CDF fit to group \mu','CDF fit to group \sigma','Location','northwest','FontSize',12);
        % title(sprintf('Detection Rate: %s analyzed %s\nusing %s',filename,char(datetime),nameofthismfile),'FontSize',8);
        if strcmp(xscale,'log')
            lps.Location = 'northwest';
            lps.Position = [0.1495 0.855 0.05 0.05];
            lps.FontSize = 12.05;
            lps.Location = 'southeast';
        elseif strcmp(xscale,'linear')
            lps.Location = 'southeast';
        end
        title('Pooled Data for All Participants','FontSize',14)
%         title('Linear Scale','FontSize',14)
%         xlabel('Stimulation Pause Duration (ms)','FontSize',14)
%         ylabel('Fraction Detected','FontSize',14)
        set(gcf,'color','w');
        box on
        axis([3 1.1*max(xdataALL) -0.05 1.05]);
    end
end

set(ha(2:4),'XTickLabel',[])
set(ha(6:8),'YTickLabel',[])
% xlabel(ha(4:6),'Stimulation Pause Duration (ms)','FontSize',14)
% ylabel(ha([1 4]),{'Fraction Detected'},'FontSize',14)
for i = 1:length(patients)+1
    ha(i).XAxis.FontSize = 13;
    ha(i).YAxis.FontSize = 13;
end

annotation('textbox',[0.417 .02 1 .02],'String','Stimulation Pause Duration (ms)','FontSize',14,...
    'HorizontalAlignment','left','EdgeColor','none');
annotation('textarrow',[0.02,0.02],[0.5,1],'String','Fraction Detected',...
    'TextRotation',90,'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle',...
    'LineStyle','none','HeadStyle','none');
% For insert figure
% axes(ha(1))
% title('Linear Scale','FontSize',14)
% xlabel('Stimulation Pause Duration (ms)','FontSize',14)
% ylabel('Fraction Detected','FontSize',14)
% axes(ha(5))
% title('')
% axes(ha(6))
% title('')
% axes(ha(7))
% title('')


%% Population RLs, and tiles around
fig3 = figure(3);
xscale = 'log';
set(fig3,'Color',[1,1,1],'Position',[115   172   978   1000])
ha = gobjects(1,8);
ha_pos = NaN(8,4);
restplots = [4 8 12 13 14 15 16];
cont = -1;
for i = 1:8
    cont = cont + 1;
    if i == 1 
        ha(i) = subplot(4,4,[1 2 3 5 6 7 9 10 11]);
    else
        ha(i) = subplot(4,4,restplots(cont));
    end
    ha_pos(i,:)=get(ha(i),'Position'); 
end
xmin = 0.06;
xmax = 0.99;
xspac = 0.015;
ymin = 0.08;
ymax = 0.95;
yspac = 0.03;
ywid = (ymax-ymin-3*yspac)/4;
xwid = (xmax-xmin-3*xspac)/4;
xpos = repmat(xmin:(xwid+xspac):xmax,1,4);
ypos = flip(reshape(repmat(ymin:(ywid+yspac):ymax,4,1),[],1)',2);

yspacBIG = 0.03;
xspacBIG = 0.02;

for i = 1:8
    if i == 1
        ha(i).Position = [xpos(i),ypos(end-4)+yspacBIG,3*xwid + 2*xspac - xspacBIG ,3*ywid + 2*yspac - yspacBIG];
    elseif i>1 && i<5
        ha(i).Position = [xpos(4),ypos(4*(i-1)),xwid,ywid];
    else
        ha(i).Position = [xpos(8+i),ypos(8+i),xwid,ywid];
    end
end
set(ha,'box','on')
annotation('textbox',[0.417 .02 1 .02],'String','Stimulation Pause Duration (ms)','FontSize',14,...
    'HorizontalAlignment','left','EdgeColor','none');
annotation('textarrow',[0.02,0.02],[0.5,1],'String','Response Latency (s)',...
    'TextRotation',90,'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle',...
    'LineStyle','none','HeadStyle','none');


xdataPOOLED = [0 4 9 19 29 39 59 79 99 198 298]; % same as xdataALL but binning similar pauses like 297 ms and 298 ms for easier calculations
rlPOOLED = zeros(length(patients),length(xdataPOOLED)); % response latencies
order = [5:8 4:-1:2];
for i = 1:length(patients)+1
    if i <= length(patients)
        load([patients{i},'.mat'])
        sub_id = str2num(patients{i}(4:6));
        rlPOOLED(i,:) = detectionsummary(:,5)';
        axes(ha(order(i)))
        if i == 1
            errorbar(xdata,detectionsummary(:,5),detectionsummary(:,6),'k','Marker',sub_mark(sub_id),'MarkerFaceColor','none','LineWidth',1)
        else
            errorbar(xdata,detectionsummary(:,5),detectionsummary(:,6),'k','Marker',sub_mark(sub_id),'MarkerFaceColor','none','LineWidth',1)        
        end
        hold on
%         errorbar(xdata,detectionsummary(:,5),detectionsummary(:,6),'k','Marker','none','LineWidth',1.5)
        hla(1,2) = plot([0.01 400],ones(1,2),'k:','LineWidth',1.5);
        hla(1,1) = errorbar(NaN,10,'k','Marker','*','MarkerFaceColor','k','LineWidth',1);
        %         plot(x(1),0.5,'x','MarkerSize',10,'LineWidth',3,'Color',[.5 .5 .5])
        %         plot(x(1) + [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',[.5 .5 .5])
        %         plot(x(1) - [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',[.5 .5 .5])
%         plot(x(1),0.5,'x','MarkerSize',10,'LineWidth',3,'Color',meancolor)
%         plot(x(1) + [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',stdcolor)
%         plot(x(1) - [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',stdcolor)
%         hl(order(i),1) = plot(NaN,NaN,'x','MarkerSize',10,'LineWidth',3,'Color',meancolor);
%         hl(order(i),2) = plot(NaN,NaN,':','LineWidth',2,'Color',stdcolor);
%         grid on
        grid minor
        %legend('Data',sprintf('Cum Gaussian Dist Fit (mean %.3g.2,sd %.3g)',x(1),x(2)));
%         legend(ha(order(i)),hl(order(i),:),['\mu ' sprintf('= %.3g ms',x(1))],['\sigma ' sprintf('= %.3g ms',x(2))],'Location','southeast','FontSize',12)
        % title(sprintf('Detection Rate: %s analyzed %s\nusing %s',filename,char(datetime),nameofthismfile),'FontSize',8);
        legend(ha(order(1)),hla,'Subject \mu and \sigma','Latency = 1 s','Location','southwest','FontSize',12)
        title(patients{i},'FontSize',14)
        set(gcf,'color','w');
        box on
        set(gca,'Xscale',xscale)
        if strcmp(xscale,'linear')
            axis([0 1.1*max(xdataPOOLED) 0.21 1.3]);
            jitf = 3;
            jitf = 0;
        elseif strcmp(xscale,'log')
            axis([3.5 1.1*max(xdataPOOLED) 0.205 1.3]);
            jitf = [0 0.3 0.5 1 1.3 3 5 7 10 12 15];
            jitf = 0;
        end
%         cALL(i,:) = fun(x,times);
    else
        times = linspace(xdataPOOLED(1),xdataPOOLED(end));
        axes(ha(1))
        means(i) = x(1);
        stds(i) = x(2);
        for j = 1:size(rlPOOLED,1)
            sub_id = str2num(patients{j}(4:6));
            if j == 1
                plot(xdataPOOLED + jitf,rlPOOLED(j,:),'k','Marker',sub_mark(sub_id),'MarkerFaceColor','k','LineWidth',0.7,'LineStyle','none')
            else
                plot(xdataPOOLED + jitf,rlPOOLED(j,:),'k','Marker',sub_mark(sub_id),'MarkerFaceColor','none','LineStyle','none')
            end
            hold on
        end
%         errorbar(xdataPOOLED,mean(rlPOOLED,1,'omitnan'),std(rlPOOLED,1,'omitnan'),'k','Marker','.','MarkerFaceColor','k','LineWidth',2);
        plot(xdataPOOLED,mean(rlPOOLED,1,'omitnan'),'k','Marker','.','MarkerFaceColor','k','LineWidth',2);
        xdataPOOLED(1) = 0.01;
        upp = mean(rlPOOLED,1,'omitnan')+std(rlPOOLED,1,'omitnan')./sum(~isnan(rlPOOLED));
        loww = fliplr(mean(rlPOOLED,1,'omitnan')-std(rlPOOLED,1,'omitnan')./sum(~isnan(rlPOOLED)));
        fill([xdataPOOLED(1:5) fliplr(xdataPOOLED(1:5))], [upp(1:5) loww(end-4:end)],[.5 .5 .5],'facealpha',.1,'LineStyle','none');
        fill([xdataPOOLED(5:end) fliplr(xdataPOOLED(5:end))], [upp(5:end) loww(1:end-4)],[.5 .5 .5],'facealpha',.4,'LineStyle','none');
        fill([xdataPOOLED(1:end) fliplr(xdataPOOLED(1:end))], [upp(1:end) loww(1:end)],[.5 .5 .5],'facealpha',0);
        plot([0.01 400],ones(1,2),'k:','LineWidth',1.5)
        %         plot(x(1),0.5,'x','MarkerSize',10,'LineWidth',3,'Color',[.5 .5 .5])
        %         plot(x(1) + [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',[.5 .5 .5])
        %         plot(x(1) - [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',[.5 .5 .5])
        hl(1,1) = errorbar(NaN,NaN,'k','Marker','.','MarkerFaceColor','k','LineWidth',2);
        hl(1,2) = plot(NaN,NaN,'k:','LineWidth',1.5);
        hl2(1) = fill(NaN,NaN,[.5 .5 .5],'facealpha',.25);
        grid on
        grid minor
        %legend('Data',sprintf('Cum Gaussian Dist Fit (mean %.3g.2,sd %.3g)',x(1),x(2)));
        legend(ha(1),[hl(1,1) hl2(1) hl(1,2)],'Group \mu','\pmSEM','Latency = 1 s','Location','southwest','FontSize',12)
        % title(sprintf('Detection Rate: %s analyzed %s\nusing %s',filename,char(datetime),nameofthismfile),'FontSize',8);
        title('Pooled Data for All Participants','FontSize',14)
        set(gcf,'color','w');
        box on
        set(gca,'Xscale',xscale)
        if strcmp(xscale,'linear')
            axis([0 1.1*max(xdataPOOLED) 0.35 1.1]);
        elseif strcmp(xscale,'log')
            axis([3.5 1.1*max(xdataPOOLED) 0.35 1.1]);
        end
    end
end
for i = 1:length(patients)+1
    ha(i).XAxis.FontSize = 13;
    ha(i).YAxis.FontSize = 13;
end
set(ha(2:4),'XTickLabel',[])
set(ha(6:8),'YTickLabel',[])
% xlabel(ha(4:6),'Stimulation Pause Duration (ms)','FontSize',14)
% ylabel(ha([1 4]),{'Fraction Detected'},'FontSize',14)

%% Side-by-side comparison for MVI007R765
figc = figure;
set(figc,'Color',[1,1,1],'Position',[440  378  1000   450])
ha = gobjects(1,2);
for i = 1:2
   ha(i) = subplot(1,2,i); 
end
xmin = 0.05;
xmax = 0.99;
xspac = 0.02;
ymin = 0.08;
ymax = 0.95;
yspac = 0.03;
ywid = (ymax-ymin);
xwid = repmat((xmax-xmin-xspac)/2,1,2);
xpos = [xmin xmin+xwid(1)+xspac];
ypos = flip(reshape(repmat(ymin:(ywid+yspac):ymax,1,2),[],1)',2);

for i = 1:2
    ha(i).Position = [xpos(i),ypos(i),xwid(i),ywid];
    axes(ha(i))
    if i == 1
        load('MVI007R765NoDistractions.mat')
        title('MVI007R765 without Distractions','FontSize',14)
    else
        load('MVI007R765WithDistractions.mat')
        title('MVI007R765 with Distractions','FontSize',14)
    end
    hold on
    fun = @(x,xdata)cdf('Normal',xdata,x(1),x(2));
    x0 = [1,20,10];
    x = lsqcurvefit(fun,x0,xdata,ydata);
    times = linspace(xdata(1),xdata(end));
    plot(xdata,ydata,'ko','MarkerFaceColor','k','Marker',sub_mark(7),'MarkerSize',7)
    plot(times,fun(x,times),'k-','LineWidth',2)
    hs(1) = plot(x(1),0.5,'*','MarkerSize',10,'LineWidth',1.5,'Color',meancolor);
    hs(2) = plot(x(1) + [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',stdcolor);
    plot(x(1) - [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',stdcolor)
    grid on
    grid minor
    ha(i).XAxis.FontSize = 13;
    ha(i).YAxis.FontSize = 13;
    %legend('Data',sprintf('Cum Gaussian Dist Fit (mean %.3g.2,sd %.3g)',x(1),x(2)));
    legend(ha(i),hs,['\mu ' sprintf('= %.3g ms',x(1))],['\sigma ' sprintf('= %.3g ms',x(2))],'Location','southeast','FontSize',12)
%     ylabel('Fraction Detected','FontSize',14);
%     xlabel('Stimulation Pause Duration (ms)','FontSize',14);
    % title(sprintf('Detection Rate: %s analyzed %s\nusing %s',filename,char(datetime),nameofthismfile),'FontSize',8);
    set(gcf,'color','w');
    box on
    axis([0 1.01*max(xdataALL) 0 1.05]);
    % set(gca,'Xscale','log')
end
annotation('textbox',[0.417 .02 1 .02],'String','Stimulation Pause Duration (ms)','FontSize',14,...
    'HorizontalAlignment','left','EdgeColor','none');
annotation('textarrow',[0.015,0.02],[0.5,1],'String','Fraction Detected',...
    'TextRotation',90,'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle',...
    'LineStyle','none','HeadStyle','none');
set(ha(2),'YTickLabels','')


%% Alternative tiled layout of psychometric curves

fig1 = figure(1);
clf;
set(fig1,'Color',[1,1,1],'Position',[115   172   978   1000])
ha = gobjects(1,6);
for i = 1:6
    ha(i) = subplot(3,2,i);
end
xmin = 0.06;
xmax = 0.99;
xspac = 0.01;
ymin = 0.08;
ymax = 0.92;
yspac = 0.02;
ywid = (ymax-ymin-yspac)/2;
xwid = (xmax-xmin-2*xspac)/3;
xpos = repmat(xmin:(xwid+xspac):xmax,1,2);
ypos = reshape(repmat([(ymax-ywid),ymin],3,1),[],1)';
for i = 1:6
    ha(i).Position = [xpos(i),ypos(i),xwid,ywid];
end
set(ha,'box','on')
% annotation('textbox',[0 .9 1 .1],'String',sub,'FontSize',20,...
%     'HorizontalAlignment','left','EdgeColor','none','FontWeight','bold');

for i = 1:length(patients)
    if i == 1
        load([patients{i},'.mat'])
        fun = @(x,xdata)cdf('Normal',xdata,x(1),x(2));
        x0 = [1,20,10];
        x = lsqcurvefit(fun,x0,xdata,ydata);
        times = linspace(xdata(1),xdata(end));
        means(i) = x(1);
        stds(i) = x(2);
    else
        axes(ha(i-1))
        load([patients{i},'.mat'])
        fun = @(x,xdata)cdf('Normal',xdata,x(1),x(2));
        x0 = [1,20,10];
        x = lsqcurvefit(fun,x0,xdata,ydata);
        times = linspace(xdata(1),xdata(end));
        means(i) = x(1);
        stds(i) = x(2);
        plot(xdata,ydata,'ko','MarkerFaceColor','k')
        hold on
        plot(times,fun(x,times),'k-','LineWidth',2)
        plot(x(1),0.5,'x','MarkerSize',10,'LineWidth',3,'Color',[.5 .5 .5])
        plot(x(1) + [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',[.5 .5 .5])
        plot(x(1) - [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',[.5 .5 .5])
        hl(i-1,1) = plot(NaN,NaN,'x','MarkerSize',10,'LineWidth',3,'Color',[.5 .5 .5]);
        hl(i-1,2) = plot(NaN,NaN,':','LineWidth',2,'Color',[.5 .5 .5]);
        grid on
        grid minor
        %legend('Data',sprintf('Cum Gaussian Dist Fit (mean %.3g.2,sd %.3g)',x(1),x(2)));
        legend(ha(i-1),hl(i-1,:),['\mu ' sprintf('= %.3g ms',x(1))],['\sigma ' sprintf('= %.3g ms',x(2))],'Location','southeast')
        % title(sprintf('Detection Rate: %s analyzed %s\nusing %s',filename,char(datetime),nameofthismfile),'FontSize',8);
        title(patients{i})
        set(gcf,'color','w');
        box on
        axis([0 1.1*max(xdata) 0 1.05]);
    end
end

set(ha(1:3),'XTickLabel',[])
set(ha([2,3,5,6]),'YTickLabel',[])
xlabel(ha(4:6),'Stimulation Pause Duration (ms)')
ylabel(ha([1 4]),{'Fraction Detected'})
%% Alternative with curves stacked on top of one another
figure
col_order = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F'};
for i = 1:length(patients)
    load([patients{i},'.mat'])
    fun = @(x,xdata)cdf('Normal',xdata,x(1),x(2));
    x0 = [1,20,10];
    x = lsqcurvefit(fun,x0,xdata,ydata);
    times = linspace(xdata(1),xdata(end));
    plot(xdata,ydata,'o','MarkerEdgeColor',col_order{i},'MarkerFaceColor',col_order{i})
    hold on
    plot(times,fun(x,times),'-','LineWidth',2,'Color',col_order{i})
    plot(x(1),0.5,'x','MarkerSize',10,'LineWidth',3,'Color',[.5 .5 .5])
%     plot(x(1) + [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',[.5 .5 .5])
%     plot(x(1) - [x(2) x(2)]/2,[0 1.5],':','LineWidth',2,'Color',[.5 .5 .5])
end
grid on
grid minor
%legend('Data',sprintf('Cum Gaussian Dist Fit (mean %.3g.2,sd %.3g)',x(1),x(2)));
% legend('Data','CDF Fit',['\mu ' sprintf('= %.3g ms',x(1))],['\sigma ' sprintf('= %.3g ms',x(2))],'Location',[0.685, 0.15, .15, .1])
ylabel('Fraction Detected');
xlabel('Stimulation Pause Duration (ms)');
% title(sprintf('Detection Rate: %s analyzed %s\nusing %s',filename,char(datetime),nameofthismfile),'FontSize',8);
title('')
set(gcf,'color','w');
box on
axis([0 1.1*max(xdata) 0 1.05]);

%% Boxplots for Summary

for i = 1:length(patients)
    leg1_plist{i} = num2str(str2num(patients{i}(4:6)));
    leg1_mlist(i) = sub_mark(str2num(patients{i}(4:6)));
end

line_norm = 0.5;
line_bold = 2;
line_med = 1.5;
mark_size_big = 16;
mark_size_med = 5;
fig4 = figure(4);
set(fig4,'Color',[1,1,1],'Position',[115   172   978   1000])
ha = gobjects(1,3);
for i = 1:3
   ha(i) = subplot(1,3,i); 
end
xmin = 0.06;
xmax = 0.99;
xspac = 0.01;
xspacBIG = 0.07;
ymin = 0.08;
ymax = 0.95;
yspac = 0.03;
ywid = (ymax-ymin);
xwid = repmat((xmax-xmin-xspac -xspacBIG)/3,1,3);
xwid(3) = xwid(3)*0.7;
xpos = [xmin xmin+xwid(1)+xspacBIG xmin+2*xwid(1)+xspacBIG+xspac];
ypos = flip(reshape(repmat(ymin:(ywid+yspac):ymax,1,3),[],1)',2);

axes(ha(1))
offs2 = 0.35*rand(length(patients),1)+1.1;
plot(NaN,NaN)
hold on
boxp = boxplot(means,'Colors','k','Symbol','','whisker',1000);
set(boxp,'LineWidth',line_med)
for i = 1:length(patients)
    plot(offs2(i),means(i),'Marker',sub_mark(str2num(patients{i}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor',meancolor,'LineWidth',line_med)
    %     errorbar(offs2(i),means(i),stds(i),'Marker','o','MarkerSize',7,'MarkerFaceColor',meancolor,'MarkerEdgeColor',meancolor,'Color',stdcolor,'CapSize',15,'LineWidth',2);
    plot(4+offs2(i),means(i),'k:','LineWidth',line_norm);
%     plot([0 10], zeros(1,2),'k--','LineWidth',1)
    h1(i) = plot(NaN,NaN,'Marker',sub_mark(str2num(patients{i}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',line_med,'LineStyle','none');
end
% %%%plot(offs2(1),means(1),'Marker',sub_mark(str2num(patients{1}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor',meancolor,'LineWidth',2)
%plot((4:5)+offs2(i),params(:,4:5),'k-','LineWidth',line_norm);
% plot mean and std
hold off

axes(ha(2))
offs2 = 0.4*rand(length(patients),1)+1.2;
plot(NaN,NaN)
hold on
boxp = boxplot([mRL(1,:)' MRL(1,:)'],'Colors','k','Symbol','','whisker',1000);
set(boxp,'LineWidth',line_med)
for i = 1:length(patients)
    plot(offs2(i),mRL(1,i),'Marker',sub_mark(str2num(patients{i}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',line_med)
    plot(1+offs2(i),MRL(1,i),'Marker',sub_mark(str2num(patients{i}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',line_med)
    plot([offs2(i), 1+offs2(i)],[mRL(1,i) MRL(1,i)],'k','LineWidth',1,'LineStyle',':')
    %     errorbar(offs2(i),means(i),stds(i),'Marker','o','MarkerSize',7,'MarkerFaceColor',meancolor,'MarkerEdgeColor',meancolor,'Color',stdcolor,'CapSize',15,'LineWidth',2);
%     plot(4+offs2(i),diff_params(i,1),'k:','LineWidth',line_norm);
end
plot([0 10], ones(1,2),'k:','LineWidth',1.5)
h2 = plot(NaN,NaN,'k:','LineWidth',1.5);
% %%%plot([offs2(1), 1+offs2(1)],[mRL(1,1) MRL(1,1)],'Marker',sub_mark(str2num(patients{1}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',2,'LineStyle','none')
plot([offs2(1), 1+offs2(1)],[mRL(1,1) MRL(1,1)],'k','LineWidth',1,'LineStyle',':')
%plot((4:5)+offs2(i),params(:,4:5),'k-','LineWidth',line_norm);
% plot mean and std
hold off

axes(ha(3))
offs2 = 0.35*rand(length(patients),1)+1.1;
plot(NaN,NaN)
hold on
boxp = boxplot(DRL(1,:),'Colors','k','Symbol','','whisker',1000);
set(boxp,'LineWidth',line_med)
for i = 1:length(patients)
    plot(offs2(i),MRL(1,i) - mRL(1,i),'Marker',sub_mark(str2num(patients{i}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',line_med)
    %     errorbar(offs2(i),means(i),stds(i),'Marker','o','MarkerSize',7,'MarkerFaceColor',meancolor,'MarkerEdgeColor',meancolor,'Color',stdcolor,'CapSize',15,'LineWidth',2);
%     plot(4+offs2(i),diff_params(i,1),'k:','LineWidth',line_norm);
end
plot([0 10], zeros(1,2),'k--','LineWidth',2)
h3 = plot(NaN,NaN,'k--','LineWidth',2);
% %%%plot(offs2(1),MRL(1,1) - mRL(1,1),'Marker',sub_mark(str2num(patients{1}(4:6))),'MarkerSize',mark_size_big,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',2,'LineStyle','none')
%plot((4:5)+offs2(i),params(:,4:5),'k-','LineWidth',line_norm);
% plot mean and std
hold off

annotation('textarrow',[0.02,0.02],[0.5,1],'String','Time (ms)',...
    'TextRotation',90,'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle',...
    'LineStyle','none','HeadStyle','none');
annotation('textarrow',[xmin+xwid(1)+0.035,0.02],[0.5,1],'String','Time (s)',...
    'TextRotation',90,'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle',...
    'LineStyle','none','HeadStyle','none');

annotation('textbox',[0.40 .9613 1 .02],'String','Response Latencies for Shortest and Longest Detected Stimulation Pauses','FontSize',14,...
    'HorizontalAlignment','left','EdgeColor','none','FontWeight','bold');

set(ha(1),'XTickLabel','Threshold','FontSize',14)
title(ha(1),'Detection Thresholds for','Stimulation Pause Duration','FontSize',14,'FontWeight','bold')
% xlabel(ha(2),'Walking Speed (mph)')
set(ha(1:3),'Xgrid','on','YGrid','on')
set(ha(1),'XLim',[0.8 1.5],'YLim',[28 87])
set(ha(2),'XLim',[0.8 2.7],'YLim',[0.45 1.1])
set(ha(3),'XLim',[0.8 1.5],'YLim',[-0.4 0.0727])
set(ha(2),'XTick',1:2,'XTickLabel',{'Shortest SP','Longest SP'},'XTickLabelRotation',0,'box','on','FontSize',14)
set(ha(3),'XTick',1:5,'XTickLabel',{'Change'},'XTickLabelRotation',0,'box','on','FontSize',14)
% set(ha(3),'YAxisLocation','right')
% title(ha(2),'Response Latencies corresponding to','FontSize',14)
% title(ha(3),'Shortest and Longest Detected Stimulation Pauses','FontSize',14)
% leg3_reord = [1,6,2,7,3,8,4,9,5,10];
leg1 = legend(ha(1),h1,leg1_plist,'NumColumns',7,'Location','north');
leg1.ItemTokenSize(1) = 8;
title(leg1,'Subjects')
leg2 = legend(ha(2),h2,'Latency = 1s','NumColumns',1,'Location','north');
leg2.ItemTokenSize(1) = 45;
leg3 = legend(ha(3),h3,'No difference','NumColumns',1,'Location','north');
leg3.ItemTokenSize(1) = 45;
% set(ha(3),'XLim',[3.3 4.5],'YLim',YLim3)
set(ha(3),'XTick',1,'XTickLabel',{'Difference'},'XTickLabelRotation',0,'box','on')
set(ha(3),'YAxisLocation','right')
% title(ha(3),'D. \DeltaDVAS Relative to Tonic')
ylabel(ha(3),'\DeltaTime (s)')
% leg5.Position = [0.6966    0.8737    0.2040    0.0820];
ha(3).XAxis.FontSize = 15;
ha(3).XAxis.FontSize = 15;

for i = 1:3
    ha(i).Position = [xpos(i),ypos(i),xwid(i),ywid];
end


% set(ha(3),'Position',[xpos(2),ypos(1),xwid2,ymax-ymin])
% set(ha(3),'Position',[xpos(3),ypos(1),xwid2,ymax-ymin])
% set(ha(1:2),'XTick',0:0.5:3,'XLim',XLim,'YLim',YLim,'XGrid','on','YGrid','on','box','on')


%% Now make a separate figure to plot distribution means and standard deviations
l = length(means);

fs = figure;
col_order = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F'};
plot([0 10],max(means)*ones(1,2),'-','LineWidth',1,'Color',0.5*ones(1,3))
hold on
plot([0 10],min(means)*ones(1,2),'-','LineWidth',1,'Color',0.5*ones(1,3))
fill([-1 10 10 -1],[min(means) min(means) max(means) max(means)],0.8*ones(1,3), 'FaceAlpha', 0.5)
for i = 1:l
    e = errorbar(i,means(i),stds(i));
    e.Marker = 'o';
    e.MarkerSize = 7;
    e.MarkerFaceColor = 'k';
    e.Color = 'k';
    e.CapSize = 15;
    e.LineWidth = 2;
end
plot([0 10],median(means)*ones(1,2),'k--','LineWidth',2)
hs(1) = errorbar(NaN,NaN,'Marker','o','MarkerSize',7,'MarkerFaceColor','k','Color','k','CapSize',15,'LineWidth',2);
% hs(2) = plot(NaN,NaN,'-','LineWidth',1,'Color',0.5*ones(1,3));
hs(2) = fill(NaN,NaN,0.8*ones(1,3), 'FaceAlpha', 0.5);
hs(3) = plot(NaN,NaN,'k--','LineWidth',2);
grid on
grid minor
ls = legend(hs,'Individual \mu and \sigma','Pooled \mu Range','Pooled \mu Median','Location','northwest');
ls.ItemTokenSize(1) = 30;
%legend('Data',sprintf('Cum Gaussian Dist Fit (mean %.3g.2,sd %.3g)',x(1),x(2)));
% legend('Data','CDF Fit',['\mu ' sprintf('= %.3g ms',x(1))],['\sigma ' sprintf('= %.3g ms',x(2))],'Location',[0.685, 0.15, .15, .1])
ylabel('Stimulation Pause Duration (ms)');
xlabel('Subject');
% xticklabels({'1','2','4','5','6','9','10'})
% title(sprintf('Detection Rate: %s analyzed %s\nusing %s',filename,char(datetime),nameofthismfile),'FontSize',8);
title('Distribution means and standard deviations')
set(gcf,'color','w');
box on
xlim([0.75 i+.25]);
ylim([min(means-stds)-10, max(means+stds)+10])
xticklabels({'',patients{1}(1:6),patients{2}(1:6),patients{3}(1:6),patients{4}(1:6),patients{5}(1:6),patients{6}(1:6),patients{7}(1:6)})
xtickangle(0)


%% Plotting distribution properties together with response latency properties

fig2 = figure(1);
clf;
set(fig2,'Color',[1,1,1],'Position',[115   172   978   1000])
ha = gobjects(1,3);
for i = 1:3
    ha(i) = subplot(3,1,i);
end
xmin = 0.08;
xmax = 0.99;
xspac = 0.01;
ymin = 0.08;
ymax = 0.94;
yspac = 0.03;
ywid = (ymax-ymin-2*yspac)/3;
xwid = (xmax-xmin-0*xspac)/1;
xpos = repmat(xmin:(xwid+xspac):xmax,1,3);
ypos = reshape(repmat([(ymax-ywid),ymax-2*ywid-yspac,ymin],1,1),[],1)';
for i = 1:3
    ha(i).Position = [xpos(i),ypos(i),xwid,ywid];
end
set(ha,'box','on')
% set(ha([2,3,5,6]),'YTickLabel',[])
% xlabel(ha(4:6),'Stimulation Pause Duration (ms)')
% ylabel(ha([1 4]),{'Fraction Detected'})

% First plot distribution properties
axes(ha(1))
l = length(means);
plot([0 10],max(means)*ones(1,2),'-','LineWidth',1,'Color',0.5*ones(1,3))
hold on
plot([0 10],min(means)*ones(1,2),'-','LineWidth',1,'Color',0.5*ones(1,3))
fill([-1 10 10 -1],[min(means) min(means) max(means) max(means)],0.8*ones(1,3), 'FaceAlpha', 0.5)
plot([0 10],median(means)*ones(1,2),'k--','LineWidth',2)
for i = 1:l
    e = errorbar(i,means(i),stds(i));
    e.Marker = 'o';
    e.MarkerSize = 7;
    e.MarkerFaceColor = meancolor;
    e.MarkerEdgeColor = meancolor;
    e.Color = stdcolor;
    e.CapSize = 15;
    e.LineWidth = 2;
end
hs(1) = errorbar(NaN,NaN,'Marker','o','MarkerSize',7,'MarkerFaceColor',meancolor,'MarkerEdgeColor',meancolor,'Color',stdcolor,'CapSize',15,'LineWidth',2);
% hs(2) = plot(NaN,NaN,'-','LineWidth',1,'Color',0.5*ones(1,3));
hs(2) = fill(NaN,NaN,0.8*ones(1,3), 'FaceAlpha', 0.5);
hs(3) = plot(NaN,NaN,'k--','LineWidth',2);
grid on
grid minor
ls = legend(hs,'Individual \mu and \sigma','Pooled \mu Range','Pooled \mu Median','Location','northwest','FontSize',12);
ls.ItemTokenSize(1) = 30;
%legend('Data',sprintf('Cum Gaussian Dist Fit (mean %.3g.2,sd %.3g)',x(1),x(2)));
% legend('Data','CDF Fit',['\mu ' sprintf('= %.3g ms',x(1))],['\sigma ' sprintf('= %.3g ms',x(2))],'Location',[0.685, 0.15, .15, .1])
ylabel('Stimulation Pause Duration (ms)','FontSize',14);
% xticklabels({'1','2','4','5','6','9','10'})
% title(sprintf('Detection Rate: %s analyzed %s\nusing %s',filename,char(datetime),nameofthismfile),'FontSize',8);
title('Distribution means and standard deviations','FontSize',14)
set(gcf,'color','w');
box on
ylim([min(means-stds)-10, max(means+stds)+10])
% xlim([0 l+1]);
xlim([0.75 l+0.25]);

% Then plot maximum response latency for each subject
axes(ha(2))
errorbar(1:l,MRL(1,:),MRL(2,:),'Marker','o','MarkerSize',7,'MarkerFaceColor','k','Color','k','CapSize',15,'LineWidth',2,'LineStyle','none')
hold on
plot([0 10],ones(1,2),'k--','LineWidth',2)
hm(1) = errorbar(NaN,NaN,'Marker','o','MarkerSize',7,'MarkerFaceColor','k','Color','k','CapSize',15,'LineWidth',2);
lm = legend(hm,'Maximum Response Latency (\mu and \sigma)','Location','southwest','FontSize',12);
lm.ItemTokenSize(1) = 30;
ylabel('Maximum Response Latency (s)','FontSize',14);
% xlim([0 l+1]);
xlim([0.75 l+0.25]);
grid on
grid minor
title('Response Latency Properties','FontSize',14)

% Finally plot difference between the response latency for the maximum
% detected stim pause duration and that of the minimum detected stim pause
% duration
axes(ha(3))
plot(1:l,DRL(1,:),'Marker','o','MarkerSize',7,'MarkerFaceColor','k','Color','k','LineStyle','none')
hold on
plot([0 10],zeros(1,2),'k--','LineWidth',2)
hd(1) = plot(NaN,NaN,'Marker','o','MarkerSize',7,'MarkerFaceColor','k','Color','k','LineStyle','none');
ld = legend(hd,'Difference in Response Latencies (\mu)','Location','southwest','FontSize',12);
ld.ItemTokenSize(1) = 30;
ylabel({strcat('(Response Latency_{MaxStimPause} ',char(8211)),'Response Latency_{MinStimPause}) (s)'},'FontSize',14)
grid on
grid minor
% ylim([-0.4 0.4])

xlabel('Subject','FontSize',14);
xticklabels({'',patients{1}(1:6),patients{2}(1:6),patients{3}(1:6),patients{4}(1:6),patients{5}(1:6),patients{6}(1:6),patients{7}(1:6)})
xtickangle(0)
set(ha(1:2),'XTickLabel',[])
xlim([0.75 l+0.25]);