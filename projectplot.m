clear all;close all;
rng('default')

video_name='PF_Gaussian1.mp4';
record_video=0;

rng(100);
%
% Two cases to consider: select either baseline or swervy
scenario_type='baseline';

[Uacc,Uomega]=get_controlinputs(scenario_type);

%state vector: x,y 2D position, velocity, heading
%simulate no noise system
[Xnonoise,n,t,dt,nt]=simulate_2Dcar(Uacc,Uomega);
nk=nt;
%
%plots
plot_birdseyeview(Xnonoise,[],[],'Truth: Birds Eye View');
ii_plot=[3 4];
plot_estimator(t,Xnonoise,[],[],ii_plot,'Truth: Velocity/Heading States');

%% Simulation Function calls

function [Uacc,Uomega]=get_controlinputs(scenario_type);
%
%scenario_type='baseline';
%scenario_type='swervy';
%
%   generates clean acceleration and rotational rate control inputs
%
%
%define inputs: acceleration and heading rate
if strcmp(scenario_type,'baseline'),
    Uacc=[ones(1,60) zeros(1,130) -ones(1,60) zeros(1,100)];
    Uomega=[ones(1,50)*pi/2/2 zeros(1,50) ones(1,50)*pi/2/4 zeros(1,50) -ones(1,50)*pi/2/4 zeros(1,50) -ones(1,50)*pi/2/2];
elseif strcmp(scenario_type,'swervy'),
    Uacc=[ones(1,60) zeros(1,230) -ones(1,60)];
    Uomega=[zeros(1,60) ones(1,10)*pi/2/2 -ones(1,20)*pi/2/2 ones(1,20)*pi/2/2 ...
        -ones(1,20)*pi/2/2 ones(1,20)*pi/2/2 -ones(1,10)*pi/2/2 zeros(1,30) ...
        ones(1,10)*pi/2/2 -ones(1,20)*pi/2/2 ones(1,20)*pi/2/2 ...
        -ones(1,20)*pi/2/2 ones(1,20)*pi/2/2 -ones(1,10)*pi/2/2 ...
        zeros(1,60)];  
else,
    return;
end
%
end

function [Xnonoise,n,t,dt,nt]=simulate_2Dcar(Uacc,Uomega);
%
%   simulates a 2D car using acceleration and rotational rate control inputs
%   state vector: x,y 2D position, velocity, heading
%
dt=0.1;
nt=length(Uacc);
t=[0:dt:dt*(nt-1)];
n=4;
Xnonoise=zeros(n,nt);
for k=1:(nt-1),
    Vk=Xnonoise(3,k);
    Tk=Xnonoise(4,k);    
    Xnonoise(:,k+1) = Xnonoise(:,k) +...
        dt*[Vk*cos(Tk);Vk*sin(Tk);Uacc(k);Uomega(k)];   
end
end
%% Plotting functions

function plot_estimator(t,x1,x2,P2,ii_plot,title_name);
% x1 is the true value or reference comparison
% x2,P2 is the estimator state and covariance
% ii_plot: 2x1 vector of which states to plot
%
axis_names={'North (m)','East (m)','Velocity (m/sec)','Heading (rad)','Accel Bias (m/sec^2)','RG Bias (rad/sec)'};
figure;subplot(122);
ii_x1=[];ii_x2=[];ii_P2=[]; %for legend
%
for i=1:length(ii_plot),
    ii=ii_plot(i);
    subplot(1,2,i);
    hold on;
    if ~isempty(x1),
        plot(t,x1(ii,:),'color',[0 0.5 0]);ii_x1=1;
    end  
    if ~isempty(x2),
        plot(t,x2(ii,:),'b-');ii_x2=2;
    end  
    if ~isempty(P2)
        plot(t,x2(ii,:)'-2*sqrt(squeeze(P2(ii,ii,:))),'b:');
        plot(t,x2(ii,:)'+2*sqrt(squeeze(P2(ii,ii,:))),'b:');ii_P2=3;
    end
    hold off
    xlabel('time (sec)');ylabel(axis_names(ii));grid;
    xlim([0 35]);set(gca,'xtick',[0:5:35]);
end
legend_names={'true state','estimate','2\sigma bound'};
legend(legend_names{ii_x1},legend_names{ii_x2},legend_names{ii_P2},'Location','South');
%
sgtitle(title_name);
PrepFigPresentation(gcf);
end

function plot_estimator_error(t,x1,x2,P2,ii_plot,title_name);
% x1 is the true value or reference comparison
% x2,P2 is the estimator state and covariance
% ii_plot: 2x1 vector of which states to plot
%
axis_names={'North (m)','East (m)','Velocity (m/sec)','Heading (rad)','Accel Bias (m/sec^2)','RG Bias (rad/sec)'};
figure;subplot(122);
%
for i=1:length(ii_plot),
    ii=ii_plot(i);
    subplot(1,2,i);
    err=x2(ii,:)-x1(ii,:);
    plot(t,err,'b-');
    hold on;
    if ~isempty(P2)
        plot(t,err'-2*sqrt(squeeze(P2(ii,ii,:))),'b:');
        plot(t,zeros(length(t),1),'r--');
        plot(t,err'+2*sqrt(squeeze(P2(ii,ii,:))),'b:');
    end
    hold off
    xlabel('time (sec)');ylabel(axis_names(ii));grid;
    xlim([0 35]);set(gca,'xtick',[0:5:35]);
end
legend('estimator error','2\sigma bound','zero error','Location','South');
%
sgtitle(title_name);
PrepFigPresentation(gcf);
end

function plot_birdseyeview(x1,x2,P2,title_name);
% x1 is the true value or reference comparison
% x2,P2 is the estimator state and covariance
% 
ii_x1=[];ii_x2=[];ii_P2=[]; %for legend
figure;
if ~isempty(x1),
    plot(x1(1,:),x1(2,:),'color',[0 0.5 0]);ii_x1=1;
end
hold on;
if ~isempty(x2),
    plot(x2(1,:),x2(2,:),'b-');ii_x2=2;
end
if ~isempty(P2),
    iell=[2 10:10:350];
    for i=1:length(iell),
        ii=iell(i);
        [Xe,Ye] = calculateEllipseCov(x2([1 2],ii),P2([1 2],[1 2],ii),3);
        plot(Xe,Ye,'m-');
    end
    ii_P2=3;
end
xlabel('North (m)');ylabel('East (m)');grid;
hold off;
legend_names={'true trajectory','estimated trajectory','3\sigma bound'};
legend(legend_names{ii_x1},legend_names{ii_x2},legend_names{ii_P2},'Location','South')
%
title(title_name);
PrepFigPresentation(gcf);
end

function [Xe,Ye] = calculateEllipseCov(X, P, nsig, steps) 
    %# This functions returns points to draw an ellipse 
    %# 
    %#  @param X     x,y coordinates 
    %#  @param P     covariance matrix 
    %# 
 
    error(nargchk(2, 3, nargin)); 
    if nargin<3, nsig = 1; end 
    if nargin<4, steps = 36; end 
    
    [U,S,V]=svd(P);
    s1=sqrt(S(1,1));s2=sqrt(S(2,2));angle=acos(U(1,1))*180/pi;
    x=X(1);
    y=X(2);

    %scale by nsig
    s1=nsig*s1;
    s2=nsig*s2;

    beta = angle * (pi / 180); 
    sinbeta = sin(beta); 
    cosbeta = cos(beta); 
 
    alpha = linspace(0, 360, steps)' .* (pi / 180); 
    sinalpha = sin(alpha); 
    cosalpha = cos(alpha); 
 
    Xe = x + (s1 * cosalpha * cosbeta - s2 * sinalpha * sinbeta); 
    Ye = y + (s1 * cosalpha * sinbeta + s2 * sinalpha * cosbeta); 
 
end 

function PrepFigPresentation(fignum);
%
% prepares a figure for presentations
%
% Fontsize: 14
% Fontweight: bold
% LineWidth: 2
% 

figure(fignum);
fig_children=get(fignum,'children'); %find all sub-plots

for i=1:length(fig_children),
    
    set(fig_children(i),'FontSize',16);
    set(fig_children(i),'FontWeight','bold');
    
    fig_children_children=get(fig_children(i),'Children');
    set(fig_children_children,'LineWidth',2);
end
end