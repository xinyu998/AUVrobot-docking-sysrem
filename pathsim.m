%MAE 6760 PROJECT
%Path planning particle filter
clear all;close all;
rng('default')

video_name='PF_Gaussian1.mp4';
record_video=0;

rng(100);
%% Define map
myMap=binaryOccupancyMap(8,8,5)
walls=zeros(40,40);
walls(1,:)=1;
walls(end,:)=1;
walls(:,1)=1;
walls(:,end)=1;
walls(10:20,10)=1;
walls(1:10,20)=1;
walls(30:40,20)=1;
walls(20,1:10)=1;
walls(10,10:20)=1;
walls(30,20:40)=1;
walls(10,30:40)=1;
walls(1:10,30)=1;
setOccupancy(myMap,[1 1],walls,"grid")
show(myMap)
xlabel('X [Miles]');
ylabel('Y [Miles]');
title('AUV Workspace')
