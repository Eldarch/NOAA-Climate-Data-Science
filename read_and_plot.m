function out = read_and_plot()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File paths
localDir = 'D:\RECCS';    % location of matlab script and data
addpath([localDir,'/colormaps'])    % location of colormap for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------------------------------------------------------------------------------------------------------
% Input flags
% ------------------------------------------------------------------------------------------------------------------------
readVar='mslp';      % Variable to be read in. Current choices are: mean sea level pressure (mslp)
year=1997;

% ------------------------------------------------------------------------------------------------------------------------
% Read in data
% ------------------------------------------------------------------------------------------------------------------------
% Sea-level pressure data
data=ncread([localDir,'/data/',readVar,'.JRA-55.daily.grid2-intp2.',num2str(year),'.nc'],readVar);   

% ------------------------------------------------------------------------------------------------------------------------
% Explore the data and create a simple plot
% ------------------------------------------------------------------------------------------------------------------------
% Print the dimensions of the data to the command line
%size(data)



%180 corresponds to longitude
%91 corresponds to latitude
%365 corresponds to time

% Create a simply plot, which requires a grid
% Longitude grid
lons=[0:2:358]';
% Latitude grid
lats=[-90:2:90];
t=1;    % Day to plot

numbCont=20;                                      % number of contours to use when plotting 
contInt=(max(data(:))-min(data(:)))/numbCont;     % contour interval
plotScale=[min(data(:)):contInt:max(data(:))];     % plot scale (check to see that this matches the colorbar on the plot)

figure(1)
contourf(lons,lats,squeeze(data(:,:,t))',plotScale)
colormap(bluetored(numbCont,1));
caxis([min(plotScale) max(plotScale)])
colorbar
title(['Mean Sea Level Pressure for 1997'])
xlabel('Longitude')
ylabel('Latitude')


%squeeze takes out all of the dimensions with length 1
%mean sea level pressure was plotted at different coordinates

% TO DO: 
% ------------------
% - Instead of creating lat and lon grid, read it in from the original data file
% Look at the metadata of the netCDF file again (recall: run 'ncdisp(['name of netCDF file'])' from the Matlab command line). Look at the names of the variables 
% inside, they are: mslp, latitude, longitude, and initial_time0_hours. In a similar manner to to what is done above, use the ncread command and read in the latitude and 
% longitude data. It should be something like the following:
 filename='mslp.JRA-55.daily.grid2-intp2.1997.nc';
 ncdisp(filename);
 lons_newGrid=ncread(filename,'lon');
 lats_newGrid=ncread(filename,'lat');

% - Verify that it matches the grid you created above by creating new plot with two panels (use the Matlab 'subplot' function). 
% Call this 'figure (2)'. In the first panel, plot the first time step of 'data' using the original lon/lat grid. In the second panel, plot the same thing, but
% use your new grid. 
% - Save this new figure as both a '.fig' file and as a '.jpg' file
figure(2)
contourf(lons_newGrid,lats_newGrid,squeeze(data(:,:,t))',plotScale)
colormap(bluetored(numbCont,1));
caxis([min(plotScale) max(plotScale)])
colorbar
title(['Mean Sea Level Pressure for 1997'])
xlabel('Longitude')
ylabel('Latitude')
saveas(figure(2),'MeanSeaLevelPressureFor1997.jpg');
saveas(figure(2),'MeanSeaLevelPressureFor1997.fig');
% ------------------------------------------------------------------------------------------------------------------------
% Create time averages via 'for loop'
% ------------------------------------------------------------------------------------------------------------------------
% Find the dimensions of the data
[Ly,Lx,Lt]=size(data);
% First do it for the whole time period
tAvg1=zeros(Ly,Lx);
for i=1:Ly
    for j=1:Lx
        for t=1:Lt
            tAvg1(i,j)=data(i,j,t)+tAvg1(i,j);
        end
    end
end
tAvg1=tAvg1/Lt;
size(tAvg1);


% Plot the result
numbCont=10;                                        % number of contours to use when plotting 
contInt=(max(tAvg1(:))-min(tAvg1(:)))/numbCont;     % contour interval
plotScale=[min(tAvg1(:)):contInt:max(tAvg1(:))];     % plot scale (check to see that this matches the colorbar on the plot)


% TO DO: 
% ------------------
% -create a similar time average plot, but do it using the Matlab function 'mean'
% - plot the result as figure(4) and compare it to figure(3), they should match exactly
% - add axis labels and a title to your plot (look through Matlab help to figure out how to do this)
% - save the two figures as '.fig' file and as a '.jpg' files
size(tAvg1)
size(lons)
size(lats)
figure(3)
contourf(lons,lats,tAvg1',plotScale)
colormap(bluetored(numbCont,1));
caxis([min(plotScale) max(plotScale)])
colorbar
title(['Time Mean of MSLP for 1997'])
xlabel('Longitude')
ylabel('Latitude')
saveas(figure(3),'TimeMeanOfMSLPFor1997.jpg');
saveas(figure(3),'TimeMeanOfMSLPFor1997.fig');

data_mean = squeeze(mean(data(:,:,:),3));

%size(data_mean)

figure(4)
contourf(lons,lats,data_mean',plotScale)
colormap(bluetored(numbCont,1));
caxis([min(plotScale) max(plotScale)])
title(['Time Mean of MSLP for 1997'])
xlabel('Longitude')
ylabel('Latitude')
colorbar
clearvars data_mean



% To begin, read in the time variable
time = ncread([localDir,'/data/',readVar,'.JRA-55.daily.grid2-intp2.',num2str(year),'.nc'],'initial_time0_hours');

% To understand what the numbers ('serial date numbers') in the 'time' variable represent, begin by looking at the metadata of the netCDF file again (recall: run 'ncdisp(['name of netCDF file'])'
% QUESTION: What are the 'units' of time in this file?

%Hours since Jan 1st, Year 0

% Matlab uses the 'ConvertFrom' argument in the 'datetime' function to covert 'serial date numbers' to dates. However, it does not have a routine that
% explicitly translates from the 'hours since 1800-01-01 00:00' convention used in these netCDF files, so a conversion to days and an offset needs to be added
% to convert the serial data number so that the 'datenum' argument in the Matlab 'datetime' function can be used (matlab uses number of days since 0-Jan-0000)
% Read about the Matlab 'datetime' function here (https://www.mathworks.com/help/matlab/ref/datetime.html)

% Convert to hours to days
time = time/24;
time = time+657438;

time_dates = datetime(time,'ConvertFrom','datenum');
% Print the dates to the command line
%display(time_dates)

% Calculate and plot time averages of the data for specific dates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First create two dates (notice I am using a different functionality of the matlab 'datetime' function here compared to how it was used above)
t1 = datetime(1997,6,5);
t2 = datetime(1997,7,15);

% Find what index these dates correspond to in the 'time' variable. This requires converting the dates to serial dates using the matlab 'datenum' function
x1 = find(time==datenum(t1));
x2 = find(time==datenum(t2));

% It is a good idea to clear 'dummy' variables when you are done using them (for memory reasons and because in big sections of code you may want to reuse those dummies)
clearvars t1 t2 

% Verify that the dates are what you think there are
%display(time_dates(x1))
%isplay(time_dates(x2))

% Plot the two dates and add a title to the plots
numbCont=20;                                      % number of contours to use when plotting 
contInt=(max(data(:))-min(data(:)))/numbCont;     % contour interval
plotScale=[min(data(:)):contInt:max(data(:))];     % plot scale (check to see that this matches the colorbar on the plot)

figure(5)
subplot(2,1,1)
contourf(lons,lats,squeeze(data(:,:,x1)'),plotScale)
colormap(bluetored(numbCont,1));
caxis([min(plotScale) max(plotScale)])
title(['Mean Sea Level Pressure Plot for ',datestr(time_dates(x1))])
xlabel('Longitude')
ylabel('Latitude')
subplot(2,1,2)
contourf(lons,lats,squeeze(data(:,:,x2)'),plotScale)
colormap(bluetored(numbCont,1));
caxis([min(plotScale) max(plotScale)])
title(['Mean Sea Level Pressure Plot for ',datestr(time_dates(x2))])
xlabel('Longitude')
ylabel('Latitude')
saveas(figure(5),'MeanSeaLevelPressurePlotForJune5AndJuly15_1997.jpg');
saveas(figure(5),'MeanSeaLevelPressurePlotForJune5AndJuly15_1997.fig');

% Plot the time average between the two dates
data_mean = squeeze(mean(data(:,:,x1:x2),3));

%size(data_mean)

figure(6)
contourf(lons,lats,data_mean',plotScale)
colormap(bluetored(numbCont,1));
caxis([min(plotScale) max(plotScale)])
title(['Time Mean of MSLP for ',datestr(time_dates(x1)),' to ',datestr(time_dates(x2))])
xlabel('Longitude')
ylabel('Latitude')
clearvars data_mean
saveas(figure(6),'TimeMeanOfMSLPForJune5ToJuly15_1997.jpg');
saveas(figure(6),'TimeMeanOfMSLPForJune5ToJuly15_1997.fig');


% TO DO:
% ------------------
% - Create a time average for two different time periods: (1) the month (January), and 
%   (2) the final six months of the year (July-December)
% - do the above step using both the 'for loop' method and the 'mean' method
% - plot the results as new figures

% First create two dates (notice I am using a different functionality of the matlab 'datetime' function here compared to how it was used above)
t1 = datetime(1997,1,1);
t2 = datetime(1997,1,31);
t3 = datetime(1997,7,1);
t4 = datetime(1997,12,31);

% Find what index these dates correspond to in the 'time' variable. This requires converting the dates to serial dates using the matlab 'datenum' function
x1 = find(time==datenum(t1));
x2 = find(time==datenum(t2));
x3 = find(time==datenum(t3));
x4 = find(time==datenum(t4));

% It is a good idea to clear 'dummy' variables when you are done using them (for memory reasons and because in big sections of code you may want to reuse those dummies)
clearvars t1 t2 t3 t4

% Verify that the dates are what you think there are
%display(time_dates(x1))
%display(time_dates(x2))
%display(time_dates(x3))
%display(time_dates(x4))

% Plot the two dates and add a title to the plots
numbCont=20;                                      % number of contours to use when plotting 
contInt=(max(data(:))-min(data(:)))/numbCont;     % contour interval
plotScale=[min(data(:)):contInt:max(data(:))];     % plot scale (check to see that this matches the colorbar on the plot)

% Plot the time average between the two dates
data_mean = squeeze(mean(data(:,:,x1:x2),3));

%size(data_mean)

figure(7)
contourf(lons,lats,data_mean',plotScale)
colormap(bluetored(numbCont,1));
caxis([min(plotScale) max(plotScale)])
title(['Time mean of MSLP for ',datestr(time_dates(x1)),' to ',datestr(time_dates(x2))])
xlabel('Longitude')
ylabel('Latitude')
clearvars data_mean
saveas(figure(7),'TimeMeanOfMSLPForJanuary1997.jpg');
saveas(figure(7),'TimeMeanOfMSLPForJanuary1997.fig');

% Plot the time average between the two dates
data_mean = squeeze(mean(data(:,:,x3:x4),3));

size(data_mean)

figure(8)
contourf(lons,lats,data_mean',plotScale)
colormap(bluetored(numbCont,1));
caxis([min(plotScale) max(plotScale)])
title(['Time mean of MSLP for ',datestr(time_dates(x3)),' to ',datestr(time_dates(x4))])
xlabel('Longitude')
ylabel('Latitude')
clearvars data_mean
saveas(figure(8),'TimeMeanOfMSLPFfromJulyToDecember1997.jpg');
saveas(figure(8),'TimeMeanOfMSLPFfromJulyToDecember1997.fig');


% TO DO:
% ------------------
% - open multiple years of data and create January average for multiple years (i.e., average the January average for 1997-1999)
% - plot the results
% - save the data that you plot for this step as a netCDF file (see ncreate function in the Matlab help)

years=[1997:1:1999];
curYear=1997;
lats=ncread('mslp.JRA-55.daily.grid2-intp2.1997.nc','lat');
lons=ncread('mslp.JRA-55.daily.grid2-intp2.1997.nc','lon');
Ny=length(lats);
Nx=length(lons);
avgMSLP=zeros(Nx,Ny,31);


for t=1:length(years)
    %lats=ncread([localDir,'/data/',readVar,'.JRA-55.daily.grid2-intp2.',num2str(curYear),'.nc'],'lat');
    %lons=ncread([localDir,'/data/',readVar,'.JRA-55.daily.grid2-intp2.',num2str(curYear),'.nc'],'lon');
    JanData=ncread([localDir,'/data/',readVar,'.JRA-55.daily.grid2-intp2.',num2str(curYear),'.nc'],'mslp');
    JanData=squeeze(JanData(:,:,1:31)); %Getting the data from just January
    avgMSLP=avgMSLP+JanData; %adding all the January data together
    curYear=curYear+1; %increment the year
    clear JanData; %clear the temp variable
end
%avgMSLP
%length(years)
avgMSLP=avgMSLP/length(years); %Have an average MSLP for every day in January
%avgMSLP
%size(avgMSLP)
avgMSLP=squeeze(mean(avgMSLP,3));
%size(avgMSLP)

figure(9)
contourf(lons,lats,avgMSLP',plotScale)
colormap(bluetored(numbCont,1));
caxis([min(plotScale) max(plotScale)])
title(['January Climatology from ',num2str(years(1)),' to ',num2str(curYear - 1)])
xlabel('Longitude')
ylabel('Latitude')
delete("avgMSLP.nc");
nccreate('avgMSLP.nc','avgMSLP');
saveas(figure(9),'JanuaryClimatology1997To1999.jpg');
saveas(figure(9),'JanuaryClimatology1997To1999.fig');


% TO DO:
% ------------------
% - Create a January climatology and three January anomalies
% This proceeds in two steps:
% STEP 1: Create the January climatology for 1997-1999, which you already did above
% STEP 2: Create and plot three January anomalies, one each for 1997, 1998, 1999. To do this, subtract the January climatology from the
% the mean for the specific year 
% Read more about what an anomaly is here, which talks about anomalies in terms of temperature, but the concept is more general https://en.wikipedia.org/wiki/Temperature_anomaly

curYear=1997;
for t=1:length(years)
    JanData=ncread([localDir,'/data/',readVar,'.JRA-55.daily.grid2-intp2.',num2str(curYear),'.nc'],'mslp');
    JanData=squeeze(JanData(:,:,1:31)); %Getting the data from just January
    JanData=squeeze(mean(JanData,3));%averaging the data in January for 1 average MSLP for the month
    curYear=curYear+1; %increment the year

    figure(9+t)
    contourf(lons,lats,JanData',plotScale)
    colormap(bluetored(numbCont,1));
    caxis([min(plotScale) max(plotScale)])
    title(['January Anomaly for ',num2str(years(t))])
    xlabel('Longitude')
    ylabel('Latitude')
    saveas(figure(9+t),'JanuaryAnomaly1999.jpg');
    saveas(figure(9+t),'JanuaryAnomaly1999.fig');
    clear JanData; %clear the temp variable
end


end