%% Program Presets
clear all;
close all;
clc;

% To use the program:
%1. select LMA_filename of the day you want to study the sources
%2. Set the "examined_day" variable as the same day.
%3. Set timestep
%4. Set plotting options: histograms/plotting (sources in
%time)/sources+events in time
%5.RUN

disp("########## START OF EXECUTION ##########");
%{
sources_data_file = 'C:\Users\icar\Google Drive\LIGHTNING\10minPeriodData\1710181700\171018_1700.txt';
events_data_file = 'C:\Users\icar\Desktop\LIS_HDF_files\Ebre\Ebre_txt\ISS_LIS_20171018_1701_1701_events.txt';
read_workspace_dir = 'C:\Users\icar\Desktop\LIS_HDF_files\Ebre\Ebre_other_info\workspace_40.7_0.7_171018_180713.mat';
corrected_file='C:\Users\icar\Desktop\txt4database\corrected_file.txt';
%}
LIS_total_filename = '/Users/Icar/Google Drive/TFG/LMA_LIS_DATA/ISS_LIS.txt';
LMA_filename = '/Users/Icar/Google Drive/TFG/LMA_LIS_DATA/LMA_FOV_20171018_pass1.txt';
examined_day = datetime(2017,10,18); %used to get the LMA date correctly
 

disp("###Remember to set the examined_day to the one from where the LMA is coming");

timestep = 2e-3; %sec

%{
% addpath('C:\Users\icar\Desktop\LIS_HDF_files\Ebre\Ebre_txt\');
% addpath('C:\Users\icar\Desktop\Data_from_scilab\');

%starttimes/endtimes
%ebre
%[18-Oct-2017 10:32:46    18-Oct-2017 10:34:31]
%[18-Oct-2017 17:01:25 17:01:28]

%barranca
%[18-Jul-03 6:40:14 6:41:18]
%[18-Jun-08 3:56:04 3:56:38]

%lisdata = import_lis(events_data_file);
correcting=0; %1 if the database is used to correct postion
toinput=0; %1 if the user wants to input the start time and endtime of fov. If 0 the program will take the min(listime) and max(endtime). THIS IS BAD
storedinput=1; %the same as above but previously saved in fovtime:
%}

%Plotting Options
histograms =0;
plotting=1;
sources_and_events = 0; %plot sources and events over time. COMPUTATION LOADING



comparing_length_section=0;


%Other Options
save_workspace=0; %if you want to write workspace. SAFEMODE=0
savingcsvfile=0;



disp(' ');
disp('-----------------------------------');
%% Get LMA data
ouput_data = get_LMA_data(LMA_filename);

disp("Reading LMA data...");


sources_data = get_LMA_data(LMA_filename);

lma.flash = sources_data(:,1);
lma.sources_time = examined_day + seconds(sources_data(:,2));
lma.slats = sources_data(:,4);
lma.slons = sources_data(:,5);
lma.salts = sources_data(:,6);
lma.pwr   = sources_data(:,7);%power
lma.LIS_FOV = sources_data(:,8);
lma.noise_fID = sources_data(:,9);

%eliminate those sources that were not under LIS FOV:
indexes = find(lma.LIS_FOV ==0);

                firstsize = length(lma.sources_time);
                secondsize = length(indexes);
                disp(['Eliminating ' 'sources outside the LIS FOV time...']);
                disp([num2str(firstsize-secondsize) ' elements within the FOV']);

lma.flash(indexes) = [];
lma.sources_time(indexes) = [];
lma.slats(indexes) = [];
lma.slons(indexes) = [];
lma.salts(indexes) = [];
lma.pwr(indexes) = [];
lma.LIS_FOV(indexes) =[];
lma.noise_fID(indexes) = [];

%The starttime of that day will be the minum time of the sources that were
%inside the FOV, that day.

starttime = min(lma.sources_time);

endtime = max(lma.sources_time);




%% Get LIS data

   lisdata=import_LIS_total_data(LIS_total_filename);

        disp('Reading LIS data...');

        %LIS data
        day_string = num2str(lisdata(:,1));
        date_day = datetime(day_string,'InputFormat','yyyyMMdd');
        [lis.events_time,I] = sort(date_day+seconds(lisdata(:,2))); 
        
        lis.elats = lisdata(I,3);
        lis.elons = lisdata(I,4);
        lis.erad = lisdata(I,5);
        lis.group = lisdata(I,6);
        lis.glat = lisdata(I,7);
        lis.glon = lisdata(I,8);
        lis.flash = lisdata(I,9);
        lis.flat = lisdata(I,10);
        lis.flon = lisdata(I,11);
        lis.area = lisdata(I,12);
        lis.alat = lisdata(I,13);
        lis.alon = lisdata(I,14);
        lis.aobstime = lisdata(I,15);
        lis.xpixel = lisdata(I,16);
        lis.ypixel = lisdata(I,17);
        lis.bgrad = lisdata(I,18);
        lis.dis2lma = lisdata(I,19);

    

   lis_starttime_address=find(lis.events_time>=starttime,1,'first'); %first source time where LIS had it under its FOV
   lis_endtime_address=find(lis.events_time<=endtime,1,'last');
   
    
    erase_indexes = [find(lis.events_time<starttime); find(lis.events_time>endtime)];
    

        lis.events_time(erase_indexes)=[];
        %lma.seconds_flash(erase_indexes)=[];
        lis.elats(erase_indexes)=[];
        lis.elons(erase_indexes)=[];
        lis.erad(erase_indexes)=[];
        lis.group(erase_indexes)=[];
        lis.glat(erase_indexes)=[];
        %lma.cstp(erase_indexes)=[];
        lis.glon(erase_indexes)=[];
        lis.flash(erase_indexes)=[];
        lis.flat(erase_indexes)=[];
        lis.flon(erase_indexes)=[];
        lis.area(erase_indexes)=[];
        lis.alat(erase_indexes)=[];
        lis.alon(erase_indexes)=[];
        lis.aobstime(erase_indexes)=[];
        lis.xpixel(erase_indexes)=[];
        lis.ypixel(erase_indexes)=[];
        lis.bgrad(erase_indexes)=[];
        lis.dis2lma(erase_indexes)=[];

 disp(['Start time: ']);
 disp(starttime);
 disp(['End time: ']);
 disp(endtime);
 
%% Make chunksinfo structures

%{
%load(read_workspace_dir,'fovinfo');
%[min_starttime_on_area, max_endtime_on_area]=lisboundaries(fovinfo,lma);
%}

%{
%Let's set a function that leaves only the data that should be seen from
%both sensors simultaniously

%Analize all the LMA detections and eliminate those that are outside the
%field of view. Create new lma struct and do the following:
%}


timeperiod = seconds(endtime-starttime);

%ftimestep = timeperiod/(lma.total_nflashes.*10); %an adimensional number. 
                                   %The flashes may not fit in the same chunk but
                                   %but w/ this criteria we got a avalue
                                   %that changes with no of flashes

                                   
%etimestep = timeperiod/(lis.nevents*2); %check if inside a eventstime chunk therees a LMA detection!
etimestep = timestep;
                                   
disp('Dividing time in chunks....');
% Precence
[chunksinfo,chunksdata]=check_detections(etimestep,starttime,endtime,lis,lma);

chunksinfo=e_vs_s_presence(chunksdata,chunksinfo);

chunksinfo.annotation='.times gives the minutes coordinates, from the hour we are regarding, of all chunks that divide the detection period.  Other .chunk_ guive the ADDRESS of the chunks.';


%physical properties of chunks where events and/or sources where detected
[secproperties, scproperties]=chunk_physical_properties(chunksinfo,lis,lma);
disp('Info about the time chunks saved');


chunksinfo

%% Statistics Mean,Medians...
%power
secmaxpwrarray=array_max_power(secproperties);
scmaxpwrarray=array_max_power(scproperties);

semeanmaxpower=mean(secmaxpwrarray);
semedianmaxpower=median(secmaxpwrarray);

smeanmaxpower=mean(scmaxpwrarray);
smedianmaxpower=median(scmaxpwrarray);

secpwrcentroid=array_power_centroid(secproperties);
scpwrcentroid=array_power_centroid(scproperties);

semeanpowercentroid=mean(secpwrcentroid);
semedianpowercentroid=median(secpwrcentroid);

smeanpowercentroid=mean(scpwrcentroid);
smedianpowercentroid=median(scpwrcentroid);


%height
%secmedian, secmen, ... are vectors taht contain the mean height of all the
%sources inside each chunk, one value per chunk.
secmedian=array_h_median(secproperties);
scmedian=array_h_median(scproperties);

secmean=array_h_mean(secproperties);
scmean=array_h_mean(scproperties);

sheightmedianmean=mean(scmedian);
seheightmedianmean=mean(secmedian);

seheightmeanmean=mean(secmean);
sheightmeanmean=mean(scmean);

seheightmeanmedian=median(secmean);
sheightmeanmedian=median(scmean);

%number of sources
%nsources/nesources is the vector that contain information about the number
%of sources in each chunk that contains only/with events sources

[nesources,nsources]=nsourcesxchunks(chunksinfo);
senumsourcesmean=mean(nesources);
senumsourcesmedian=median(nesources);

snumsourcesmean=mean(nsources);
snumsourcesmedian=median(nsources);

%Sum of powers in each chunk
%sesumpower/ssumpower are the vectors that contain the simple sum of all
%the powers in each chunk with_events/only sources
[sesumpower,ssumpower]=sumpowersinchunks(chunksinfo,lma);


sesumpowermean=mean(sesumpower);
ssumpowermean=mean(ssumpower);

sesumpowermedian=median(sesumpower);
ssumpowermedian=median(ssumpower);

if save_workspace==1
disp('saving workspace...');
save([erase(sources_data_file,'.txt') '.mat']);
disp('done');

end





%% Plotting representative values
n=1;
if plotting==1
disp('Plotting over time started...');
close all;
%x axis limits:
[h,m,s]=hms([starttime-seconds(etimestep) endtime+seconds(etimestep)]);
lims=minutes(minutes(m)+seconds(s));

% Power may be related to the detection from LIS

%maximum power
figure(n);
n=n+1;

scatter(minutes(chunksinfo.times(chunksinfo.chunks_w_both)),secmaxpwrarray,'*b');
hold on;
scatter(minutes(chunksinfo.times(chunksinfo.chunks_only_sources)),scmaxpwrarray,'.k');
grid on;
title('Maximum power from the sources registered in each time chunk');
xlabel('Time from the current hour [min.]');
ylabel('Max power [dW]');
legend('Chunk w/ events+sources','Chunk w/ only sources');
%xlim(lims);
colormap(jet);

figure(n);
n=n+1;
scatter(minutes(chunksinfo.times(chunksinfo.chunks_w_both)),secpwrcentroid,'*b');
hold on;
scatter(minutes(chunksinfo.times(chunksinfo.chunks_only_sources)),scpwrcentroid,'.k');
grid on;

title('Height of power centroid in the chunk');
xlabel('Time from the current hour [min.]');
ylabel('Height [m]');
legend('Chunk w/ events+sources','Chunk w/ sources');
%xlim(lims);
colormap(jet);

%.-----------------------------------------------------------------
figure(n);
n=n+1;
[or_times, indexv]=sort(lis.events_time,1);
or_pixels(:,1)=lis.xpixel(indexv);
or_pixels(:,2)=lis.ypixel(indexv);

or_times=datenum(or_times);
scatter(or_pixels(:,1),or_pixels(:,2),1,or_times,'.');
axis equal;
grid on;
axis([0 128 0 128]);
title('Pixels excited on the CCD. Colored by time');
colormap(jet);

%-------------------------------------------------------------------

disp('Continuing plotting...');
% Comparing Altitudes in detection
%Median of the height
figure(n);
n=n+1;



scatter(minutes(chunksinfo.times(chunksinfo.chunks_w_both)),secmedian,'*b');
hold on;
scatter(minutes(chunksinfo.times(chunksinfo.chunks_only_sources)),scmedian,'.k');
grid on;
title('Median of sources height in each time chunk');
xlabel('Time from the current hour [min.]');
ylabel('Height [m]');
legend('Chunk w/ events+sources','Chunk w/ sources');
xlim(lims);
colormap(jet);



figure(n);
n=n+1;

scatter(minutes(chunksinfo.times(chunksinfo.chunks_w_both)),secmean,'*b');
hold on;
scatter(minutes(chunksinfo.times(chunksinfo.chunks_only_sources)),scmean,'.k');
grid on;
title('Mean of sources height in each time chunk');
xlabel('Time from the current hour [min.]');
ylabel('Height [m]');
legend('Chunk w/ events+sources','Chunk w/ sources');
%xlim(lims);

colormap(jet);

end

if sources_and_events ==1
disp("Plotting events & sources over time...");
% plot events with sources
figure(n);
n=n+1;
k=0;
grid on;
hold on;

    [C,ia,ic]=unique(lis.flash);
    colorvalues(1)=5;
    for j=2:length(C)
        
        colorvalues(1,j)=colorvalues(j-1)+30;
        
    end
    
    colorvector=colorvalues(ic)';


    for i=1:length(chunksinfo.events)
        
        local_addresses=chunksinfo.events{i};
        if ~isempty(local_addresses)                  
            scatter(lis.events_time(local_addresses),linspace(500,500,length(local_addresses)),lis.erad(local_addresses),colorvector(local_addresses),'x');
             k=k+length(local_addresses);

        end
    
    
    end

    [C,ia,ic]=unique(lma.flash);
    colorvalues(1)=5;
    for j=2:length(C)
        
        colorvalues(1,j)=colorvalues(j-1)+30;
        
    end
    
    colorvector=colorvalues(ic)';
    
    for i=1:length(chunksinfo.sources)
        
        local_addresses=chunksinfo.sources{i};
        if ~isempty(local_addresses)
            try
            scatter(lma.sources_time(local_addresses),lma.salts(local_addresses),abs(lma.pwr(local_addresses)),colorvector(local_addresses),'+');
            catch
                warning("Unable to plot Events & Sources");
            end
        end
    
    
    end
    ylabel('Height [m]');
    title('Events and sources printed over time');
    legend('x Events. Size==radiance. Color==flash.','+ Sources. Size==power. Color==flash.');

%etimestep=seconds(etimestep);
%timev=(starttime-etimestep/2):etimestep:(endtime+etimestep/2);
%xticks(timev);

 xlim([starttime-seconds(etimestep) endtime+seconds(etimestep)]);
 
end
if savingcsvfile==1
% Write CSV file to open with excel
disp('Writing csv file...');
M=[semeanmaxpower semedianmaxpower smeanmaxpower smedianmaxpower semeanpowercentroid semedianpowercentroid smeanpowercentroid smedianpowercentroid sheightmedianmean seheightmedianmean seheightmeanmean sheightmeanmean seheightmeanmedian sheightmeanmedian senumsourcesmean senumsourcesmedian snumsourcesmean snumsourcesmedian sesumpowermean ssumpowermean sesumpowermedian ssumpowermedian];
csvwrite([erase(sources_data_file,'.txt') '.csv'],M);
disp('Done.');

end

%% Plot histograms
if histograms == 1
  disp('Plotting histograms...');

    
    
 %maximum power   
figure(n);
n=n+1;

histogram(secmaxpwrarray,30);
hold on;
histogram(scmaxpwrarray,30);
title('Bins max power histogram');
xlabel('Max power in the bin [dbW]');
ylabel('Counts');
legend('Sources + events bins','only sources bins');


%Heights    

figure(n);
n=n+1;
histogram(secmedian,30);
hold on;
histogram(scmedian,30);

title('Bins median height histogram');
ylabel('Counts');
xlabel('Height [m]');
legend('sources+events bins','only sources bins');

    
figure(n);
n=n+1;

histogram(secmean,30);
hold on;

histogram(scmean,30);
title('Bins height mean histogram');
ylabel('Counts');
xlabel('Height [m]');
legend('sources+events bins','only sources bins');




%sources densities per bin
figure(n);
n=n+1;
histogram(nesources);
hold on;
histogram(nsources);
title('Density of sources per bin');
ylabel('Counts');
xlabel('Number of sources per bin');
legend('Bins with sources+events','Bins with only sources');

    
figure(n);
n=n+1;
histogram(secpwrcentroid,30);
hold on;
histogram(scpwrcentroid,30);
title('Bins power centroid histogram');
ylabel('Counts');
xlabel('Height [m]');
legend('sources+events bins','only sources bins');


figure(n);
n=n+1;
histogram(sesumpower)
hold on;
histogram(ssumpower)
title('Distributions for the sum of power in each bin');
xlabel('Power [dbW]');
ylabel('Counts');
legend('sources+events bins','only sources bins');

 %Power population distribution
 figure(n);
 n=n+1;
 
formatOut='yyyy-mm-dd HH:MM';
histogram(lma.pwr)
title( ['Power Histogram for ' datestr(starttime,formatOut)]);
xlabel('Power [dbW]');
ylabel('Counts');
 
 figure(n);
 n=n+1;
 
histogram(lis.erad)
title(['Radiance histogram for ' datestr(starttime,formatOut)]);
ylabel('Counts');
xlabel('Radiance [J/(m^2 * sr *m)]');
 
     
    
    
figure(n);
n=n+1;
histogram(lma.salts);
title(['Heights histograph for ' datestr(starttime,formatOut)]);
ylabel('Counts');
xlabel('Height [m]');
disp('plotting ended');
    

end

%% Power of events that had a source at 2500m associated

DISCHARGE_ALT = 2500; %m
tolerance = 50; %m. Tolerance for the search 

[rad_dischargealt] = typical_radiance_at_height(DISCHARGE_ALT,tolerance,lma,lis,chunksinfo);

    disp(['Radiances at discharge height: ' num2str(DISCHARGE_ALT) ...
        'm with tolerance: ' num2str(tolerance) 'm:']);
    disp(['Mean: ' num2str(mean(rad_dischargealt.averages)) '[\mu J/sr/m^2/\mu m]']);
    disp(['Median: ' num2str(median(rad_dischargealt.medians)) '[\mu J/sr/m^2/\mu m]']);

figure(n);
n= n+1;

histogram(rad_dischargealt.medians,30);
hold on;
histogram(rad_dischargealt.averages,30);
title(['Typical values for the radiance of events that had sources associated at ' num2str(DISCHARGE_ALT) ...
    'm with tolerance: ' num2str(tolerance) 'm:']);

xlabel('Radiance [\mu J/sr/m^2/\mu m] ');
ylabel('Time bin counts');

text = {strcat(['Altitude: ' num2str(rad_dischargealt.alt) 'm ']); ...
    strcat(['Tolerance: ' num2str(rad_dischargealt.tol) 'm ']);...
    strcat(['Total n of bins: ' num2str(rad_dischargealt.nbins)])};

dim = [.2 .55 .3 .3];
annotation('textbox',dim,'String',text,'FitBoxToText','on');

legend('Radiance median in the bin','Radiance average in the bin');


%% Comparing lightning length
if comparing_length_section == 1
disp('comparing flash lengthes...');
% GET THE LENGTH OF THE FLASHES. WE HAVE TO ASSIGN TIMINGS TO FLASHES AND
% THEN ASSIGN FLASHES OF LMA TO FLASHES OF LIS, SO WE COMPARE THE SAME
% FLASH. The LIS and LAST

    %get the LIS flashes
    %start times (and positions)
    [flashnst,ilistst]=unique(lis.flash,'first'); %"flash number starts / index list starts]
    %end times (and positions)
    [flashnend,ilistend]=unique(lis.flash,'last'); %"flash number ends / index list ends"
    
    lisflashprops.nflash=flashnst;
    lisflashprops.ilistst=ilistst;
    lisflashprops.ilistend=ilistend;
    lisflashprops.times(:,1)=lis.events_time(ilistst);
    lisflashprops.times(:,2)=lis.events_time(ilistend);
    
    
    %get the LMA flashes
    %start times (and positions)
    [flashnst,ilistst]=unique(lma.flash,'first'); %"flash number starts / index list starts]
    %end times (and positions)
    [flashnend,ilistend]=unique(lma.flash,'last'); %"flash number ends / index list ends"
    
    lmaflashprops.nflash=flashnst;
    lmaflashprops.ilistst=ilistst;
    lmaflashprops.ilistend=ilistend;
    lmaflashprops.times(:,1)=lma.sources_time(ilistst);
    lmaflashprops.times(:,2)=lma.sources_time(ilistend);
    
    
    %there might be flashes that lis did not detect (usually the case).
    %Associate the flashes:
    
    flashesrelated=zeros(size(lisflashprops.nflash,1),4);
    flashesrelated(:,1)=lisflashprops.nflash;
    
    
    %now we have to associate the LIS flashes to LMA flashes
    
    %4 options: 1   start time of LIS inside LMA
               %2   end time of LIS inside LMA (the first two cases alredy
               %comprehend the case of a LIS flash inside a LMA flash
               %3   LMA flash inside LIS flash
    %the starttime of lis flash has to be higher or equal 
    for i=1:size(lisflashprops.times,1)
        
        lmaindex=intersect(find(lmaflashprops.times(:,1)<=lisflashprops.times(i,1)), find(lmaflashprops.times(:,2)>=lisflashprops.times(i,1))); %the start of LIS is comprehended in the LMA flash
        
        if isempty(lmaindex) %in case that it is not the start, but the end of a LIS flash that is comprehended inside a LMA flash
            
            lmaindex=intersect(find(lmaflashprops.times(:,1)<=lisflashprops.times(i,2)), find(lmaflashprops.times(:,2)>=lisflashprops.times(i,2)));
            
            if isempty(lmaindex) %repeat the process for in case that the LMA flash is comprehended inside LIS flash
                
                
                lmaindex=intersect(find(lmaflashprops.times(:,1)>=lisflashprops.times(i,1)), find(lmaflashprops.times(:,2)<=lisflashprops.times(i,2))); %the start of LIS is comprehended in the LMA flash
                
                
            end
            
            
        end
        
        if isempty(lmaindex)
            
            disp('error on comparing flashes times. lma index not assigned');
            disp('It also can be that there are non-simultaneous flashes. (This can be seen in the data');
            
        else
            
            flashesrelated(i,2)=lmaflashprops.nflash(lmaindex);
            
        end
            
        
    end
   
    
    for i=1:size(flashesrelated,1)
        
        if flashesrelated(i,2) ~=0
           
            locallisflash=flashesrelated(i,1);
            locallmaflash=flashesrelated(i,2);
            
            indexlis=find(lisflashprops.nflash==locallisflash);
            
            timeslis=(lisflashprops.times(indexlis,:));
            
            durationlis=second(timeslis(2))-second(timeslis(1));
            
            indexlma=find(lmaflashprops.nflash==locallmaflash);
            
            timeslma=(lmaflashprops.times(indexlma,:));
            
            durationlma=second(timeslma(2))-second(timeslma(1));
            
            
            flashesrelated(i,3)=durationlis;
            flashesrelated(i,4)=durationlma;
            
        end
        
    end
    
    diffdurations=flashesrelated(:,4)-flashesrelated(:,3);
    
    
    figure(n);
    n=n+1;
    histogram(flashesrelated(:,3),30);
    hold on;
    histogram(flashesrelated(:,4),30);
    legend('LIS duration count','LMA duration count');
    ylabel('Counts');
    xlabel('Duration [sec.]');
    
    
    
    
    
    
    
    disp('END at line 660');
    disp(' '); disp(' '); disp(' ');
    %flashesrelated content: /LIS#/ LMA#/ LIS durat./ LMA durat./
    
 


end


%% END
disp("################END OF EXECUTION #####################");
disp(" ");
%% Functions


function [chunksinfo,chunksdata] = check_detections(timestep,starttime,endtime,lis,lma)
timestep=seconds(timestep);


timev=(starttime-timestep/2):timestep:(endtime+timestep/2);
%the way the vector is constructed will assure that the first and last
%event/source is within the time perdiod. This vector has the values of the
%surrounding positions of timechunkns. It has thee end and final time in
%it.

central_chunk_times=starttime:timestep:endtime; %will have the same size as the chunk vector

%now we should assign events and sources to each chunk. We will assign, to
%each chunk, what events and sources index there are. In the first row we
%will put the events and in the second one the sources

%{
chunksdata{2,length(timev)-1}=[];

    for i = 1:lis.nevents %check in what chunks events were detected
        
        
        subind=find(lis.events_time(i)>timev,1,'last');
        
        if ~isempty(subind)
            
            chunksdata{1,subind} = [chunksdata{1,subind}, i];
 
        end
        
    end
    
    for i = 1:lma.nsources %check in what chunks sources were detected
        
        subind=find(lma.sources_time(i)>timev,1,'last'); %notice here the vector is the "timev". There is no >= because in case that the source_time would be equal to timev(end) it would assign to the last chunk, indexed throught the left time limit of the chunk
        
           if ~isempty(subind) %if it does find something
            
            chunksdata{2,subind} = [chunksdata{2,subind}, i];
 
           end
        
    end
%}
    %for each time chunk check if event and/or source was recorded
    
    %lets store only the minutes and seconds of the chunk since the time is
    %obvious
    [h,m,s]=hms(central_chunk_times);
    times=minutes(m)+seconds(s);
    
    
lis_addresses=knnsearch(seconds(central_chunk_times-datetime(1993,1,1))',seconds(lis.events_time-datetime(1993,1,1)));
lma_addresses=knnsearch(seconds(central_chunk_times-datetime(1993,1,1))',seconds(lma.sources_time-datetime(1993,1,1)));%to make the knn search we need numeric arrays.
%what we do is to set those datetime vectors as duration vectors (with the
%same time reference) to afterwards pass them to numeric seconds with
%seconds() function

chunksdata{2,length(central_chunk_times)}=[];

for i=1:length(lis_addresses)
    
    chunksdata{1,lis_addresses(i)}=[chunksdata{1,lis_addresses(i)} i];
end

for i=1:length(lma_addresses)
    
    chunksdata{2,lma_addresses(i)}=[chunksdata{2,lma_addresses(i)} i];
end

chunksinfo.times=times;
chunksinfo.events=chunksdata(1,:);
chunksinfo.sources=chunksdata(2,:);    
chunksinfo.timestep=timestep;
chunksinfo.timeperiod=[starttime, endtime];


end

function chunksinfo=e_vs_s_presence(chunksdata,chunksinfo)

%this function only checks if each chunk is full of sources/events. It is,
%therefore, realted to its presence on the chunks.

    chunks_w_events=find(~cellfun(@isempty,chunksdata(1,:)));%logical vectors. 1==has events/sources. 0==empty
    chunks_w_sources=find(~cellfun(@isempty,chunksdata(2,:)));
    
    chunks_w_both=intersect(chunks_w_events,chunks_w_sources);
    
    chunks_only_events=setdiff(chunks_w_events,chunks_w_sources);
    chunks_only_sources=setdiff(chunks_w_sources,chunks_w_events);
    
    empty_chunks=setdiff(1:1:length(chunksdata),[chunks_w_events chunks_w_sources]);
    
    chunksinfo.chunks_w_events=chunks_w_events;
    chunksinfo.chunks_w_sources=chunks_w_sources;
    chunksinfo.chunks_w_both=chunks_w_both;
    chunksinfo.chunks_only_events=chunks_only_events;
    chunksinfo.chunks_only_sources=chunks_only_sources;
    chunksinfo.empty_chunks=empty_chunks;
    
  
end

function [secproperties, scproperties]=chunk_physical_properties(chunksinfo,lis,lma)

secproperties=[];
scproperties=[];


s_and_e=chunksinfo.sources(chunksinfo.chunks_w_both); %sources address that appeared alongside with events

%see what's the maximum height of each chunk and its median. (Half sources
%will be up that value and half will be under that value).

%lets store some physical properties of each chunk where sources and events
%where detected

%secproperties stands for "sources with events chunk properties"
%each "i" relates to a chunk of time!
    for i=1:length(s_and_e)
        
        local_addresses=s_and_e{1,i};
        local_heights=lma.salts(local_addresses);
        secproperties(i).mean_s_h=mean(local_heights);
        secproperties(i).median_s_h=median(local_heights);
        
        local_pwr=lma.pwr(local_addresses);
        secproperties(i).max_pwr=max(local_pwr);
        
        %influence of the power on the detection
        local_pwr = 10.^(local_pwr/10); % temporal conversion to linear units to weight the centroid
        secproperties(i).pwr_centroid=sum(local_pwr.*local_heights)/sum(local_pwr);
        
        %properties coming from LIS detection
        
        
    end

os=chunksinfo.sources(chunksinfo.chunks_only_sources); 

%store the same physical variables of chunks where only sources were
%detected (and so events should have been detected also)

%secproperties stands for "sources chunk properties"

    for i=1:length(os)
        
        local_addresses=os{1,i};
        local_heights=lma.salts(local_addresses);
        scproperties(i).mean_s_h=mean(local_heights);
        scproperties(i).median_s_h=median(local_heights);
        
        local_pwr=lma.pwr(local_addresses);
        scproperties(i).max_pwr=max(local_pwr);
        
        %influence of the power on the detection
        
        local_pwr = 10.^(local_pwr/10); % temporal conversion to linear units to weight the centroid
        scproperties(i).pwr_centroid=sum(local_pwr.*local_heights)/sum(local_pwr);
       
    end


   disp('Physical Poperties of chunks were events and/or sources were detected done');
end


function property=array_power_centroid(struct)

property=zeros(1,length(struct));

maximum = 0;
indexmax = 0;
    for i=1:length(struct)
        
        property(i)=struct(i).pwr_centroid;
        
%         if property(i)>maximum
%             maximum = property(i);
%             indexmax = i;
%         end
    end
% disp(maximum);
% disp(indexmax);
    
    
    
end

function property=array_max_power(struct)
    
    property=zeros(1,length(struct));

    for i=1:length(struct)
        
        property(i)=struct(i).max_pwr; %the array property will have the maximum power of each chunk

    end

end

function property=array_h_mean(struct)
    
    property=zeros(1,length(struct));

    for i=1:length(struct)
        
        property(i)=struct(i).mean_s_h;

    end

end

function property=array_h_median(struct)
    
    property=zeros(1,length(struct));

    for i=1:length(struct)
        
        property(i)=struct(i).median_s_h;

    end

end

function write_txt_files(lma,write_dir)
    
    
addpath(write_dir);

%processing sources times to let them in a OK format for Paulino

fullfilename=fullfile(write_dir,'sources2database.txt');

disp('Opening .txt writed file...');
stringtimes=datestr(lma.sources_time,'yyyy-mm-dd;HH:MM:SS');
fileID=fopen(fullfilename,'w');

if fileID == -1
    
    disp('Could not open file.');
else
    
    disp('File opened. Starting to print');

    for i=1:length(lma.sources_time)
        
        fprintf(fileID, '%s;%f;%f\r\n',stringtimes(i,:),lma.slats(i),lma.slons(i));
    
    
    end

fclose(fileID);

disp('Printing ended and file closed.');
end
end

function [nesources,nsources]=nsourcesxchunks(chunksinfo)


global_addresses=chunksinfo.sources(chunksinfo.chunks_w_both);

nesources=zeros(1,length(global_addresses));

    for i=1:length(global_addresses)
        
        local_addresses=global_addresses{i};
        
        nesources(i)=length(local_addresses);
        
    end

global_addresses=chunksinfo.sources(chunksinfo.chunks_only_sources);

nsources=zeros(1,length(global_addresses));

   for i=1:length(global_addresses)
        
        local_addresses=global_addresses{i};
        
        nsources(i)=length(local_addresses);
        
    end

end

function [sesumpower,ssumpower]=sumpowersinchunks(chunksinfo,lma)



global_addresses=chunksinfo.sources(chunksinfo.chunks_w_both);

sesumpower=zeros(1,length(global_addresses));

    for i=1:length(global_addresses)
        local_addresses=global_addresses{i};
        
        pwr_sumable = 10.^(lma.pwr(local_addresses)/10); %from dB to W
        
        sesumpower(i)=10*log10((sum(pwr_sumable))); %from W to dB
        
    end
global_addresses=chunksinfo.sources(chunksinfo.chunks_only_sources);

ssumpower=zeros(1,length(global_addresses));

    for i=1:length(global_addresses)
        local_addresses=global_addresses{i};
        
         pwr_sumable = 10.^(lma.pwr(local_addresses)/10); %from dB to W
        
        ssumpower(i)=10*log10((sum(pwr_sumable))); %from W to dB
        
    end
end


%Reading data
function sources_data = import_sourcesdata(filename)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   SOURCES_DATA = IMPORTFILE(FILENAME) Reads data from text file FILENAME
%   for the default selection.
%
%   SOURCES_DATA = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   sources_data = importfile('2017_10_06_10_30.txt', 2, 4856);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/07/27 12:20:39

% Initialize variables.
delimiter = ' ';
if nargin<= 2
    startRow = 2;
    endRow = inf;
end

% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%*s%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block = 2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col = 1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

% Close the text file.
fclose(fileID);

% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

% Create output variable
sources_data = [dataArray{1:end-1}];
end
function sources_header = import_sources_header(filename)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   SOURCES_DATA = IMPORTFILE(FILENAME) Reads data from text file FILENAME
%   for the default selection.
%
%   SOURCES_DATA = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   sources_data = importfile('2017_10_06_10_30.txt', 1, 1);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/07/27 10:45:10

% Initialize variables.
delimiter = ' ';
if nargin<= 2
    startRow = 1;
    endRow = 1;
end

% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block = 2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col = 1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

% Close the text file.
fclose(fileID);

% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

% Create output variable
sources_header = [dataArray{1:end-1}];
end
function lisdata = import_lis(filename)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   ISSLIS2017101916101610EVENTS = IMPORTFILE(FILENAME) Reads data from
%   text file FILENAME for the default selection.
%
%   ISSLIS2017101916101610EVENTS = IMPORTFILE(FILENAME, STARTROW, ENDROW)
%   Reads data from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   ISSLIS2017101916101610events = importfile('ISS_LIS_20171019_1610_1610_events.txt', 2, 7);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/07/27 12:09:22

% Initialize variables.
delimiter = ' ';
if nargin<= 2
    startRow = 2;
    endRow = inf;
end

% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block = 2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col = 1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

% Close the text file.
fclose(fileID);

% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

% Create output variable
lisdata = [dataArray{1:end-1}];
end


function [min_starttime_on_area, max_endtime_on_area]=lisboundaries(fovinfo,lma)

%which file are we lookin in?

k=1; %i know "manually that in this case is the first interesting file
    coordinates=fovinfo.fov_coordinates;

    fovlats=coordinates(1,:);
    fovlons=coordinates(2,:);

    minlat=lma.minlat; 
    minlon=lma.minlon;
    maxlat=lma.maxlat;
    maxlon=lma.maxlon;

    disp('Checking which points are inside the interesting area...');
    inside_range_index=intersect(intersect(find(fovlats>minlat),find(fovlats<maxlat)),intersect(find(fovlons>minlon),find(fovlons<maxlon)));

    %check at when the lis will start to see and leave the area. This is not
    %exact due to the fact that maybe the centroids are outside but the fov
    %cell has some part inside or viceversa.

    interestingfovtimes=[fovinfo(k).fovstart(inside_range_index); fovinfo(k).fovend(insiderange_index)];

    min_starttime_on_area=min(interestingfovtimes(1,:)); %minimum time where the lis is sensing the area
    max_endtime_on_area=max(interestingfovtimes(2,:)); %last time when LIS saw info on the area
    
    min_starttime_on_area=datetime(1993,1,1)+seconds(min_starttime_on_area);
    max_endtime_on_area=datetime(1993,1,1)+seconds(max_endtime_on_area);
    
    %we are assuming that LIS will record until the last moment on the
    %interesting area. Truly, during this last time the LIS only would be
    %able to see the little, last, only "franja" of the area.
end




function [date,time,lat,lon] = importcorrected(filename)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [DATE1,TIME,LAT,LON] = IMPORTFILE(FILENAME) Reads data from text file
%   FILENAME for the default selection.
%
%   [DATE1,TIME,LAT,LON] = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads
%   data from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   [date1,time,lat,lon] = importfile('paulino.txt.txt',1, 1713);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/07/31 13:20:32

% Initialize variables.
delimiter = ';';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

% Format for each line of text:
%   column1: datetimes (%{yyyy-MM-dd}D)
%	column2: datetimes (%{HH:mm:ss}D)
%   column3: double (%f)
%	column4: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%{yyyy-MM-dd}D%{HH:mm:ss}D%f%f%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

% Close the text file.
fclose(fileID);

% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

% Allocate imported array to column variable names
date = dataArray{:, 1};
time = dataArray{:, 2};
lat = dataArray{:, 3};
lon = dataArray{:, 4};

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% date1=datenum(date1);
% time=datenum(time);
end  


function [dir_list]=check_for_folders(read_dir,dir_list)
disp('Looking for folders with HDF4 files inside');
addpath(read_dir); %current reading directory
folderinfoprev=dir(read_dir);
folderinfo=folderinfoprev(~ismember({folderinfoprev.name},{'.','..','.DS_Store'})); %.DS_Store is a metadata file created by iOS environment

aretherefolders=cell2mat({folderinfo.isdir});

if any(aretherefolders)==true %check if there are folders inside the folder
    
    namesarray={folderinfo(aretherefolders).name};
    %gives a cell array but in char
    
    %foldernames=convertCharsToStrings(namesarray); %lets convert it to string  

    for i=1:length(namesarray)
   % new_read_dir(i)=fullfile(read_dir,foldernames(i));
    new_read_dir=(fullfile(read_dir,(namesarray{i})));
    
    [dir_list]=check_for_folders(new_read_dir,dir_list); %if the folder contains more folders, re-check
    
    end
    
else
    
        dir_list=[dir_list; read_dir]; %if the folder is a file folder save its directory and make it travel through the function

     
end

disp('Read directories will be:');
    for i=1:size(dir_list,1)
        disp(dir_list(i,:));
    end
    disp(' ');
end






%NEW VERSION FUNCTIONS
function ouput_data = get_LMA_data(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   OUPUT_DATA = IMPORTFILE(FILENAME) Reads data from text file FILENAME
%   for the default selection.
%
%   OUPUT_DATA = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   ouput_data = importfile('total_LMA_file_noNoise.txt', 1, 4178);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2019/03/16 10:47:29

%% Initialize variables.
delimiter = '\t';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
ouput_data = [dataArray{1:end-1}];
end
function output_data = import_LIS_total_data(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   OUTPUT_DATA = IMPORTFILE(FILENAME) Reads data from text file FILENAME
%   for the default selection.
%
%   OUTPUT_DATA = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   output_data = importfile('ISS_LIS.txt', 1, 216837);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2019/03/16 11:00:14

%% Initialize variables.
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
%   column19: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
output_data = [dataArray{1:end-1}];
end

function [rad_dischargealt] = typical_radiance_at_height(DISCHARGE_ALT,tolerance,lma,lis,chunksinfo)


    rad_dischargealt.medians = [];
    rad_dischargealt.averages = [];
    counter = 0;
    for i = 1:length(chunksinfo.chunks_w_both)

        chunk_address = chunksinfo.chunks_w_both(i);

        local_sources = chunksinfo.sources{chunk_address};


        if any(lma.salts(local_sources)>DISCHARGE_ALT-tolerance) ...
                && any(lma.salts(local_sources)<DISCHARGE_ALT +tolerance) %The heights of the sources is comprised in the tolerance 

            counter = counter +1;

            local_events = chunksinfo.events{chunk_address}; %access the content of the cell at that address;

            local_rads = lis.erad(local_events);

            median_rads = median(local_rads);
            average_rads = mean(local_rads); 

            rad_dischargealt.medians(counter) = median_rads;
            rad_dischargealt.averages(counter) = average_rads;

        end

        rad_dischargealt.alt = DISCHARGE_ALT;
        rad_dischargealt.tol = tolerance;
        rad_dischargealt.nbins = counter;
    end

end