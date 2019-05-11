%////////////////////////////////////////////
%// LIS HDF4 data processor                //
%//       MATLAB converter                 //
%// Universitat Politècnica de Catalunya   //
%//   icar.fontcuberta@gmail.com           //
%//                                        //
%// Requires: MATLAB R2017b at least       //
%//                                        //
%//                                          //
%////////////////////////////////////////////





clear all;


% #################USER MODIFIABLE VARIABLES ####################
read_dir='C:\Users\icar\Desktop\LIS_HDF_files\Barrancabermeja\Barranca_HDF';
write_dir_scilab='C:\Users\icar\Desktop\LIS_HDF_files\Barrancabermeja\Barranca_txt_Jun_Jul';
write_dir='C:\Users\icar\Desktop\LIS_HDF_files\Barrancabermeja\Barranca_txt_Jun_Jul'; 
read_dir_txtfiles=write_dir; %atm we read the txt files from the same place we store them
urls_file='C:\Users\icar\Google Drive\PRACTIQUES\HDF_reading\reading_txt_files\GHRC_URLs.txt';
corrected_urls_file_write_dir='C:\Users\icar\Desktop\LIS_HDF_files\Barrancabermeja\'; %here will also go the file twith urls from the website names file
relevant_orbits_file='C:\Users\icar\Google Drive\PRACTIQUES\HDF_reading\reading_txt_files\relevant_orbits.txt';
website_filename='C:\Users\icar\Desktop\LIS_HDF_files\Barrancabermeja\filenames_website.txt';
map_folder='C:\Users\icar\Google Drive\PRACTIQUES\LMAvsLIS (sci oscar)\Sci_program\mapfolder';

write_other_info_dir='C:\Users\icar\Desktop\LIS_HDF_files\Ebre\Ebre_other_info\'; %specify this to save the workspace and other info there


%coordinates info
deltebre=[40.7212388 0.7176492];
santamarta=[11.2403547 -74.2110227];
barranca=[7.06878 -73.744418];

%scanning_area specification
centroid=barranca; %LAT/LON (remember, on the plot, this would be y and x)
range=60*sqrt(2); %range in km 
%time interval
starttime=datetime(2018,6,1,0,0,0);
endtime=datetime(2018,7,3,8,0,0);
timerange=[starttime endtime];


%###############################################################

%correction for change in radius. extracted from https://rechneronline.de/earth-radius/


B=centroid(1);
r1= 6378.137;%radius at the equator 
r2=6356.752;%radius at the poles
earth_radius=sqrt(((((r1^2)*cos(B))^2)+((r2^2)*sin(B))^2)/((r1*cos(B))^2+(r2*sin(B))^2)); %radius at your location

ang_range=range/earth_radius*360/(2*pi); %range in degrees


disp('Do you want to: ');
disp('0: Exit Program');
disp('1: Write general event txt files');
disp('2: Write general event txt and plot them');
disp('3: Write event txt files for scilab (TO VERIFY) ');
disp('4: Plot events in interesting time-space from HDF4 files (TO VERIFY)');
disp('5: Plot events in interesting time-space from .txt files ');
disp('6: Correct GHRCs URLs and generate a new URLs txt file');
disp('7: Process website filenames to list of interesting URLs');
disp('8: Check only for interesting files and save the workspace.');
disp(' ');
n=input('Enter option: ');



isok=false;

while isok==false
switch n
    case 0
        interestingfiles=[];
        isok=true;
    case 1
        
            dir_list=[];
            dir_list=check_for_folders(read_dir,dir_list);
            interestingfiles=[]; a=1;
            corruptfiles=[]; b=1;
            
        for i=1:size(dir_list,1)
            
            


            local_read_dir=dir_list(i,:);
            [interestingfiles, corruptfiles,a,b]=select_interestingfiles(local_read_dir,centroid, timerange, ang_range,interestingfiles,corruptfiles,a,b,write_other_info_dir);

        end

        w_txtfiles(interestingfiles,write_dir,centroid,ang_range,timerange,n); %in this function we will calibate the events because we don't have the scilab post-processing program
        isok=true;
    case 2
        
            dir_list=[];
            dir_list=check_for_folders(read_dir,dir_list);
            interestingfiles=[]; a=1;
            corruptfiles=[]; b=1;
        
        for i=1:size(dir_list,1)

            local_read_dir=dir_list(i,:);
            [interestingfiles, corruptfiles,a,b]=select_interestingfiles(local_read_dir,centroid, timerange, ang_range,interestingfiles,corruptfiles,a,b,write_other_info_dir);

        end


         
        w_txtfiles(interestingfiles,write_dir,centroid,ang_range,timerange,n); %in this function we will calibate the events because we don't have the scilab post-processing program
        plot_coastline(map_folder,centroid,ang_range);
        
        isok=true;
    case 3
        
            dir_list=[];
            dir_list=check_for_folders(read_dir,dir_list);
            interestingfiles=[]; a=1;
            corruptfiles=[]; b=1;
        
        for i=1:size(dir_list,1)

             local_read_dir=dir_list(i,:);
             [interestingfiles, corruptfiles,a,b]=select_interestingfiles(local_read_dir,centroid, timerange, ang_range,interestingfiles,corruptfiles,a,b,write_other_info_dir);
           
        end


        
        w_txtfiles_4scilab(interestingfiles,write_dir_scilab);
        isok=true;
    case 4
        
            dir_list=[];
            dir_list=check_for_folders(read_dir,dir_list);
            interestingfiles=[]; a=1;
            corruptfiles=[]; b=1;
        
        for i=1:size(dir_list,1)

            local_read_dir=dir_list(i,:);
             [interestingfiles, corruptfiles,a,b]=select_interestingfiles(local_read_dir,centroid, timerange, ang_range,interestingfiles,corruptfiles,a,b,write_other_info_dir);

        end

            read_hdf4_and_plot(interestingfiles,centroid,ang_range,timerange);
            plot_coastline(map_folder,centroid,ang_range);
        
        isok=true;
        
    case 5 %disp('5: Plot events in interesting time-space from .txt files (TO VERIFY)');

            dir_list=[];
            dir_list=check_for_folders(read_dir_txtfiles,dir_list);
            interestingfiles=[]; a=1;
            corruptfiles=[]; b=1;
            
            for i=1:size(dir_list,1)
                
                local_read_dir=dir_list(i,:);
                
                read_txt_and_plot(local_read_dir,centroid,ang_range,timerange);
                
            end
            plot_coastline(map_folder,centroid,ang_range);
            disp('All events in the interesting sace-time dominum plotted');
          
            isok=true;  
    case 6
        
        isok=true;
        correct_hdf_urls(urls_file,relevant_orbits_file,corrected_urls_file_write_dir);
        %if you have a list of relevant orbits computed with the py program
    case 7
        
        isok=true;
        
        disp('Have you entered a correct space-time domain?');
        disp('Write dir for the txt file?');
        webfilenames2urls(corrected_urls_file_write_dir,website_filename);
        %if you have a file taken manually from the data from LIS website
        
    case 8
        
        isok=true;
        
           dir_list=[];
            dir_list=check_for_folders(read_dir,dir_list);
            interestingfiles=[]; a=1;
            corruptfiles=[]; b=1;
        
        for i=1:size(dir_list,1)

            local_read_dir=dir_list(i,:);
             [interestingfiles, corruptfiles,a,b]=select_interestingfiles(local_read_dir,centroid, timerange, ang_range,interestingfiles,corruptfiles,a,b,write_other_info_dir);

        end
    
    otherwise
        n=input('Bad entry. Rechoice option: ');
        
end

end

disp(' ');
disp('Writed files at:');
disp(['.txt files for scilab: ' write_dir_scilab]);
disp(['.txt files general events: ' write_dir]);
disp('Bye!');


%--------------------------------------------------------------------------
% ------------------------------FUNCTIONS----------------------------------
%--------------------------------------------------------------------------
%Found folders where the files are
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
%Select interesting files from those folders
function [interestingfiles, corruptfiles,a,b]=select_interestingfiles(read_dir,centroid, timerange, ang_range,interestingfiles,corruptfiles,a,b,write_other_info_dir)

%-----------------------------NOTATION------------------------------------%
%k=index of current file. all files inside folder, interesting files or not
%a=index of file in the interesting files subgroup
%b=index for possible cannot-read files inside the folder
%i = flash index inside the HDF file
%-------------------------------------------------------------------------%
disp('Selecting the interesting files from the interesting folders. This might take a while, depending on the nº of HDF4 files.');
addpath(read_dir);
folderinfoprev=dir(read_dir);
folderinfo=folderinfoprev(~ismember({folderinfoprev.name},{'.','..','.DS_Store'})); %.DS_Store is a metadata file created by iOS environment

totalnoffiles=size(folderinfo,1);
for k=1:size(folderinfo,1)
%get the names of all files inside your data directory
filename=folderinfo(k).name;
    
    fileinfo=hdfinfo(filename);
    %try to read, the file may be corrupt/empty
    try
        event_vdata=hdfread(fileinfo.Vgroup.Vgroup.Vgroup.Vdata(4));

    catch
        
        %{
        disp(['The file ' filename ' with index ' num2str(k) ' could not be read.']);
        disp(' ');
        %}
        corruptfiles(b).Filename=filename;
        corruptfiles(b).File_index=k;
        b=b+1;
        
        continue
    end
      viewtime_vdata=hdfread(fileinfo.Vgroup.Vgroup.Vdata(2));
        fovinfo(k).fov_coordinates=viewtime_vdata{1};
        fovinfo(k).fovstart=viewtime_vdata{2};
        fovinfo(k).fovend=viewtime_vdata{3};
        
      
    event.coordinates=event_vdata{3};
    event.tai_time=event_vdata{1};
    
          time=datetime(1993,1,1)+seconds(event.tai_time);
          lats=event.coordinates(1,:);
          lons=event.coordinates(2,:);
          
          ang_distances=sqrt((lats-centroid(1)).^2+(lons-centroid(2)).^2);
          
          space_ind=find(ang_distances<ang_range);%get indexes of space dominum
          time_ind=intersect(find(time>timerange(1)),find(time<timerange(2)));
          %get indexes of time dominum. Intersect gives common values, so
          %here its ok because finds gives indexes and we want to find
          %common indexes
          
          [view_ind]=intersect(space_ind,time_ind); %relevant elements of each file
          %intersect both indexes in order to find positios that meet
          %space-time dominum
          
          if ~isempty(view_ind)  
              
              interestingfiles(a).Filename=filename;
              interestingfiles(a).File_index=k;
        
            a=a+1;
              
          end
              
    
    %{
    
    i=1; 
    eventinsiderange=false; eventinsidetime=false;
    while i<=size(event.coordinates,2) && ~(eventinsiderange && eventinsidetime)
       eventinsiderange=false; eventinsidetime=false;
       
       %check space domain
       coordinates=event.coordinates(:,i);
       c_distance=sqrt((coordinates(1)-centroid(1))^2+(coordinates(2)-centroid(2))^2);

       if c_distance<ang_range
           eventinsiderange=true;
       end
       
       %check time domain
       time=datetime(1993,1,1)+seconds(event.tai_time(i));
       if time<timerange(2) && time>timerange(1)
           eventinsidetime=true;
       end
       
       %if both are positive, this file is interesting, and we dont need to
       %further loop
        if eventinsidetime==true && eventinsiderange==true
        
            interestingfiles(a).Filename=filename;
            interestingfiles(a).File_index=k;
        
            a=a+1;
        
        end
    
       i=i+1;
       
    end
    %}
    disp([num2str(k*100/totalnoffiles) '% of the interesting files checking inside ' read_dir ' is done.']);

end

    %save FOV info from these interesting files
    disp(' ');
    disp('Saving Workspace...');
    save([write_other_info_dir 'workspace_' num2str(centroid(1)) '_' num2str(centroid(2)) '_' datestr(timerange(1),'yymmdd') '_' datestr(timerange(2),'yymmdd') '.mat'] ,'fovinfo', 'interestingfiles', 'corruptfiles' ,'centroid', 'ang_range','timerange');

if ~isempty(interestingfiles)
    disp([num2str(length(interestingfiles)) ' interesting files have been found']);
else
    disp(['No interesting files were found in ' read_dir]);
end
    
end
%Make .txt for scilab files from the interesting files
function w_txtfiles_4scilab(interestingfiles,write_dir)
addpath(write_dir);

 
    for k=1:length(interestingfiles)
    
         filename=interestingfiles(k).Filename;
         fileinfo=hdfinfo(filename);
       
         event_vdata=hdfread(fileinfo.Vgroup.Vgroup.Vgroup.Vdata(4));

         filetextname=[erase(filename,'.hdf'),'_event4scilab','.txt'];
         fullfilename=fullfile(write_dir,filetextname);   
         
        tai_time=event_vdata{1};       %get the instant, for comparison pruposes
        observe_time=event_vdata{2};
        location=event_vdata{3};
        radiance=event_vdata{4};
        footprint=event_vdata{5};
        address=event_vdata{6}+1; %add +1 so it doesn't start w/ 0 adrexx
        parent_address=event_vdata{7}+1;
        x_pixel=event_vdata{8};
        y_pixel=event_vdata{9};
        bg_value=event_vdata{10};
        bg_radiance=event_vdata{11};
        amplitude=event_vdata{12};
        sza_index=event_vdata{13};
        glint_index=event_vdata{14};
        approx_threshold=event_vdata{15};
        alert_flag=event_vdata{16};
        cluster_index=event_vdata{17};
        density_index=event_vdata{18};
        noise_index=event_vdata{19};
        bg_value_flag=event_vdata{20};
        grouping_sequence=event_vdata{21};

        
        fileID=fopen(fullfilename,'w');
         if fileID==-1
             disp('Could not open file!');
             
         else
         fprintf(fileID, 'TAI93_time observe_time latitude longitude radiance footprint address parent_address x_pixel ypixel bg_value bg_radiance amplitude sza_index glint_index approx_threshold alert_flag cluster_index density_index noise_index bg_value_flag grouping_sequence\r\n');
         for i=1:length(tai_time)  
         fprintf(fileID, '%-18.16E %2d %7.3f %7.3f %4.1d %3.1d %4u %3u %3u %3u %4u %3u %3u %3u %3u %3u %1u %2u %2u %3d %1u %6u\r\n',tai_time(i),observe_time(i),location(1,i),location(2,i),radiance(i),footprint(i),address(i),parent_address(i),x_pixel(i),y_pixel(i),bg_value(i),bg_radiance(i),amplitude(i),sza_index(i),glint_index(i),approx_threshold(i),alert_flag(i),cluster_index(i),density_index(i),noise_index(i),bg_value_flag(i),grouping_sequence(i)); 
         end
         
         fclose(fileID);
         disp(['Printing file ' filename ' @ ' write_dir]);
         end
           
    end
             
 

end
function w_txtfiles(interestingfiles,write_dir,centroid, ang_range,timerange,n)
disp('Printing general events .txt files...');


ninterestingfiles=length(interestingfiles); %added constant for % use
for k=1:length(interestingfiles) %We will go through all interestingfiles
                                 %and search for interestingevents
    clear area_vdata flash_vdata group_vdata event_vdata
    
         filename=interestingfiles(k).Filename;
         fileinfo=hdfinfo(filename);
        
         event_vdata=hdfread(fileinfo.Vgroup.Vgroup.Vgroup.Vdata(4));
         event.location=event_vdata{3};
         event.tai_time=event_vdata{1};
         
         area_vdata=hdfread(fileinfo.Vgroup.Vgroup.Vgroup.Vdata(1));
         flash_vdata=hdfread(fileinfo.Vgroup.Vgroup.Vgroup.Vdata(2));
         group_vdata=hdfread(fileinfo.Vgroup.Vgroup.Vgroup.Vdata(3));
         event_vdata=hdfread(fileinfo.Vgroup.Vgroup.Vgroup.Vdata(4));
         
        event.location=event_vdata{3};
        event.tai_time=event_vdata{1};       
        event.observe_time=event_vdata{2};
        event.radiance=event_vdata{4};
        event.address=event_vdata{6}+1; %add +1 so it doesn't start w/ 0 adress
        event.parent_address=event_vdata{7}+1;
        event.bg_radiance=event_vdata{11};
        
        event.x_pixel=event_vdata{8};
        event.y_pixel=event_vdata{9};
        event.bg_rad=event_vdata{11};
        
        group.parent_address=group_vdata{7}+1; %no need to store group.address as they are ordered by number inside each file
        group.location=group_vdata{3};
        
        flash.parent_address=flash_vdata{8}+1;
        flash.location=flash_vdata{4};
        
        area.location=area_vdata{4};
        area.observetime=area_vdata{3};
    
        orbit_vdata=hdfread(fileinfo.Vgroup.Vdata(1));
        orbit_id=orbit_vdata{1};
        
        
    
        %We could size this info for printing and plotting, but better not
        %to change what works
          time=datetime(1993,1,1)+seconds(event.tai_time);
          lats=event.location(1,:);
          lons=event.location(2,:);
          
          ang_distances=sqrt((lats-centroid(1)).^2+(lons-centroid(2)).^2);
          
          space_ind=find(ang_distances<ang_range);%get indexes of space dominum
          time_ind=intersect(find(time>timerange(1)),find(time<timerange(2)));
          %get indexes of time dominum. Intersect gives common values, so
          %here its ok because finds gives indexes and we want to find
          %common indexes
          
          [view_ind]=intersect(space_ind,time_ind); %relevant elements of each file
          %intersect both indexes in order to find positios that meet
          %space-time dominum
          starttime=time(min(view_ind));
          starttime=datestr(starttime,'HHMM');
          endtime=time(max(view_ind));
          endtime=datestr(endtime,'HHMM');
          %-----------------------------------------------------------
          
    %prepare for writing interestingevents
    q=erase(filename,'.hdf');
    q=erase(q,'_NQC');
    erasethis=['_' num2str(orbit_id,'%05i')];
    q=erase(q,erasethis);
    q=erase(q,'_SC_P0.2');
    filetextname=[q,'_' starttime, '_', endtime,'_events','.txt'];
    fullfilename=fullfile(write_dir,filetextname); 
    fileID=fopen(fullfilename,'w');
    
    event_properties{1,17}=[];
    if fileID==-1
             disp('Could not open file!');
             
    else
    
    %disp(['Printing file nº' num2str(k) ', ' filetextname ' @ ' write_dir]);
    disp([num2str(k*100/ninterestingfiles) '% of .txt files printed']);
       
       %Write Header of the file
       
       fprintf(fileID, 'TAI93_time e_lat e_lon e_radiance group g_lat g_lon flash f_lat f_lon area a_lat a_lon a_observe_time x_pixel y_pixel bg_radiance\r\n');
       
          
          
          
        
       
        j=1; %nº of intersting events in each file
 for i=1:length(event.tai_time) %HERE THE PROGRAM PRINTS
            
            eventinsidetime=false;
            eventinsiderange=false;
          
    %Is this event interesting?
    %----------------------------------------------------------------------
         coordinates=event.location(:,i);
         c_distance=sqrt((coordinates(1)-centroid(1))^2+(coordinates(2)-centroid(2))^2);
         if c_distance<ang_range
           eventinsiderange=true;
         end
         
         %check time domain
         time=datetime(1993,1,1)+seconds(event.tai_time(i));
         if time<timerange(2) && time>timerange(1)
           eventinsidetime=true;
         end
    %----------------------------------------------------------------------    
           
    if eventinsidetime==true && eventinsiderange==true
  
    
        %Make new clear info with time event, its radiance, adress, group
        %adress, flash adress and area adress        

            tai_time=event.tai_time(i);
            event_lat=event.location(1,i);
            event_lon=event.location(2,i);
            radiance=event.radiance(i);
            bg_rad=event.bg_rad(i);
            x_pixel=event.x_pixel(i);
            y_pixel=event.y_pixel(i);
            
            egroup=event.parent_address(i);%his group address
            group_lat=group.location(1,egroup);
            group_lon=group.location(2,egroup);
            
            eflash=group.parent_address(egroup);%his flash address
            flash_lat=flash.location(1,eflash);
            flash_lon=flash.location(2,eflash);
            
            earea=flash.parent_address(eflash);%his area address
            area_lat=area.location(1,earea);
            area_lon=area.location(2,earea);
            area_observetime=area.observetime(earea);
        
         fprintf(fileID, '%-18.16E %7.3f %7.3f %4.1d %u %7.3f %7.3f %u %7.3f %7.3f %u %7.3f %7.3f %u %u %u %u\r\n',tai_time,event_lat,event_lon,radiance,egroup,group_lat,group_lon,eflash,flash_lat,flash_lon,earea,area_lat,area_lon,area_observetime,x_pixel,y_pixel,bg_rad); 
    
         a={ tai_time,event_lat,event_lon,radiance,egroup,group_lat,group_lon,eflash,flash_lat,flash_lon,earea,area_lat,area_lon,area_observetime,x_pixel,y_pixel,bg_rad};
         event_properties(j,:)=a;
         %the group, event flash and area address will not be ordered nor
         %listed from 1 to X, maybe from 543 to 987. That's because the
         %parent address of an event is reset in each file (1-to-end), but
         %when we recheck if an event is inside the spacew-time dominum we
         %will eliminate, for example, all events that were part of de
         %groups 1 to 453.
         
         if isempty(event_properties{j,4})
             disp('Radiance empty. code line 450 aprox');
         end
         
         j=j+1;
         
         
    end
       
   
  end
       
        fclose(fileID);
        

    end
    
    
    if n==2 
        disp('Plotting events form the current file...');
        plot_events(event_properties,centroid,ang_range);
       
        
    end
end 
    
end

function read_hdf4_and_plot(interestingfiles,centroid,ang_range,timerange)

disp('Reading data from interestingfiles and plotting...');
for k=1:length(interestingfiles) %We will go through all interestingfiles
                                 %and search for interestingevents
    
         eventinsiderange=false;
         eventinsidetime=false; 
         
         
         filename=interestingfiles(k).Filename;
         fileinfo=hdfinfo(filename);
        
         event_vdata=hdfread(fileinfo.Vgroup.Vgroup.Vgroup.Vdata(4));
         event.location=event_vdata{3};
         event.tai_time=event_vdata{1};
         
         area_vdata=hdfread(fileinfo.Vgroup.Vgroup.Vgroup.Vdata(1));
         flash_vdata=hdfread(fileinfo.Vgroup.Vgroup.Vgroup.Vdata(2));
         group_vdata=hdfread(fileinfo.Vgroup.Vgroup.Vgroup.Vdata(3));
         event_vdata=hdfread(fileinfo.Vgroup.Vgroup.Vgroup.Vdata(4));
         
        event.location=event_vdata{3};
        event.tai_time=event_vdata{1};
        event.tai_time=event_vdata{1};       
        event.observe_time=event_vdata{2};
        event.radiance=event_vdata{4};
        event.address=event_vdata{6}+1; %add +1 so it doesn't start w/ 0 adress
        event.parent_address=event_vdata{7}+1;
        event.bg_radiance=event_vdata{11};
        
        
        group.parent_address=group_vdata{7}+1; %no need to store group.address as they are ordered by number inside each file
        group.location=group_vdata{3};
        
        flash.parent_address=flash_vdata{8}+1;
        flash.location=flash_vdata{4};
        
        area.location=area_vdata{4};
        area.observetime=area_vdata{3};
    
    %prepare for writing interestingevents

       
       %Write Header of the file
       

      event_properties{1,14}=[];
        for i=1:length(event.tai_time) 
          
    %Is this event interesting?
    %----------------------------------------------------------------------
         coordinates=event.location(:,i);
         c_distance=sqrt((coordinates(1)-centroid(1))^2+(coordinates(2)-centroid(2))^2);
         if c_distance<ang_range
           eventinsiderange=true;
         end
         
         %check time domain
         time=datetime(1993,1,1)+seconds(event.tai_time(i));
         if time<timerange(2) && time>timerange(1)
           eventinsidetime=true;
         end
    %----------------------------------------------------------------------    
            if eventinsidetime==true && eventinsiderange==true
  
    
        %Make new clear info with time event, its radiance, adress, group
        %adress, flash adress and area adress        

            tai_time=event.tai_time(i);
            event_lat=event.location(1,i);
            event_lon=event.location(2,i);
            radiance=event.radiance(i);
            
            egroup=event.parent_address(i);%his group address
            group_lat=group.location(1,egroup);
            group_lon=group.location(2,egroup);
            
            eflash=group.parent_address(egroup);%his flash address
            flash_lat=flash.location(1,eflash);
            flash_lon=flash.location(2,eflash);
            
            earea=flash.parent_address(eflash);%his area address
            area_lat=area.location(1,earea);
            area_lon=area.location(2,earea);
            area_observetime=area.observetime(earea);

             a={ tai_time,event_lat,event_lon,radiance,egroup,group_lat,group_lon,eflash,flash_lat,flash_lon,earea,area_lat,area_lon,area_observetime};
             event_properties(i,:)=a;
            end
            
             eventinsiderange=false;
             eventinsidetime=false; 
        end
       
    plot_events(event_properties,centroid,ang_range);
end

 

end

function read_txt_and_plot(read_dir,centroid,ang_range,timerange)

addpath(read_dir);
folderinfoprev=dir(read_dir);
folderinfo=folderinfoprev(~ismember({folderinfoprev.name},{'.','..','.DS_Store'})); %.DS_Store is a metadata file created by iOS environment

      for k=1:size(folderinfo,1)
          filename=folderinfo(k).name;
          
          data=import_event_LIS_txtfiles(filename);
          event_properties=[];
          time=datetime(1993,1,1)+seconds(data(:,1));
          lats=data(:,2);
          lons=data(:,3);
          
          ang_distances=sqrt((lats-centroid(1)).^2+(lons-centroid(2)).^2);
          
          space_ind=find(ang_distances<ang_range);
          time_ind=intersect(find(time>timerange(1)),find(time<timerange(2)));
          
          view_ind=intersect(space_ind,time_ind); %relevant elements of each file
          
          event_properties(:,:)=data(view_ind,:);
          
          plot_events(event_properties,centroid,ang_range); 

          %{
          for i=1:size(data,1)
                eventinsiderange=false; %for each element we initialise inside checking
                eventinsidetime=false;
              
             coordinates=[data(i,2),data(i,3)];
             c_distance=sqrt((coordinates(1)-centroid(1))^2+(coordinates(2)-centroid(2))^2);
            if c_distance<ang_range
                 eventinsiderange=true;
            end
         
         %check time domain
             time=datetime(1993,1,1)+seconds(data(i,1));
             if time<timerange(2) && time>timerange(1)
                 eventinsidetime=true;
             end
         
            if eventinsiderange==true && eventinsidetime==true
                 event_properties(i,:)=data(i,:);
            end
             
         
          end
          
         
          %call this function for each file. If the .txt file contains
          %interesting events or not doesn't matter. If event_properties is
          %empty the function plot_events will simply do nothing.
          %This doesn't happen with HDF4 files because in that case we
          %first search for interesting files that they already contain
          %interesting information.
         %}
     
      end
end

function correct_hdf_urls(urls_file,relevant_orbits_file,corrected_urls_file_dir)

disp('Reading urls and relevant orbits txt file...');
relevant_orbits=import_relevant_orbits(relevant_orbits_file);
allurls=import_urls_times(urls_file);
urls=import_GHRCURLs(urls_file);

disp('Reading URLs timing...');
for i=1:length(allurls)
    
    date=num2str(allurls(i,1));
    year=str2double(date(1:4));
    month=str2double(date([5 6]));
    day=str2double(date([7 8]));
    orbit_days_from_urls(i)=datetime(year,month,day);

    
end

days_that_LIS_passes_the_zone=relevant_orbits(:,1);    
        
unique_days_LIS_passes=table2cell(unique(days_that_LIS_passes_the_zone));


disp('Checking URLs coincidence with relevant orbits...');
eliminate_indexes=[];
for i=1:size(orbit_days_from_urls,2)
    
    finding=false;
    for j=1:size(unique_days_LIS_passes,1) 
        if orbit_days_from_urls(i) == unique_days_LIS_passes{j,1}
            
            finding=true;
            
        end
    end
    
    if finding==false
        try
        
        eliminate_indexes=[eliminate_indexes; i];
        
        catch
            
            disp('error line 673');
        end
    end
    
    disp([num2str(i*100/size(orbit_days_from_urls,2)) '% checked.']);

end
 
 disp(['Found ' num2str(length(eliminate_indexes)) ' non-relevant URLs']);
 urls(eliminate_indexes,:)=[];

 fullfilename=fullfile(corrected_urls_file_dir,'interesting_URLs.txt');
 writetable(urls,fullfilename);
 disp(['Urls .txt file written at:  ' corrected_urls_file_dir]);
 

end
%plotters
function plot_coastline(read_dir,centroid,ang_range)

load coastlines;
addpath(read_dir);
disp('Plotting coastline...');


latlim=[ centroid(1)-ang_range centroid(1)+ang_range ];
lonlim=[centroid(2)-ang_range centroid(2)+ang_range];

tf=ingeoquad(coastlat,coastlon,latlim,lonlim);
lats=coastlat(tf); 
lons=coastlon(tf);

hold on;
plot(lons,lats,'k-');


end
function plot_events(event_properties,centroid,ang_range)
        
hold on;

    if ~isempty(event_properties)

        if iscell(event_properties)
      
        size=cell2mat(event_properties(:,4));
        color=cell2mat(event_properties(:,5));
        lats=cell2mat(event_properties(:,2));
        lons=cell2mat(event_properties(:,3));
       
    
     
        else
  
        size=event_properties(:,4);
        color=event_properties(:,5);
        lats=event_properties(:,2);
        lons=(event_properties(:,3));
        end
        
        plot_interesting_area(centroid,ang_range);
        scatter(lons,lats,size,color,'x');
        xlabel('Longitude [º]');
        ylabel('Latitude [º]');
        grid on;
        title('Events detected. Size -> Radiance. Color -> group');
        axis equal;


    else
        disp('Event_properties vector is empty!');
    end   

    
end       
function plot_interesting_area(centroid,r)
x=centroid(2);
y=centroid(1);

hold on;
grid on;

th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
plot(xunit, yunit,'k--','LineWidth',0.01);
%plot(x,y,'kx','LineWIdth',0.01);
axis equal;

end

%txt reading. Made by MATLAB
function data = import_event_LIS_txtfiles(filename)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   DATA = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   DATA = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   data = importfile('ISS_LIS_SC_P0.2_20180429_NQC_07514_event.txt', 2, 47);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/07/25 13:04:15

% Initialize variables.
delimiter = {' ',' '};
if nargin<=2
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
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

% Close the text file.
fclose(fileID);

%Post processing for unimportable data.
%No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

% Create output variable
data = [dataArray{1:end-1}];
end

function GHRCURLs = import_urls_times(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   GHRCURLS = IMPORTFILE(FILENAME) Reads data from text file FILENAME for
%   the default selection.
%
%   GHRCURLS = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   GHRCURLs = importfile('GHRC_URLs.txt', 1, 433);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/07/23 16:23:22

% Initialize variables.
delimiter = {'_','.'};
if nargin<=2
    startRow = 1;
    endRow = inf;
end

% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%*s%*s%*s%*s%*s%*s%*s%s%s%s%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%Close the text file.
fclose(fileID);

% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


% Exclude columns with non-numeric cells
I = ~all(cellfun(@(x) (isnumeric(x) || islogical(x)) && ~isnan(x),raw),1); % Find columns with non-numeric cells
raw(:,I) = [];

 %Initialize column outputs.
columnIndices = cumsum(~I);

% Create output variable
GHRCURLs = cell2mat(raw);

end

function relevantorbits = import_relevant_orbits(filename)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   RELEVANTORBITS = IMPORTFILE(FILENAME) Reads data from text file
%   FILENAME for the default selection.
%
%   RELEVANTORBITS = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   relevantorbits = importfile('relevant_orbits.txt', 1, 1021);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/07/23 15:13:23

% Initialize variables.
delimiter = ' ';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

% Format for each line of text:
%   column1: datetimes (%{yyyy-MM-dd}D)
%	column2: datetimes (%{HH:mm:ss}D)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%{yyyy-MM-dd}D%{HH:mm:ss}D%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
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

% Create output variable
relevantorbits = table(dataArray{1:end-1}, 'VariableNames', {'year_day','hour'});

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% relevantorbits.year_day=datenum(relevantorbits.year_day);
% relevantorbits.hour=datenum(relevantorbits.hour);
end 
function GHRCURLs = import_GHRCURLs(filename)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   GHRCURLS = IMPORTFILE(FILENAME) Reads data from text file FILENAME for
%   the default selection.
%
%   GHRCURLS = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   GHRCURLs = importfile('GHRC_URLs.txt', 1, 433);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/07/23 15:39:58

% Initialize variables.
delimiter = {''};
if nargin<=2
    startRow = 1;
    endRow = inf;
end

% Format for each line of text:
%   column1: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%[^\n\r]';

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
    dataArray{1} = [dataArray{1};dataArrayBlock{1}];
end

% Close the text file.
fclose(fileID);

%Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

% Create output variable
GHRCURLs = table(dataArray{1:end-1}, 'VariableNames', {'httpsghrcnsstcnasagovpublisissdatasciencenqchdf20170401ISS_LIS_'});

end

function  webfilenames2urls(write_dir,read_filename)

website='https://ghrc.nsstc.nasa.gov/pub/lis/iss/data/science/nqc/hdf/';

files= import_web_filenames(read_filename);
%urls generation
disp('Creating url vector...');

write_filename=fullfile(write_dir,'interesting_URLs.txt');
fileID=fopen(write_filename,'w');

disp(['Printing @ ' write_dir]);
    for i=1:size(files,1)
        
        name=char(files{i,1});
        year=name([17:20]);
        month=name([21 22]);
        day=name([23 24]);         
        url(i,:)=[website year '/' month day '/' name];
        fprintf(fileID,[url(i,:) '\r\n']);
        
    end
    fclose(fileID);
end
function files = import_web_filenames(filename)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   FILES = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   FILES = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   files = importfile('filenames_website.txt', 2, 30);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/07/24 10:17:13

% Initialize variables.
delimiter = {'\t','?','[',']',' '};
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%Format for each line of text:
%   column1: text (%s)
%	column2: datetimes (%{MMM}D)
%   column3: datetimes (%{dd}D)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%{MMM}D%{dd}D%*s%*s%*s%*s%*s%*s%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
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

% Create output variable
files = table(dataArray{1:end-1}, 'VariableNames', {'name','month','month_day'});

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% files.month=datenum(files.month);
% files.month_day=datenum(files.month_day);

end
