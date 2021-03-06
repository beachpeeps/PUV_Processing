function [res,vars] = swash_loadMatFile(file)
%--------------------------------------------------------------------------
% 
% swash_loadMatFile(file)
%
% file :  Filename of matlab file as generated by swan / SWASH.
%
% -- DESCRIPTION --
%
% Function to load time dependent block data generated by swan. The
% function collects all the matrices of a single variable and puts them
% into a single 3d array. So if HSig is given as output in file then the
% mat file contains matrices
% 
% HSig_date1_time1
% HSig_date2_time2
% HSig_date3_time3
% ...          ...
% HSig_dateNT_timeNT
%
% these are collected and stored as
%
% res.HSig(NY,NX,NT)
%
% furthermore a time vector is generated which stores the time past is
% seconds measured from the date/time of the first output value
%
% res.time
%
% variables that have no time stamp (e.g. Xp, Yp) are treated using NT = 1
%
% DATE              AUTHOR              DESCRIPTION
% 30-3-2010         P.B.Smit            New function
%
%--------------------------------------------------------------------------



% Check for file existence
  if ~( exist(file,'file') == 2 )
      error('swan:io:load_SwanMatFile','File %s does not exist',file);
  end  
  
% load variables from file
  DATA = load(file);
  
% Get variable names in file
  varnames = fieldnames(DATA); %who('-file',file);

  isprocessed = false;
  for ii=1:length(varnames)
     %
     if strcmpi( varnames{ii}, 'isprocessed')
        isprocessed = true;
        break
     end
     %
  end

  if isprocessed
      %
      res  = load( file );
      vars = varnames;
      return
      %
  end
  


  
  kv     = 0;         %Number of variables stored
  vars   = cell(0,0); %Cell array containing variable names/times/info etc.
  
% Loop over all varnames
  h = waitbar(0,'Processing matfile...');
  nn = length(varnames);
  iupd=0;
  for ivars = 1:nn
   %  ivars
     iupd = iupd+1;
     if (iupd == floor(nn/100))
         %
         iupd = 0;
         waitbar( ivars/length(varnames) , h);
         %
     end
     
     nstr = length( varnames{ivars} );
    % keyboard
     if nstr > 11 && strcmpi( varnames{ivars}(nstr -10) , '_' )
        %
        %datetime = varnames{ivars}(nstr -10:nstr );
          name     = varnames{ivars}(1:nstr -11);
        
         %Remove leading _ in _date_time
         %datetime = datetime(2:length(datetime));
         %Split datetime in [date , _time]
         
         %[datestr,timestr]  = strtok(datetime,'_');   

         
         %datestr
         datestr = varnames{ivars}(nstr -9:nstr-4 );
         %timestr
         timestr  =varnames{ivars}(nstr -2:nstr );
         %pause()
         %Remove leading _ in _time
         %timestr = timestr(2:length(timestr));
         %date = [Y,M,D] datestr = YYYYMMDD         
         date = [str2double(datestr(1:2)),str2double(datestr(3:4)),str2double(datestr(5:6))];
         %time = [Hour,Min,Sec] timestr = HHMMSS
         time = str2double(timestr(1:3));
         %Convert time to seconds
         sec  = date(1)*3600 + date(2) * 60 + date(3)  + time/1000;
         %Convert day to number of days since 1 jan 0000  00:00:00
        % sec=0;
         daten = sec;
         timedep = true;        
        %
     else
        % 
        datetime = [];
        name     =  varnames{ivars};
        
        %If datetime is empty, not time dependent
        date = [0,0,0];
        time = 0;
        sec  = 0;
        daten = 0;
        timedep = false;        
        %
     end
     
     jv = 1;
     while (true)
        if jv>kv
            %Name is unknown, add new name
            vars{jv}.name = name;
            vars{jv}.nt   = 1;
            vars{jv}.dim  = size(DATA.(varnames{ivars}));
            
            vars{jv}.fullname      = char(zeros(nn,length(varnames{ivars})));
            vars{jv}.fullname(1,:) = varnames{ivars};
            
            %vars{jv}.date = zeros(nn,3);            
            
            
            vars{jv}.date = zeros(nn,3);            
            vars{jv}.date(1,:) = date;
            vars{jv}.t   = zeros(nn,1);
            vars{jv}.t(1)    = sec;
            vars{jv}.datenum = zeros(nn,1);
            vars{jv}.datenum(1) = daten;
            kv = kv+1;
            break;
        end
        
        if strcmpi(vars{jv}.name,name)
            if ~timedep
                break
            end
            
            %Variable is already known, add time strings etc.
            nt            = vars{jv}.nt + 1;
            vars{jv}.nt   = nt;
            vars{jv}.date(nt,:) = date;
            vars{jv}.t(nt)    = sec;
            vars{jv}.fullname(nt,:) = varnames{ivars};            
            vars{jv}.datenum(nt) = daten;
            break
        end
        jv = jv + 1;
    end
    
  end  
  close(h)
  
  % Loop over all varnames
  h = waitbar(0,'Saving variables...');
  res.time = [];
  for iv = 1:length(vars)
      %
      waitbar( iv/length(vars) , h);
      varname = vars{iv}.name;
      %
      eval([varname,' = zeros([vars{iv}.dim,vars{iv}.nt]);']);
      datetime   = vars{iv}.datenum(1:vars{iv}.nt);
      [t,I]      = sort(datetime);
      %

      if isempty(time) || length(time) < length(vars{iv}.t)
         
         time = vars{iv}.t(1:vars{iv}.nt);
         time = time(I);%+(t-t(1));
      end

      for in=1:vars{iv}.nt
          jn = I(in);
          eval([varname,'(:,:,jn) = DATA.(vars{iv}.fullname(in,:));']);
      end
      
      if iv==1
          isProcessed = true;
          save( file, 'isProcessed','-v7.3' );
      end
      save( file, varname ,'-append' );
      %
  end
  % Write time only once, should be same for all variables and is slow if
  % overwritten each time.........
  save( file, 'time' ,'-append' );
  
  clear('DATA');
  
  res = load(file);
  %Overwrite file for future fast loading
  save( file , '-struct', 'res','-v7.3' );

  close(h);
end


