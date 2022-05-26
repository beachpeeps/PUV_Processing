function res = Swash_Lay2Arr(res,vartype)

if vartype == 0
    %
    % Nprs, non-hydrostatic pressure
    namebase = 'Nprs_k'; 
    nameend  = '';
    outname  = 'Pnh';
    im       = 0;
    %
elseif vartype == 1
    %
    % Nprs, non-hydrostatic pressure
    namebase = 'Pres_k'; 
    nameend  = '';
    outname  = 'P';
    im       = 0;
    %  
elseif vartype == 2
    %
    % Nprs, non-hydrostatic pressure
    namebase = 'vel_k'; 
    nameend  = '_x';
    outname  = 'U';
    im       = 0;
    %    
elseif vartype == 3
    %
    % Nprs, non-hydrostatic pressure
    namebase = 'vel_k'; 
    nameend  = '_y';
    outname  = 'V';    
    im       = 0;
    % 
elseif vartype == 4
    %
    % Nprs, non-hydrostatic pressure
    namebase = 'w'; 
    nameend  = '';
    outname  = 'W';    
    im       = 1;
    % 
elseif vartype == 5
    %
    % Nprs, non-hydrostatic pressure
    namebase = 'Mvel_k'; 
    nameend  = '_x';
    outname  = 'Um';    
    im       = 0;
    %        
    % 
elseif vartype == 6
    %
    % Nprs, non-hydrostatic pressure
    namebase = 'Tke'; 
    nameend  = '';
    outname  = 'Tke';    
    im       = 0;
elseif vartype == 7
    %
    % Nprs, non-hydrostatic pressure
    namebase = 'MEps'; 
    nameend  = '';
    outname  = 'eps';    
    im       = 0;    
    % 
elseif vartype == 8
    %
    % Nprs, non-hydrostatic pressure
    namebase = 'zk'; 
    nameend  = '';
    outname  = 'zk';    
    im       = 0;  
elseif vartype == 9
    %
    % Nprs, non-hydrostatic pressure
    namebase = 'omega'; 
    nameend  = '';
    outname  = 'omega';
    im       = 1;      
end



fnd = true;
kmax=1;
fmt = '%02.0f';
names = fieldnames(res);

looping = true;
while (true)
    while ( fnd )
        %
        name = [namebase,num2str(kmax-im,fmt),nameend];

        foun = find(strcmpi(name,names),1);
        if ~isempty(foun)
             ny = size( res.(name),1 );
             nx = size( res.(name),2 );
             nt = size( res.(name),3 );        
             fnd = true;
             kmax = kmax+1;
        else       
             kmax = kmax-1;
             break;
        end
        %
    end

    if (looping)
        if kmax<1
            fmt = '%1.0f';
            kmax=1;
            looping = false;
        else
            %looping = false;
            break;
        end
    else
       break 
    end
end



res.(outname) = zeros( ny , nx , nt , kmax );
for ii = 1:kmax
    name = [namebase,num2str(ii-im,fmt),nameend];
    dat = res.(name);
    res.(outname)(:,:,:,ii) = dat;
end

res.(outname) = squeeze(res.(outname));

for ii = 1:kmax
    name = [namebase,num2str(ii-im,fmt),nameend];
    res = rmfield(res,name);
end
 