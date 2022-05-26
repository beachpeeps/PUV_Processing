
function dFuncdx = derivative( func , dx )
    %
    nt = size(func,1);

    if ndims(func) == 3
        ns = size(func,2);
        nx = size(func,3);
    else
        ns = 1;
        nx = size(func,2);
    end
            
    if ns == 1
        dFuncdx = zeros( nt,nx );
        %
        for jj=2:size(func,2)-1
            %
        
        
            dFuncdx(:,jj) = ( func( :,jj+1 ) - func(:,jj-1) ) ./ (dx(jj-1) +dx(jj));
            %
        end
        %
        dFuncdx(:,1) = nan;
        dFuncdx(:,end) = nan;
    else
        dFuncdx = zeros( nt,ns,nx );
        for jj=2:size(func,3)-1
            %        
            dFuncdx(:,:,jj) = ( func( :,:,jj+1 ) - func(:,:,jj-1) ) ./ (dx(jj-1) +dx(jj));
            %
        end
        %
        dFuncdx(:,:,1) = nan;
        dFuncdx(:,:,end) = nan;
    end
end