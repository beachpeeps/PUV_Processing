
function sigIN = get_directionality(sig,U,directionality)
        switch directionality
            case 'shoreward'
                %disp('Removing reflection')
                %% Remove reflection
                H = mean(sig)';
                u = U;%(:,1:size(sig,2))';
                fmin = 0.004; fmax=0.04;
                for i=1:size(sig,2)
                    [sigIN(:,i),sigOUT(:,i)] = RemoveReflection( sig(:,i), u(:,i), H(i),  0.5,'method','sher','flim', [fmin,fmax],'g',9.81,'trend',true,'depav',true);
                end
            case 'all'
                sigIN = sig;
                disp('Not removing reflection')
            otherwise
                warning('Directionality not correctly specified.')
        end
    end
