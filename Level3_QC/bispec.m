function [data] = bispec(x, alpha, f_bound)
    fs = 2;
    [data] = fun_compute_bispectrum_H1982_nowindow(x, fs, 50); % frequency merge 10

    alpha = interp1(f_bound, alpha, data.f);

    fcutoff = 0.25;
    
    [~,fid_min]=min(abs(data.f - 0.004)); if data.f(fid_min) == 0; fid_min = fid_min + 1;end
    [~,fid_max]=min(abs(data.f - 0.25));
    data.b_df = NaN(1, length(data.f));
    data.E_forced = NaN(1, length(data.f));
    clear f2 f2_id

    % need to normalize
    data.B = data.B.*(2./(data.df.^2));
    data.P = data.P.*(2./data.df);
    for ii = fid_min:fid_max % loop through df
        clear f2_id f2
        num = 2*nansum(data.B(ii, ii:fid_max).*data.df);
        f2=data.f(ii)+data.f(ii:fid_max);
        for ff = 1:length(f2)
            f2_id(ff) = find(round(f2(ff),5)== round(data.f,5));
        end
        denom = sqrt(2*nansum(data.P(ii:fid_max)'  .* data.P(f2_id)'.* data.P(ii) .*data.df)); 
        data.b_df(ii) = num./denom;
        data.E_forced(ii) = alpha(ii) * abs(data.b_df(ii)).^2 * data.P(ii);
    end
end 
