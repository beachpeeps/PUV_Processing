function [fc, fsp, fkurt, fp] = get_fstats(frequency,Spectrum,df)
%%
f = frequency;
y = Spectrum;
s = find(f>=0.04 & f<=0.25); % do SS freqs only
s = find(f>=0.04 & f<=0.4); % do all freqs only

%% need to write some test cases to test this function!!
%  sum ( f(iswell).* MOP.spec1D(:,iswell)' ) ./ sum ( MOP.spec1D(:,iswell),2 )' ;
% df = double(MOP.fbounds(2,:)-MOP.fbounds(1,:))';
% y(1,:) = normpdf(f,0.1,0.02);
% y(2,:) = 100*normpdf(x,1,2);
% y(3,:) =  MOP.spec1D(95019,:);


fs = f(s);
dfs = df(s);
[m,n] = size(y);
% if m>n
    ys = y(:,s);
% elseif m<n
%     ys = y(:,s);
% end

%fpeak
for i=1:m
    fpind(i) = find(ys(i,:) == max(ys(i,:)));
end

fp = fs(fpind);

% f centroid (aka mean)
fc = sum ( fs.* ys' .*dfs ) ./ sum ( ys' .* dfs) ;

% f spread (aka standard deviation)
fsp = sqrt( sum( (fs-fc).^2.* ys' .* dfs ) ./ sum (ys' .* dfs) );

% kurtosis (4th central moment normalized by fsp.^4)
fkurt = sum( (fs-fc).^4.* ys' .* dfs) ./ sum (ys' .* dfs) ./ fsp.^4 ;
