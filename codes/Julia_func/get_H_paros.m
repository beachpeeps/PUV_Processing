function [X, Hinc, Hig, nameSensors] = get_H_paros(datahour)
% returns X locations, with [Hinc/Hig, upperbound, lowerbound]

currdir = pwd;
load ~/Agate/mat/sensorlocs.mat X indP namesP
XX = X;
sensors = [2:6 8:11 13 15:18];
X = XX(indP); X= X(sensors);
X = [-X' -800 -1000 -1200 -1400];
namesP = namesP(sensors);
%%
Hz = 2;

for i=1:length(namesP)
    s = char(namesP(i));
    nsensor(i) = str2double(s(2:3));
end

for i=1:length(nsensor)
    %
    cd ~/Agate/mfiles
    %
    
    try
        [fm,Spp,H,~,dof] = parosspec_agate(nsensor(i),datahour,Hz);
        % find freqs in inc and ig
        nIG = find(fm>=0.004 & fm<=0.04);
        nINC = find(fm>=0.04 & fm<=0.25);
        df = fm(2)-fm(1);
        %
       [Hig_P(i), lbIG(i), ubIG(i), edofIG(i) ]...
            = getSWHebounds( Spp(nIG), dof, 0.95, df );
        [Hinc_P(i), lb(i), ub(i), edof(i) ]...
            = getSWHebounds( Spp(nINC), dof, 0.95, df );
    
    end
    
    cd(currdir)
    clear H fm Spp
end
%%

PUVnum = [5 6 7 8];


for ii=1:4
    n = i+ii;
    try
        [H,fm,Spp] = get_PUVspectrahour_agate(PUVnum(ii),datahour);
        nIG = find(fm>0.004 & fm<0.04);
        nINC = find(fm>0.04 & fm<0.25);
        df = fm(2)-fm(1);
       [Hig_P(n), lbIG(n), ubIG(n), edofIG(n) ]...
            = getSWHebounds( Spp(nIG), dof, 0.95, df );
        [Hinc_P(n), lb(n), ub(n), edof(n) ]...
            = getSWHebounds( Spp(nINC), dof, 0.95, df );
        
    catch
    end
end
%%

Hinc = [Hinc_P; ub; lb];
Hig = [Hig_P; ubIG; lbIG];
%%
namesPUV = [];
for i=1:4
    namesPUV{i} = sprintf('PUV%02.0f',PUVnum(i));
end

nameSensors = [namesP namesPUV];
