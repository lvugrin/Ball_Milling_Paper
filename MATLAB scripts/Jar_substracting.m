rs = Monitoring_spectra_raw(:,2:end-1);
ls = size(rs,1);
ns = size(rs,2); 
ds = Spectra.DarkSpectrum;

ss = rs;
for j=1:ns,
        ss(:,j) = rs(:,j)-ds(:,2);
end

% Removing Cosmic Rays
% ss = Cosmic(ss,484,1605);
%ss = Cosmic(ss,484,1342);

clear j

figure,plot(Monitoring_spectra_raw)
%% Removal of Cosmic Rays around peak
p = 1381; % point No. for peak maximum for the jar signal (1428)
p1 = 130; % Number of points left and right of the peak maximum
jar = Cosmic(Spectra.JarSpectrum,484,719); 
jarbgr = jar(p-p1:p+p1, 2); 
jarforscale = jarbgr;
s = size(jarbgr, 1);

for z=1:ns,  % Number of spectra it runs through
    sfs = ss(:, z);
	sfsbgr = sfs(p-p1:p+p1,1);
    sfsforscale = sfsbgr;
    for k = 1:5000,
        for x = 2:s-1,
          if sfsbgr(x,1)>((sfsbgr(x-1,1)+sfsbgr(x+1,1))/2),
              sfsbgr(x,1)=(sfsbgr(x-1,1)+sfsbgr(x+1,1))/2;
          end
          if jarbgr(x,1)>((jarbgr(x-1,1)+jarbgr(x+1,1))/2),
              jarbgr(x,1)=(jarbgr(x-1,1)+jarbgr(x+1,1))/2;
          end
       end
    end
    % position of the peak maximum is now at p1
    p2 = 6; %Number of points left and right to calculate mean scale
    b1 = jarforscale(p1-p2:p1+p2,1);
    b2 = jarbgr(p1-p2:p1+p2,1);
    b = b1-b2;								% background subtraction of jar spectrum
    
	a1 = sfsforscale(p1-p2:p1+p2,1);
    a2 = sfsbgr(p1-p2:p1+p2,1);
    a = a1-a2;								% background subtraction of exp spectrum
    
    scale = mean(a./b);				
    jarscaled = jar(:,2).*scale; 
	    
	ps(:,z) = ss(:,z)-jarscaled(:,1);		% subtracting the jar spectrum from exp spectrum
    
    New_spectra_without_jar = ps;           % there is no wavenumbers in the first column. 1st column is the first spectrum.
                                            % Work with this variable as with Monitoring_spectra_without_jar
    
    clear  a a1 a2 ans b b1 b2 jarscaled scale
    
end

clear z;

clear  jarbgr jarforscale  k  sfs sfsbgr sfsforscale x z ds  s ls ns p p1 p2 ps rs ss

figure,plot(New_spectra_without_jar);
