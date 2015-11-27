% 1st order
if isimtype == 1
    w1st_file = [caseid '.3'];
    scl_factor = rho*g; % Scaling factor   
    m = 2; % Exponent for ULEN during scaling  
elseif isimtype == 2
    w1st_file = [caseid '.4'];
    scl_factor = 1; % Scaling factor
    m = 0; % Exponent for ULEN during scaling  
else
    error('Invalid value for isimtype!')
end

[periods,incid,dof,amp,pha,Re,Im] = textread(w1st_file);

periods = unique(periods);
incid = unique(incid);

[freqs,idx] = sort(2*pi./periods);
nfreqs = length(freqs);
nincid = length(incid);
ndof = length(amp)/(nfreqs*nincid);
rep_per = 6*nincid;
amp = amp*rho*g;

w1st_amp1 = zeros(nfreqs,nincid,6);
w1st_amp2 = zeros(nfreqs,nincid,6);
w1st_pha1 = zeros(nfreqs,nincid,6);
w1st_pha2 = zeros(nfreqs,nincid,6);

cont = 0;
for k1 = 1:nfreqs
    for k2 = 1:nincid
        for k3 = 1:6
            if k1 <= 3
                w1st_amp1(nfreqs-k1+1,k2,k3) = amp(cont+k3,1)*scl_factor*ULEN^m;
                w1st_amp2(nfreqs-k1+1,k2,k3) = amp(cont+k3+6,1)*scl_factor*ULEN^m;
            else
                w1st_pha1(nfreqs-k1+1,k2,k3) = amp(cont+k3,1)*scl_factor*ULEN^(m+1);
                w1st_pha2(nfreqs-k1+1,k2,k3) = amp(cont+k3+6,1)*scl_factor*ULEN^(m+1);
            end            
            w1st_pha1(nfreqs-k1+1,k2,k3) = pha(cont+k3,1);            
            w1st_pha2(nfreqs-k1+1,k2,k3) = pha(cont+k3+6,1);
        end
        cont = cont + ndof;
    end
end

clear periods dof amp pha Re Im
% 2nd order
driftfile  = ['conjunto.9'];

[periods,incid1,incid2,dof,amp,pha,Re,Im] = textread(driftfile);

unique_periods = unique(periods);
unique_incid = unique(incid1);

[freqs,idx] = sort(2*pi./unique_periods);
nfreqs = length(freqs);
nincid = length(unique_incid);
ndof = length(amp)/(nfreqs*nincid);
amp = amp*rho*g;

w2nd_amp1 = zeros(nfreqs,nincid,6);
w2nd_amp2 = zeros(nfreqs,nincid,6);
w2nd_pha1 = zeros(nfreqs,nincid,6);
w2nd_pha2 = zeros(nfreqs,nincid,6);

cont = 0;
for k1 = 1:nfreqs
    for k2 = 1:nincid
        for k3 = 1:6
            w2nd_amp1(nfreqs-k1+1,k2,k3) = amp(cont+k3,1);
            w2nd_amp2(nfreqs-k1+1,k2,k3) = pha(cont+k3,1);
            w2nd_pha1(nfreqs-k1+1,k2,k3) = amp(cont+k3+9,1);
            w2nd_pha2(nfreqs-k1+1,k2,k3) = pha(cont+k3+9,1);
        end
        cont = cont + ndof;
    end
end

save('waveforces','w1st_amp1','w1st_amp2','w1st_pha1','w1st_pha2','w2nd_amp1','w2nd_amp2','w2nd_pha1','w2nd_pha2','incid','freqs')
