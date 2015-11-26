% 1st order
if isimtype == 1
    w1st_file = [caseid '.3'];
    scl_factor = rho*g; % Scaling factor    
elseif isimtype == 2
    w1st_file = [caseid '.4'];
    scl_factor = 1; % Scaling factor 
else
    error('Invalid value for isimtype!')
end

%% FALTA DAQUI PRA BAIXO!!!

[periods,incid,dof,amp,pha,Re,Im] = textread(w1st_file);

periods = unique(periods);
incid = unique(incid);

[freqs,idx] = sort(2*pi./periods);
nfreqs = length(freqs);
nincid = length(incid);
ndof = length(amp)/(nfreqs*nincid);
rep_per = 6*nincid;
amp = amp*rho*g;

cont = 0;
for k1 = 1:nfreqs
    for k2 = 1:nincid
        for k3 = 1:6
            FTF_amp1(nfreqs-k1+1,k2,k3) = amp(cont+k3,1);
            FTF_pha1(nfreqs-k1+1,k2,k3) = pha(cont+k3,1);
            FTF_amp2(nfreqs-k1+1,k2,k3) = amp(cont+k3+6,1);
            FTF_pha2(nfreqs-k1+1,k2,k3) = pha(cont+k3+6,1);
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

cont = 0;
for k1 = 1:nfreqs
    for k2 = 1:nincid
        for k3 = 1:6
            drift_amp1(nfreqs-k1+1,k2,k3) = amp(cont+k3,1);
            drift_pha1(nfreqs-k1+1,k2,k3) = pha(cont+k3,1);
            drift_amp2(nfreqs-k1+1,k2,k3) = amp(cont+k3+9,1);
            drift_pha2(nfreqs-k1+1,k2,k3) = pha(cont+k3+9,1);
        end
        cont = cont + ndof;
    end
end

save waveforces FTF_amp1 FTF_pha1 FTF_amp2 FTF_pha2 drift_amp1 drift_pha1 drift_amp2 drift_pha2 incid freqs
