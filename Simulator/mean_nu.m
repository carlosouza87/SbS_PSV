function [nu_m] = mean_nu(nu,wts)

lw = length(wts);
lu = length(nu);
nu_ext_beg = nu(1:(lw-1)/2);
nu_ext_end = nu(lu-(lw-1)/2+1:lu);
nu_ext = [nu_ext_beg nu nu_ext_end];

% mav_nu = conv(nu,wts,'valid');
mav_nu = conv(nu_ext,wts,'valid');

% lm = length(mav_nu);

% [p_beg] = polyfit((lw-1)/2+1:lw-1,mav_nu(1:(lw-1)/2),1);
% [p_end] = polyfit(lm-(lw-1)/2+1:lm,mav_nu(lm-(lw-1)/2+1:lm),1);
% nu_m_beg = polyval(p_beg,1:(lw-1)/2);
% nu_m_end = polyval(p_end,lm+1:lm+(lw-1)/2);  

% nu_m_beg = mean(nu(1:(lw-1)/2))*ones(1,(lw-1)/2);
% nu_m_end = mean(nu(lm+1:lm+(lw-1)/2))*ones(1,(lw-1)/2);
% nu_m = [nu_m_beg mav_nu nu_m_end];

nu_m = mav_nu;