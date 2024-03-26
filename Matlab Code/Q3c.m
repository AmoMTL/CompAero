cdf80 = 0;
cdf84 = 0;
cdf86 = 0;
cdp80 = 0;
cdp84 = 0;
cdp86 = 0;
for i=1:49
    cdf80 = cdf80 + CF80(i)*dx_2;
    cdf84 = cdf84 + CF84(i)*dx_2;
    cdf86 = cdf86 + CF86(i)*dx_2;
    
    cdp80 = cdp80 - cp80(i)*(y1(i+1) - y1(i));
    cdp84 = cdp84 - cp84(i)*(y1(i+1) - y1(i));
    cdp86 = cdp86 - cp86(i)*(y1(i+1) - y1(i));
end
cd80 = cdf80 + cdp80
cd84 = cdf84 + cdp84
cd86 = cdf86 + cdp86