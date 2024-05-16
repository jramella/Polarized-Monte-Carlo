function U2=updateU(U,phi,theta)

ux = U(1);
uy = U(2);
uz = U(3);

costheta = cos(theta);
sintheta = sqrt(1.0 - costheta*costheta);
cosphi   = cos(phi);

if phi < pi
    sinphi = sqrt(1.0 - cosphi*cosphi);
else
    sinphi = -sqrt(1.0 - cosphi*cosphi);
end



% New directional cosines. */
if (1 - abs(uz) <= 1.0e-12)
    uxx = sintheta * cosphi;
    uyy = sintheta * sinphi;
    uzz = costheta * sign(uz);
    
else
    temp = sqrt(1.0 - uz * uz);
    uxx = sintheta * (ux * uz * cosphi - uy * sinphi) / temp + ux * costheta;
    uyy = sintheta * (uy * uz * cosphi + ux * sinphi) / temp + uy * costheta;
    uzz = -sintheta * cosphi * temp + uz * costheta;
end
% Update directional cosines */
U2(1) = uxx;
U2(2) = uyy;
U2(3) = uzz;

