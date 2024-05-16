%function [linreta]=polardecomposition(muel)
%function [rvec,dvec,R,D,depol,rotation,linreta]=polardecomposition(muel)
%TEST
% Mdelta is Depolarization
% MR is Retardation
% Md is deattenuation -> this is what Jones matrix cannot do.

function [depol,linreta,orientation,diattenuation]=decompose_images(muel);
%format double
%muel=MuellerMatrix;
I=[1 0 0;
    0 1 0;
    0 0 1];

pvec=[muel(2,1),muel(3,1),muel(4,1)]*(1/muel(1,1));
dvec=[muel(1,2),muel(1,3),muel(1,4)]*(1/muel(1,1));

D=((muel(1,2)^2+muel(1,3)^2+muel(1,4)^2)^0.5)*(1/muel(1,1));

m=(1/muel(1,1))*[muel(2,2),muel(2,3),muel(2,4);
    muel(3,2),muel(3,3),muel(3,4);
    muel(4,2),muel(4,3),muel(4,4)];

D1=(1-D^2)^0.5;

if D==0
    muel0=muel/muel(1,1);
    diattenuation=NaN;
else
    mD=D1*I+(1-D1)*dvec'*dvec/D^2;

    %disp('MD')

    MD=muel(1,1)*[1,dvec;
        dvec',mD];
    
    diattenuation = ((MD(1,2)^2+MD(1,3)^2+MD(1,4)^2)^0.5)*(1/MD(1,1));
    muel0=muel*inv(MD);
end

% disp('diattenuation')
% diattenuation

m1=[muel0(2,2) muel0(2,3) muel0(2,4);
    muel0(3,2) muel0(3,3) muel0(3,4);
    muel0(4,2) muel0(4,3) muel0(4,4)];

if (isnan(mean(m1(:)))==1)
    m1(:)=0;
else
    
    end
I0=eig(m1*m1');

m0=inv(m1*m1'+((I0(1)*I0(2))^0.5+(I0(2)*I0(3))^0.5+(I0(3)*I0(1))^0.5)*I);
m00=(I0(1)^0.5+I0(2)^0.5+I0(3)^0.5)*m1*m1'+I*(I0(1)*I0(2)*I0(3))^0.5;

if det(m1)>=0
    mdelta=m0*m00;
else   
    mdelta=-m0*m00;
end
if (isnan(mean(mdelta(:)))==1)
    mdelta(:)=0;
[v,mdeltaf] = eig(mdelta);
depol=1-(abs(mdelta(1,1))+abs(mdelta(2,2))+abs(mdelta(3,3)))/3;
depol1 =1-(abs(mdeltaf(1,1))+abs(mdeltaf(2,2))+abs(mdeltaf(3,3)))/3;
else
    
    [v,mdeltaf] = eig(mdelta);
depol=1-(abs(mdelta(1,1))+abs(mdelta(2,2))+abs(mdelta(3,3)))/3;
depol1 =1-(abs(mdeltaf(1,1))+abs(mdeltaf(2,2))+abs(mdeltaf(3,3)))/3;

end

% disp('depol')
% depol
% depol1

nul=(pvec'-m*dvec')/D1^2;

%disp('Mdelta')

Mdelta=[1 0 0 0;
    nul mdelta];

Mdeltaf =[1 0 0 0;
    nul mdeltaf];

Mdinv=inv(Mdelta);
%disp('MR')
MR=Mdinv*muel0;

trmR=(MR(2,2)+MR(3,3)+MR(4,4))/2;
argu=trmR-1/2;
if abs(argu)>1
    if argu>0
        R=acos(1);
    else
        R=acos(-1);
    end
else
    R=acos(argu);
end

cssq10=(MR(2,2)+MR(3,3))^2+(MR(3,2)-MR(2,3))^2;
tanrot=(MR(3,2)-MR(2,3))/((MR(2,2))+(MR(3,3)));
de=cssq10^0.5-1;

if de>0.999999999999
    de=1;
end

if de<-0.99999999999
    de=-1;
end

linreta=acos(de); %this is linear retardance
rotation=atan(tanrot);

if tanrot<0.000000001
    rotation=rotation+pi;
end

rotation=rotation/2;

if (MR(3,2)-MR(2,3))<0.0
    if (MR(2,2)+MR(3,3))<0.0
        rotation=rotation+pi/2;
    end
end

if (MR(3,2)-MR(2,3))<0.0
    if (MR(2,2)+MR(3,3))>0.0
        rotation=rotation+pi/2;
    end
end

if abs(MR(3,2)-MR(2,3))<=0.000000001 && abs(MR(2,2)+MR(3,3))>0.0000000001
    rotation=0;
end

% disp('rotation')
% rotation %this is optical rotation

if abs(sin(R))<=0.000000001
    a3=((1+cos(linreta))/2)^0.5;
    a1=(MR(3,4)+MR(4,3))/(4*a3);
    a2=(MR(4,2)+MR(2,4))/(4*a3);
else
    D2=1/(2*sin(R));
    a1=D2*(MR(3,4)-MR(4,3));
    a2=D2*(MR(4,2)-MR(2,4));
    a3=D2*(MR(2,3)-MR(3,2));
end
rvec=[1,a1,a2,a3]';

if abs(cos(R))>=0.9999999999
    C1=MR(2,2)+MR(3,3);
    C2=MR(2,3)-MR(3,2);
    if abs(C1)<0.0000000001
        MR=MR*[1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 -1];
        rotation=0.5*acos((MR(2,2)+MR(3,3))/2);
        linreta=pi;
    end
    if C1<1.999999999
        if abs(C2)<0.0000000001
            MR=MR*[1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1];
            dum=MR(2,2)+MR(3,3)-1;
            linreta=acos(dum);
            if dum>=1
                linreta=0;
            end
            if dum<=-1
                linreta=pi;
            end
            rotation=pi/2;
        end
    end
end
orientation = 0.5*atan(real(rvec(3))/real(rvec(2)));
mr1=MR*inv((rotation)); %don't have fuction rota? what does this do?
%mr1=MR*inv(rota(rotation));

a1a=(mr1(3,4)-mr1(4,3));
a2a=(mr1(4,2)-mr1(2,4));
a3a=(mr1(2,3)-mr1(3,2));
rveca=[1,a1a,a2a,a3a]';
orientationa = 0.5*atan2(real(rvec(3)),real(rvec(2)));
mr2=inv((rotation))*MR;
%mr2=inv(rota(rotation))*MR;

a1b=(mr2(3,4)-mr2(4,3));
a2b=(mr2(4,2)-mr2(2,4));
a3b=(mr2(2,3)-mr2(3,2));
rvecb=[1,a1b,a2b,a3b]';
orientationb = 0.5*atan(real(rvec(3))/real(rvec(2)));
orientationc = 0.5*acos(MR(3,4)/sin(linreta));


if isnan(depol)==1;
    depol=0;
end

if isnan(linreta)==1
    linreta=0;
end

if isnan(orientation)==1
    orientation=0;
end

%return