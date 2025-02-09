function matrix = f_TVFF_1(pic,uilk,m,n,lh,eps)
t=0.1;
h=1;
alfa=0.5;
kez=1;
err=1;
r=5000;
for i=1:kez
   b(i)=((i+1)^(1-alfa))-(i^(1-alfa));
end
 mu=(t^alfa)*gamma(2-alfa);
n1=1;
n2=2;      
etki1=1;
etki2=2;
alfa1=0.3;
alfa2=1.5;
beta1=n1-alfa1;
beta2=n2-alfa2;
cons_11= 1^(beta1)/gamma(beta1+2);
cons_21=(lh-1)^(beta1+1)-(lh-beta1-1)*lh^beta1;
cons_12= 1^(beta2)/gamma(beta2+2);
cons_22=(lh-1)^(beta2+1)-(lh-beta2-1)*lh^beta2;
gama(1)=cons_21; gama(lh+1)=1;
teta(1)=cons_22; teta(lh+1)=1;
for i=1:(lh-1)
   gama(i+1)= [(lh-i+1)^(beta1+1)-2*(lh-i)^(beta1+1)+(lh-i-1)^(beta1+1)];
   teta(i+1)= [(lh-i+1)^(beta2+1)-2*(lh-i)^(beta2+1)+(lh-i-1)^(beta2+1)];
end
son_carpan1=gamma(beta1+1)*0.5*(lh^(alfa1-1));
son_carpan2=gamma(beta2+1)*0.5*(lh^(alfa2-2));

for i=1:lh+1
    gama(i)=gama(i)*cons_11*son_carpan1;
    teta(i)=teta(i)*cons_12*son_carpan2;
    e1(i)=gama(i)/(etki1+1);e1(2*(lh+1)-i)=e1(i);
    e2(i)=teta(i)/(etki2+1);e2(2*(lh+1)-i)=e2(i);
end
e1(lh+1)=2*(gama(lh+1)+etki1*0.5)/(etki1+1);
e2(lh+1)=2*(teta(lh+1)+etki2*0.5)/(etki2+1);
e11=zeros(1,2*lh+1);e111=zeros(1,2*lh+1);
e22=zeros(1,2*lh+1);e222=zeros(1,2*lh+1);
e11(lh+1)=e1(lh+1);e22(lh+1)=e2(lh+1);

for z=lh+2:2*lh+1
    e11(z)=2*e1(z);   
    e22(z)=2*e2(z);
end
e111(1:lh)=2*e1(1:lh); e111(lh+1)=e1(lh+1);
e222(1:lh)=2*e2(1:lh); e222(lh+1)=e2(lh+1);
exy=[0.04 0.92 0.04];
exy1=[0 0.92 0.08];
exy11=[0.08 0.92 0];
UUUU(:,:,1)=uilk;
snr_RC2_kesirli=-10+zeros(1,r+1);
snr_abs_kesirli=zeros(1,r+1);
% kk=0;
% while err>eps
% kk=kk+1
for kk=1:r
Ruson=UUUU(:,:,kk);
u00= [Ruson(1,:); Ruson(1,:); Ruson; Ruson(m,:);Ruson(m,:)];
u000=[u00(:,1) u00(:,1) u00 u00(:,n) u00(:,n)];
matris=u000;

RX=zeros(m,n); RY=zeros(m,n); 
RXX=zeros(m,n); RYY=zeros(m,n); RXY=zeros(m,n);
for i=1:m
   for j=1:n
      RX(i,j)=(matris(i+3,j+2)-matris(i+1,j+2))/(2*h);
      RY(i,j)=(matris(i+2,j+3)-matris(i+2,j+1))/(2*h);
      RXX(i,j)=(matris(i+3,j+2)-2*matris(i+2,j+2)+matris(i+1,j+2))/(h*h);
      RYY(i,j)=(matris(i+2,j+3)-2*matris(i+2,j+2)+matris(i+2,j+1))/(h*h);
RXY(i,j)=(matris(i+3,j+3)-matris(i+1,j+3)-matris(i+3,j+1)+matris(i+1,j+1))/(4*h*h);
   
   end
end

Rtopb=zeros(m,n);
if kk<=kez indis=kk-1; 
        else
          indis=kez; 
        end
for i=1:indis
    Rtopb=Rtopb+b(i)*(UUUU(:,:,kk-i+1)-UUUU(:,:,kk-i));end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  x_1
Alast_x=imfilter(RX,e1');
Blast_x=imfilter(RX,e11');
Clast_x=imfilter(RX,e111');

last_x(lh+1:m-lh,:)=Alast_x(lh+1:m-lh,:);
last_x(1:lh,:)=Blast_x(1:lh,:);
last_x(m-lh+1:m,:)=Clast_x(m-lh+1:m,:);

% % y_1

Alast_y=imfilter(RY,e1);
Blast_y=imfilter(RY,e11);
Clast_y=imfilter(RY,e111);
last_y(:,lh+1:n-lh)=Alast_y(:,lh+1:n-lh);
last_y(:,1:lh)=Blast_y(:,1:lh);
last_y(:,n-lh+1:n)=Clast_y(:,n-lh+1:n);



%%%%%%%%%%% xx_1
Alast_xx=imfilter(RXX,e2');
Blast_xx=imfilter(RXX,e22');
Clast_xx=imfilter(RXX,e222');
last_xx(lh+1:m-lh,:)=Alast_xx(lh+1:m-lh,:);
last_xx(1:lh,:)=Blast_xx(1:lh,:);
last_xx(m-lh+1:m,:)=Clast_xx(m-lh+1:m,:);
%%%%%%%%%%% yy_1
Alast_yy=imfilter(RYY,e2);
Blast_yy=imfilter(RYY,e22);
Clast_yy=imfilter(RYY,e222);
last_yy(:,lh+1:n-lh)=Alast_yy(:,lh+1:n-lh);
last_yy(:,1:lh)=Blast_yy(:,1:lh);
last_yy(:,n-lh+1:n)=Clast_yy(:,n-lh+1:n);
%%%%%%%%%%% xy_1
LAST_XY=RXY;
ALAST_XY=0.5*(imfilter(RXY,exy)+imfilter(RXY,exy'));
B1LAST_XY=imfilter(RXY,exy1);
C1LAST_XY=imfilter(RXY,exy11);
B2LAST_XY=imfilter(RXY,exy1');
C2LAST_XY=imfilter(RXY,exy11');

LAST_XY(2:m-1, 2:n-1)=ALAST_XY(2:m-1, 2:n-1);
LAST_XY(2:m-1,1)=B1LAST_XY(2:m-1,1);
LAST_XY(2:m-1,n)=C1LAST_XY(2:m-1,n);
LAST_XY(1,2:n-1)=B2LAST_XY(1,2:n-1);
LAST_XY(m,2:n-1)=C2LAST_XY(m,2:n-1);
LAST_XY=LAST_XY/2;
for i=1:m
    for j=1:n
ppay(i,j)=last_xx(i,j)*(last_y(i,j)^2+1)-2*last_x(i,j)*last_y(i,j)*LAST_XY(i,j)+last_yy(i,j)*(last_x(i,j)^2+1);
ppayda(i,j)=((last_x(i,j)^2)+(last_y(i,j)^2)+1)^(3/2);

  son_matris(i,j)=ppay(i,j)/ppayda(i,j);  
    
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r1=UUUU(:,:,kk)-Rtopb+mu*son_matris;

r3=uint8(uilk);
rsonuc=uilk;

for i=1:m
    for j=1:n
if (r3(i,j)==0 || r3(i,j)==255) rsonuc(i,j)=r1(i,j); else;

end 
    end
end
UUUU(:,:,kk+1)=rsonuc;
A= UUUU(:,:,kk);
B= UUUU(:,:,kk+1);
norm_A = norm(A, 'fro');
norm_diff = norm(A - B, 'fro');
err = round(norm_diff / norm_A, 5);
snr_abs_kesirli(kk+1)=err;
snr_RC2_kesirli(kk+1)=psnr(uint8(rsonuc), uint8(pic));
% if kk+1==1250 printf("Change value of eps");break; end
if snr_RC2_kesirli(kk+1)<snr_RC2_kesirli(kk) break; end
end
matrix=UUUU(:,:,kk+1);
   end

   