
%%project1_REE430
%%% NACA 24012 airfoil
clc;
clear;

%% Defining the constants

Vinf=10;
CL=(2*3)/20; %%2 represents the first digit of airfoil name %%lift coeff
alpha=CL/(2*pi);
p=4/20; %%4 represents the second digit of airfoil name %%x of max chamber
t=12/100; %%12 represents the forth and fifth digits of airfoil name %%max thickness
k1=6.643;
r= 0.29;
x=0:0.01:1;

   for i=1:length(x) 
     theta=acos(1-x.*2);
     theta_m=acos(1-2*r);
     
   %% Chamber line equations
     a0=0;
     a1=(k1*r^2/6)*(3-r);
     a2= (-1)*3*r*k1/6;
     a3= k1/6;
     b0 = k1*r^3/6;
     b1 = (-1)*k1*r^3/6;
     b2=0;
     b3=0;
     t10=a1+a2+9*a3/8;
     t11=-a2-3*a3/2;
     t12=3*a3/8;
     t20=b1+b2+9*b3/8;
     t21=(-1)*b2-3*b3/2;
     t22=3*b3/8;
    %I0, In, where n=1,2
     I0=pi*t20+theta_m*(t10-t20)+sin(theta_m)*(t11-t21)+(sin(2*theta_m))*(t12-t22)/2;
  I1=pi*t21/2+(theta_m)*(t11-t21)/2+(sin(theta_m))*(t10-t20+t12/2-t22/2)+(sin(2*theta_m))*(t11-t21)/4+(sin(3*theta_m))*(t12-t22)/6;
  I2=pi*t22/2+(theta_m)*(t12-t22)/2+(sin(theta_m))*(t11-t21)/2+(sin(2*theta_m))*(t10-t20)/2+(sin(3*theta_m))*(t11-t21)/6+(sin(4*theta_m))*(t12-t22)/8;
    %A0, An, where n=1,2
     A0=alpha-(I0/pi);
     A1=2*I1/pi;
     A2=2*I2/pi;
    %%Getting Gamma to find the upper and lower velocities
     Gamma=2*Vinf*(A0*(1+cos(theta))./sin(theta)+(sin(theta).*A1)+(sin(2*theta).*A2));
     Vu=Vinf+Gamma./2;
     VL=Vinf-Gamma./2;
     Cp_u=1-(Vu/Vinf).^2;
     Cp_L=1-(VL/Vinf).^2;
     
     %% Airfoil coordinates
    if x(i)<=p
         yc(i)= (k1/6)*(((x(i)).^3)-(3*r*(x(i)).^2)+(r^2)*((x(i)).*(3-r)));
         dyc(i)=(k1/6)*(((x(i)).*3)^2-(6*x(i).*r)+(r^2)*(3-r)); 
          
     elseif x(i) > p 
         yc(i)= ((k1*r^3)/6)*(1-x(i));
         dyc(i) = (-1*k1*r^3)/6;
    end
    %%thickness
    yt =5*t*(sqrt(x).*.2969-x.*.1260-.3516*(x).^2+.2843*(x).^3-0.1036*(x).^4); 
    theta_airfoil=atan(dyc(i));
    
    %coordinates of upper and lower surfaces
    Xu(i)=x(i)-yt(i).*sin(theta_airfoil);
    Yu(i)=yc(i)+yt(i).*cos(theta_airfoil);
    Xl(i)=x(i)+yt(i).*sin(theta_airfoil);
    Yl(i)=yc(i)-yt(i).*cos(theta_airfoil);
   end
    
  %% Plotting Cp vs. x and Plotting the airfoil
 hold on
 grid on
 plot(Xu,Yu,Xl,Yl,'LineWidth',2);
 plot(x,yc(i));
 plot(x,Cp_u,'LineWidth',2,'color','k');
  plot(x,Cp_L,'LineWidth',2,'color','m');
  title('Cp versus position for NACA 24012 Airfoil');
  ylabel('Cp');
  xlabel('x');
  hold off