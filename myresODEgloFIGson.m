function work=myresODEgloFIGson(A,B,C,G,timedata1,M1,initvalue,figname1,soldata)
[t,y] = ode15s(@(t,y) smodefcnON(t,y,A,B,C,G), timedata1, initvalue);

fig=figure;
highFret= (1-(y(:,1)+y(:,2)))*100; %USE FOR ON RATES
lowFret= y(:,3)*100;   
resid=soldata-highFret;
timemax=30;
%USE timemax=max(t)TO SHOW ALL DATA (NOT ONLY FROM 0 TO 30 SECONDS)
%timemax=max(t);   
plot(t,highFret,'k-',M1(:,1),M1(:,2),'b.','linewidth',3);
temp=axis;
temp(3)=-5;
temp(4)=110;
temp(1)=0;
temp(2)=timemax;
axis(temp);
zoom on;
title(['  KINETICS ' figname1]);
filename1=['TRACE_' figname1 '.pdf'];
print(filename1,'-dpdf');


        fig=figure;
        plot(t,resid,'b.');
        temp=axis;
        temp(3)=min(resid);
        temp(4)=max(resid);
        if (max(resid)-min(resid)) < 1.0;
            temp(3)=-1.0;
            temp(4)=1.0;
        end
        temp(1)=0;
        temp(2)=timemax;
        axis(temp);
        zoom on;
        title(['  Residuals ' figname1]);
        filename1=['Residuals_' figname1 '.pdf'];
        print(filename1,'-dpdf');


fig=figure;
plot(t,y(:,1),'r-',t,y(:,2),'y-',t,y(:,3),'c-','linewidth',4);
temp=axis;
temp(1)=0;
temp(2)=timemax;
temp(3)=-0.05;
temp(4)=1.05;
axis(temp);
zoom on;
title(['  KINETICS ' figname1]);
filename2=['FULL_' figname1 '.pdf'];
print(filename2,'-dpdf');
work = 1;