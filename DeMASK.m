pth=input('Directory [default=C:\\User\\tir data\\yyyy\\New Folder]  ');
	if isempty(pth);
   	pth='C:\User\tir data\yyyy\New Folder';
    end
numbootstrap=input('Number bootstrap it  ');
	if isempty(numbootstrap);
   	numbootstrap=2;
    end       
cd(pth);
disp(pth);
H=dir;
[nf,dumn]=size(H);
dateNow=date;
AnalyzeDir=zeros(nf,1);
j=1;
l=1;

%OPEN input data files
for i=1:nf ;
   if H(i).isdir == 0;
      s=H(i).name;
      if s(end-5:end) == 'ON.csv';
         disp(s);
        fname=num2str(s);
        fname1=strsplit(fname,'.');
        fname2= fname1(1);
        class(fname2);
        fname2=char(fname2);
        fnamelist=strsplit(fname2,'_');
        dsDNA1=char(fnamelist(2));
        dsDNA=str2num(dsDNA1)/1000;  %DNA concentration determined from filenames
        fnameF=[char(fnamelist(2)) char(fnamelist(4))];
        [ fname '.csv'];
         fnameread=[ fname ];
        
        if j == 1 ;
            M1=dlmread(fnameread);
            figname1=fnameF;
            [d,y]=size(M1);
            X1=M1(:,1);
            Y1=M1(:,2);
            timedata1=X1;       %THIS IS X
            soldata1=Y1;          %THIS IS Y
             Ga = dsDNA;
        end
         if j == 2 ;
             M2=dlmread(fnameread);
            figname2=fnameF;
            [d,y]=size(M2);
            X2=M2(:,1);
            Y2=M2(:,2);
            timedata2=X2;       %THIS IS X
            soldata2=Y2;          %THIS IS Y
             Gb = dsDNA;
         end
          
         if j == 3 ;
             M3=dlmread(fnameread);
            figname3=fnameF;
            [d,y]=size(M3);
            X3=M3(:,1);
            Y3=M3(:,2);
            timedata3=X3;       %THIS IS X
            soldata3=Y3;          %THIS IS Y
             Gc = dsDNA;
         end
         
         if j == 4 ;
             M4=dlmread(fnameread);
            figname4=fnameF;
            [d,y]=size(M4);
            X4=M4(:,1);
            Y4=M4(:,2);
            timedata4=X4;       %THIS IS X
            soldata4=Y4;          %THIS IS Y
             Gd = dsDNA;
         end

      
        j=j+1;
      end
      
 
      
   end 
end

param0 = [0.13 0.13 1.1];    %INITAL GUESSES (See myresODEgloOnOffFull for more details)

Lb = [0 0 0]; % LOWER BOUNDS FOR MICROSCOPIC RATE FITTING
Glist= [Ga Gb Gc Gd];  % LIST OF DNA CONCENTRATIONS, FROM FILENAMES


[optparam  resnorm]=lsqnonlin(@myresODEgloOnOffFull,param0,Lb,[],[],timedata1,soldata1,timedata2,soldata2,timedata3,soldata3,timedata4,soldata4,Ga,Gb,Gc,Gd, options);
save(['KINETICS.dat'],'optparam','-ascii');
save(['Chi2.dat'],'resnorm','-ascii');
save(['KINETICSGorder.dat'],'Glist','-ascii');

%STORE BEST FIT LINE
      
      A=optparam(1);
      B=optparam(2);
      C=optparam(3);

      
      init1=0;
      init2=1;
      initvalue=[init1 init2 0 ];

      G=Ga;
%GENERATE FIGURES WITH BEST-FIT AND INPUT DATA (calls script myresODEgloFIGson)
 myresODEgloFIGson(A,B,C,Ga,timedata1,M1,initvalue,figname1,soldata1);
 
 G=Gb;
 myresODEgloFIGson(A,B,C,Gb,timedata2,M2,initvalue,figname2,soldata2);

 G=Gc;
 myresODEgloFIGson(A,B,C,Gc,timedata3,M3,initvalue,figname3,soldata3);

 G=Gd;
 myresODEgloFIGson(A,B,C,Gd,timedata4,M4,initvalue,figname4,soldata4);
 
 close all
 
 %BOOTSTRAP
 numit=0
while ((numit) + 0) < numbootstrap ;
 [t,soldata1b]= ode15s(@(t,y) smodefcnON(t,y,A,B,C,Ga), timedata1, initvalue);
 soldata1fb= (1-(soldata1b(:,1)+soldata1b(:,2)))*100 ;
 residualbootstrap=soldata1-soldata1fb;
 numdata=1;
 
 while (numdata-1) < length(residualbootstrap);
     rand=randi(length(residualbootstrap)); %randomly select residuals (with replacement)
     soldata1bootstrapi=soldata1fb(numdata)+residualbootstrap(rand); %add residual to best-fit
     if numdata == 1;
         soldata1bootstrap=soldata1bootstrapi;
     end
     if numdata >1;
    soldata1bootstrap=cat(1,soldata1bootstrap,soldata1bootstrapi); 
     end
    numdata=numdata+1;
 end
 
[t,soldata2b]= ode15s(@(t,y) smodefcnON(t,y,A,B,C,Gb), timedata2, initvalue);
 soldata2fb= (1-(soldata2b(:,1)+soldata2b(:,2)))*100 ;
 residualbootstrap2=soldata2-soldata2fb;
 numdata=1;
 while (numdata-1) < length(residualbootstrap2);
     rand=randi(length(residualbootstrap2)); %randomly select residuals (with replacement)
     soldata2bootstrapi=soldata2fb(numdata)+residualbootstrap2(rand); %add residual to best-fit
     if numdata == 1;
         soldata2bootstrap=soldata2bootstrapi;
     end
     if numdata >1;
    soldata2bootstrap=cat(1,soldata2bootstrap,soldata2bootstrapi); 
     end
    numdata=numdata+1;
 end
 
 
 [t,soldata3b]= ode15s(@(t,y) smodefcnON(t,y,A,B,C,Gc), timedata3, initvalue);
 soldata3fb= (1-(soldata3b(:,1)+soldata3b(:,2)))*100 ;
 residualbootstrap3=soldata3-soldata3fb;
 numdata=1;
 while (numdata-1) < length(residualbootstrap3);
     rand=randi(length(residualbootstrap3)); %randomly select residuals (with replacement)
     soldata3bootstrapi=soldata3fb(numdata)+residualbootstrap3(rand); %add residual to best-fit
     if numdata == 1;
         soldata3bootstrap=soldata3bootstrapi;
     end
     if numdata >1;
    soldata3bootstrap=cat(1,soldata3bootstrap,soldata3bootstrapi); 
     end
    numdata=numdata+1;
 end
 
  [t,soldata4b]= ode15s(@(t,y) smodefcnON(t,y,A,B,C,Gd), timedata4, initvalue);
 soldata4fb= (1-(soldata4b(:,1)+soldata4b(:,2)))*100 ; 
 residualbootstrap4=soldata4-soldata4fb;
 numdata=1;
 while (numdata-1) < length(residualbootstrap4);
     rand=randi(length(residualbootstrap4)); %randomly select residuals (with replacement)
     soldata4bootstrapi=soldata4fb(numdata)+residualbootstrap4(rand); %add residual to best-fit
     if numdata == 1;
         soldata4bootstrap=soldata4bootstrapi;
     end
     if numdata >1;
    soldata4bootstrap=cat(1,soldata4bootstrap,soldata4bootstrapi); 
     end
    numdata=numdata+1;
 end

%FIT BOOTSTRAPPED DATASET
[optparam  resnorm]=lsqnonlin(@myresODEgloOnOffFull,param0,Lb,[],[],timedata1,soldata1bootstrap,timedata2,soldata2bootstrap,timedata3,soldata3bootstrap,timedata4,soldata4bootstrap,Ga,Gb,Gc,Gd);
Aboot=optparam(1);
Bboot=optparam(2);
KeqFboot=Aboot/Bboot;
optparamfinal=cat(2,optparam,KeqFboot);
optparamfinal=cat(2,optparamfinal,resnorm);
if numit==0;
    save(['KINETICSbootstrap.dat'],'optparamfinal','-ascii');
    save(['Chi2bootstrap.dat'],'resnorm','-ascii');
else
save(['KINETICSbootstrap.dat'],'optparamfinal','-append','-ascii');
save(['Chi2bootstrap.dat'],'resnorm','-append','-ascii');
end


numit=numit+1  %Iterate bootstrap and repeat loop


end
close all