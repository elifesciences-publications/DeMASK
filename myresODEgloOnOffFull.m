 function residual=myresODEgloOnOffFull(param,timedata1,soldata1,timedata2,soldata2,timedata3,soldata3,timedata4,soldata4, Ga, Gb, Gc, Gd, options)
   

      A=param(1);  %parameter A to be fit in least squares fitting
      B=param(2);  %parameter B to be fit in least squares fitting
      C=param(3);  %parameter C to be fit in least squares fitting

      
      init1=0;
      init2=1;
      initvalue=[init1 init2 0 ];

      %USE smodefcnON TO GENERATE FUNCTION TO COMPARE WITH INPUT DATA AND
      %GENERATE RESIDUALS
       [t,odesol1] = ode15s(@(t,y) smodefcnON(t,y,A,B,C,Ga), timedata1, initvalue);
       odesolf1= (1-(odesol1(:,1)+odesol1(:,2)))*100 ;
        size(odesolf1);
        [weight1,weight2]=size(timedata1);
       residual1=(soldata1-odesolf1)/weight1;

     
      %USE smodefcnON TO GENERATE FUNCTION TO COMPARE WITH INPUT DATA AND
      %GENERATE RESIDUALS
       [t,odesol2] = ode15s(@(t,y) smodefcnON(t,y,A,B,C,Gb), timedata2, initvalue);
       odesolf2= (1-(odesol2(:,1)+odesol2(:,2)))*100 ;
        size(odesolf2);
        [weight1,weight2]=size(timedata2);
       residual2=(soldata2-odesolf2)/weight1;
       
       residualv1=cat(1,residual1,residual2);  %APPEND THE RESIDUALS FROM DATASET2 TO THE RESIDUALS FROM PREVIOUS DATASETS


      %USE smodefcnON TO GENERATE FUNCTION TO COMPARE WITH INPUT DATA AND
      %GENERATE RESIDUALS
       [t,odesol3] = ode15s(@(t,y) smodefcnON(t,y,A,B,C,Gc), timedata3, initvalue);
       odesolf3= (1-(odesol3(:,1)+odesol3(:,2)))*100 ;
        size(odesolf3);
         [weight1,weight2]=size(timedata3);
       residual3=(soldata3-odesolf3)/weight1;
       
       residualv2=cat(1,residualv1,residual3);   %APPEND THE RESIDUALS FROM DATASET3 TO THE RESIDUALS FROM PREVIOUS DATASETS
       

       [t,odesol4] = ode15s(@(t,y) smodefcnON(t,y,A,B,C,Gd), timedata4, initvalue);
       odesolf4= (1-(odesol4(:,1)+odesol4(:,2)))*100 ;
        size(odesolf4);
         [weight1,weight2]=size(timedata4);
       residual4=(soldata4-odesolf4)/weight1;
       
       residual=cat(1,residualv2,residual4);  %GENERATE FINAL LIST OF RESIDUALS
