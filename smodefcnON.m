function dydt = smodefcnON(t,y,A,B,C,G)
dydt = zeros(3,1);
dydt(1) = -A*y(1)+B*y(2);
dydt(2) = -B*y(2)+A*y(1)-G*C*y(2);
dydt(3)= G*C*y(2);