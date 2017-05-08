Ts=1/128;
t=0:Ts:10;
num=[1 0.3 0.4];
den=[1 1.2 0.5 0.3 0.2];
G1=tf(num,den,Ts);
u=chirp(t,0,1,3);
u=u';
[y,t]=lsim(G1,u,t);

% 4 is the order , can be get in  https://github.com/bjlxwjs/transfer-function-Order-Estimation
n=4;

sys=dttfe(y,u,Ts,n)



  %           sys =
 
  %  -1.298e-16 z^3 + z^2 + 0.3 z + 0.4
  % -------------------------------------
  % z^4 + 1.2 z^3 + 0.5 z^2 + 0.3 z + 0.2
