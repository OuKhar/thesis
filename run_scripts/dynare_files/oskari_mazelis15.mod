% dsge model with shadow banking

% variables 

var y c pi inv k q int rd rb p lab In R Lambda varho g Z pi_dot Dp delta U F ym phi w
Nn z nw Nex x ene veta sb d	b ssb rsb v qa fs Theta shock_A shock_FS shock_G shockvar_Lambda shockvar_Ne shock_XI shock_MU;

% exogenous variables (shock sigmas)

varexo e_A e_FS e_G e_Lambda e_Ne e_I e_XI e_MU;

% parameters

parameters cbeta ch cchihh cxi cfrisch comega ctheta calpha cdelta cepsilon cgamma cgammap ckappapi ckappay crooi  
czeta cs cb comegahh cchisb ceta I_ss g_ss rho_A rho_FS rho_XI rho_MU rho_G rho_Lambda LambdaBAR rho_Ne; 

% parameter values

% households, labour
cbeta = 0.99;		    % discount rate
ch = 0.815;		        % habit formation paramateres
cfrisch = 0.276;	    % inverse frisch elasticity of labor supply
ctheta = 0.972;		    % the survival probability
calpha = 0.33;		    % capital share
cdelta = 0.02071939;	% depreciation rate
ceta =  1.728;		    % elasticity of investment adjustment cost

cb = 0.0376;		    % depreciation rate parameter
czeta = 7.2;		    % elasticity of marginal depreciation respect to utilization rate

% monetary policy
crooi = 0.8;		    % interest rate smoothing parameter
ckappapi = 1.5;		    % inflation coefficent
ckappay = -0.125;	    % output gap coefficent

% retail firms
cepsilon = 4.167;	    % elasticity of substition between goods
cgamma = 0.779;		    % calvo parameter
cgammap = 0.241;	    % price indexation parameter

% endogenous determined parameters 
cchihh = 3.409;		    % starting value for the labour utility weight
comega = .002;		    % starting value of proportional starting up funds

% shadow bank specific parameters
comegahh = 0.02;	    % divertibility, household bargaining power
cs = 4.6; 		        % fund matching efficiency
cxi = 0.87;		        % fund matching elasticity 
cchisb = 0.1;           % shadow bank sepratation rate, 

I_ss = 0.14153927;      % steady state investment
g_ss = 0.16975710;      % steady state government spending

% shock persistencies 
rho_Lambda = 0.61;	    % source: M15 posterior
rho_A = 0.95;		    % source: GK11
rho_XI = 0.66;	        % source: GK11 
rho_G = 0.95;		    % source: GK11
rho_Ne = 0.16;	        % source: M15 posterior
rho_FS = 0.76;	        % source: M15 posterior
rho_MU = 0.63;	        % source: M15 posterior
LambdaBAR = 0.381;	    % source: M15 posterior

% define the model

model;

% BASE MODEL EQUATIONS

% 1. Production output

exp(ym) = exp(shock_A)*((exp(shock_XI)*exp(k(-1)))^calpha)*exp(lab)^(1-calpha);

% 2. Euler equation of household consumption  

cbeta*exp(R)*exp(Lambda(+1)) = 1;					

% 3. Stochastic discount rate 

exp(Lambda) = exp(varho)/exp(varho(-1));

% 4. Marginal utility of consumption

exp(varho) = (exp(c)-ch*exp(c(-1)))^(-1)-cbeta*ch*(exp(c(+1))-ch*exp(c))^(-1); 

% 5. labour market clearance/equilibrium from FOC of labour (household and goods producer)

cchihh*exp(lab)^cfrisch = exp(varho)*exp(p)*(1-calpha)*exp(ym)/exp(lab);

% 6. FOC of capital (goods producer), return to capital

exp(rb) = (exp(p)*calpha*exp(ym)/exp(k(-1))+(exp(q)-exp(delta)))/exp(q(-1));

% 7. Fisher relation

exp(int) = exp(rd)*exp(pi(+1));

% 8. definition of In, net investment

In = exp(inv)*exp(shock_MU)-exp(delta)*exp(shock_XI)*exp(k(-1));

% 9. Capital & investment k(+1)

exp(k) = exp(shock_XI)*exp(k(-1))+In;

% 10. price of capital, optimal investment decision

exp(q) = 1+ceta/2*((In+I_ss)/(In(-1)+I_ss)-1)^2+
    ceta*((In+I_ss)/(In(-1)+I_ss)-1)*(In+I_ss)/(In(-1)+I_ss)-
    cbeta*exp(Lambda(+1))*ceta*((In(+1)+I_ss)/(In+I_ss)-1)
    *((In(+1)+I_ss)/(In+I_ss))^2;

% 11. resource & policy, Taylor rule

exp(int) = exp(int(-1))^crooi*((1/cbeta)*exp(pi)^ckappapi*((1/exp(p))/(cepsilon/(cepsilon-1)))^(ckappay))^(1-crooi)*exp(e_I);

% 12. government consumption

exp(g) = g_ss*exp(shock_G);

% 13. aggregate resource contraint

exp(y) = exp(c)+exp(v)+exp(inv)+(ceta/2)*((In+I_ss)/(In(-1)+I_ss)-1)^2*(In+I_ss)+exp(g);

% 14. optimal price choice

exp(F) = exp(y)*exp(p)+cbeta*cgamma*exp(Lambda(+1))*exp(pi(+1))^cepsilon*exp(pi)^(-cepsilon*cgammap)*exp(F(+1));

% 15. optimal price choice

exp(Z) = exp(y)+cbeta*cgamma*exp(Lambda(+1))*exp(pi(+1))^(cepsilon-1)*exp(pi)^(cgammap*(1-cepsilon))*exp(Z(+1));

% 16. Optimal price choice

exp(pi_dot) = cepsilon/(cepsilon-1)*exp(F)/exp(Z)*exp(pi);

% 17. definition of inflation, price index (Calvo), the New Keynesian Philips Curve

(exp(pi))^(1-cepsilon) = cgamma*exp(pi(-1))^(cgammap*(1-cepsilon))+(1-cgamma)*exp(pi_dot)^(1-cepsilon);

% 18. price dispersion

exp(Dp) = cgamma*exp(Dp(-1))*exp(pi(-1))^(-cgammap*cepsilon)*exp(pi)^cepsilon
+(1-cgamma)*((1-cgamma*exp(pi(-1))^(cgammap*(1-cgamma))*exp(pi)^(cgamma-1))/(1-cgamma))^(-cepsilon/(1-cgamma));

% 19. Depreciation rate
exp(delta) = cdelta+cb/(1+czeta)*exp(U)^(1+czeta);

% 20. Optimal capacity utilization rate

exp(p)*calpha*exp(ym)/exp(U) = cb*exp(U)^czeta*exp(shock_XI)*exp(k(-1));

% 21. arbitrage 

cbeta*exp(Lambda(+1))*exp(rb(+1))= cbeta*exp(Lambda(+1))*exp(rd);

% 22. wholesale, the retailers' output

exp(ym) = exp(y)*exp(Dp);

% 23. wages

exp(w) = exp(p)*(1-calpha)*exp(ym)/exp(lab);

% BANK SECTOR EQUATIONS

% 24. optimal leverage ratio

exp(phi) = exp(ene)/(shockvar_Lambda-exp(veta));

% 25. aggregate net worth of banks

exp(nw) = exp(Nex)+exp(Nn);

% 26. Net worth of new banks

exp(Nn) = comega*exp(q)*exp(k(-1));

% 27. Net worth of existing banks

exp(Nex) = exp(shockvar_Ne)*ctheta*exp(z)*exp(nw(-1));

% 28. Growth rate of banks' capital

exp(z) = (exp(rb)-exp(rd(-1)))*exp(phi(-1))+exp(rd(-1));

% 29. Growth rate of banks' net wealth

exp(x) = (exp(phi)/exp(phi(-1)))*exp(z);

% 30. definition of ene, value of banks' net wealth

exp(ene) = (1-ctheta)+cbeta*exp(Lambda(+1))*ctheta*exp(z(+1))*exp(ene(+1));

% 31. definition of veta, value of banks' capital

exp(veta) = (1-ctheta)*cbeta*exp(Lambda(+1))*(exp(rb(+1))-exp(rd))+cbeta*exp(Lambda(+1))*ctheta*exp(x(+1))*exp(veta(+1));

% 32. Loan portfolio

exp(sb) = (exp(phi)*exp(nw))/exp(q);

% 33. banks balance sheet

exp(d) = exp(q(-1))*exp(sb(-1))-exp(nw(-1));

% 34. weighed interest rate (from holdings of deposits and holdings of fund shares)

exp(R) = exp(rd)*(exp(d)/exp(b))+exp(rsb)*(exp(fs)/exp(b));

% SHADOW BANK EQUATIONS

% 35. definition of fund shares (fs) plus shock (FS)

exp(fs) = (1-cchisb)*exp(fs(-1))+exp(qa)*exp(v)+ shock_FS;

% 36. Euler condition for fund adverts (qa) 

1/exp(qa) = cbeta*exp(Lambda(+1))*((exp(rb)-exp(rsb))+(1-cchisb)*(1/exp(qa(+1))));

% 37. financial intermediaries profits to creditors from (k+1)
exp(k(+1)) = exp(sb) + exp(ssb);						
% TAI:
% balance sheet constraint of shadow banks 

%exp(ssb) = exp(fs)/exp(q);

% 38. interest rate for shadow banks 

exp(rsb) = (1-comegahh)*exp(rd)+comegahh*(exp(rb)+exp(Theta));

% 39. probability of shadow bank finding suitable funding 

exp(qa) = cs*(exp(Theta))^(-cxi);

% 40. Matching function of shadow banks

exp(Theta) = exp(v)/(exp(d(+1))-exp(fs));

% 41. savings definition

exp(b) = exp(fs)+exp(d);

% SHOCKS

% 42. shock A (TFP shock)

shock_A = rho_A*shock_A(-1)+e_A;

% 43. shock XI (capital quality shock)

shock_XI = rho_XI*shock_XI(-1)+e_XI;

% 44. shock G (government spending shock)

shock_G = rho_G*shock_G(-1)+e_G;

% 45. shock-variable Ne (existing bank net worth)

shockvar_Ne = rho_Ne*shockvar_Ne(-1)+e_Ne;

% 46. shock FS (shock to shadow bank fund shares)

shock_FS = rho_FS*shock_FS(-1)+e_FS;

% 47. shockvariable Lambda (divertible amount in banks' optimal leverage)

shockvar_Lambda = (1-rho_Lambda)*LambdaBAR+rho_Lambda*shockvar_Lambda(-1)+e_Lambda;

% 48. shock MU (investment shock) 

shock_MU = rho_MU*shock_MU(-1)+e_MU;

end;

% initial values for the steady state

initval;


y      		= -0.124503;
c      		= -0.587963;
pi     		= -1.66161e-15;
inv    		= -1.84693;
k      		= 1.84195;
q      		= 0;
int    		= 0.0100503;
rd     		= 0.0100503;
rb     		= 0.0100503;
p      		= -0.274412;
lab    		= -1.09305;
In     		= 0;
R      		= 0.0100494;
Lambda 		= 0;
varho  		= 0.631075;
g      		= -1.77339;
Z      		= 1.35045;
pi_dot 		= -6.1351e-15;
Dp     		= -2.98979e-15;
delta  		= -3.68888;
U      		= -0.00838709;
F      		= 1.07604;
ym     		= -0.124503;
phi    		= 0.964136;
w           = 0.293662;
Nn     		= -4.37266;
z      		= 0.0100503;
nw     		= -0.365326;
Nex    		= -0.383675;
x      		= 0.0100503;
ene    		= 3.71995e-06;
veta   		= -8.58873;
sb     		= 0.59881;
d      		= 0.118655;
b      		= 0.775746;
ssb    		= 1.5016;
rsb    		= 0.0100485;
v      		= -12.1798;
qa     		= 9.92237;
fs     		= 0.0451926;
Theta  		= -9.65093;

shock_A  		= 0;
shock_FS 		= 0;
shock_G         = 0;
shockvar_Lambda = 0.381;
shockvar_Ne     = 0;
shock_XI        = 0;
shock_MU        = 0;

end;

shocks;

var e_A = 0.01;
var e_FS = 0;
var e_G = 0;
var e_Lambda = 0;
var e_I = 0.01;
var e_Ne = 0;   
var e_XI = 0.01;
var e_MU = 0.01;

end;

resid(1);

steady;

check(qz_zero_threshold=1e-10);

stoch_simul(order=1, irf=40, nograph) c k y pi lab b int d fs rsb ssb U delta sb q inv rb rd Nn nw Nex R varho;