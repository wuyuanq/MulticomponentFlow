% This is the input file of two-component and single-phase flow. The up-rihgt 
% corner (production point) and down-left corner (injection point) are open
% and the other boundaries are closed. At the initial time, there is C3H8 in 
% the well.  We impose pressures on the production point and inject CH4 into the
% injection point. The temperature is constant in the whole process. 

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on October 23rd, 2012

path('/Users/wuy/Research/MulticomponentFlow_matlab', path);

model.numDims = 2;  
model.Nc = 2;
model.phase = 'g';

Lx = 4; % unit: meter  
Ly = 4;
timeEnd = 0.1*365*24*3600;
nx = 40;  
ny = 40; 
nt = 0.1*365*24;
model.nx = nx;   
model.ny = ny;   
model.nt = nt;
model.xs = (0:nx)*Lx/nx; 
model.ys = (0:ny)*Ly/ny; 
model.ts = (0:nt)*timeEnd/nt; 

K_const = 9.869233*1.D-15; % unit: m^2 
model.Kxx = zeros(model.nx, model.ny);  
model.Kxx(1:end, 1:end) = K_const; 
model.Kyy = zeros(model.nx, model.ny);  
model.Kyy(1:end, 1:end) = K_const; 

model.poro = zeros(model.nx, model.ny);  
model.poro(1:end, 1:end) = 0.2;

model.T = 480; % unit: K

model.gravX = 0;
model.gravY = -9.807; % m/s^2

model.ct = zeros(model.Nc, 1); % critical T
model.cp = zeros(model.Nc, 1); % critical pressure
model.cv = zeros(model.Nc, 1); % critical volume
model.af = zeros(model.Nc, 1); % acentric factor
model.mw = zeros(model.Nc, 1);
model.src = zeros(model.Nc, model.nx, model.ny); 

% the first component is methane
model.ct(1) = 190; % unit: K
model.cp(1) = 4.6*1.D6; % unit: Pa
model.cv(1) = 0.0062; % unit: m^3/kg
model.af(1) = 0.01;
model.mw(1) = 0.016; % unit: kg/mol

% the second component is propane
model.ct(2) = 370; % unit: K
model.cp(2) = 4.2*1.D6; % unit: Pa
model.cv(2) = 0.0045; % unit: m^3/kg
model.af(2) = 0.15;
model.mw(2) = 0.044; 

% the binary interaction parameters
% reference: Page 155 in A. Firoozabadi's book   
model.delta = zeros(model.Nc);
model.delta(1,2) = 0.036;
model.delta(2,1) = model.delta(1,2);

model.isDiriX = zeros(2, model.ny);   
model.isDiriY = zeros(model.nx, 2);   
model.isDiriX(2,end) = 1;
model.isDiriY(end,2) = 1;

model.xBdryX = zeros(model.Nc, 2, model.ny); 
model.xBdryY = zeros(model.Nc, model.nx, 2);  
model.xBdryX(1, 1, 1) = 1.0;
model.xBdryY(1, 1, 1) = 1.0;

model.xInit = zeros(model.Nc, model.nx, model.ny);
model.xInit(1, 1:end, 1:end) = 0.0;
model.xInit(2, 1:end, 1:end) = 1.0; 

uconst = 2*1.D-6; % unit: m/s
model.uBdryX = zeros(2, model.ny); 
model.uBdryY = zeros(model.nx, 2);  
model.uBdryX(1,1) = -uconst;
model.uBdryY(1,1) = -uconst;

% load initial pressures
initPres = zeros(1,ny);
temp = load('../../InitialReservoirData/P1.txt');
k = 1;
for i = 1 : ny
    initPres(i) = temp(k);
    k = k + 1;
end

model.pBdryX = zeros(2, model.ny); 
model.pBdryY = zeros(model.nx, 2); 
model.pBdryX(2,end) = initPres(end);
model.pBdryY(end,2) = initPres(end);

model.pInit = zeros(model.nx, model.ny); 
for j = 1 : nx
    model.pInit(j, 1:end) = initPres;
end

model.soludoc = 'case1';

RST_singlePhaseFlow( model )   

rmpath('/Users/wuy/Research/MulticomponentFlow_matlab');

