% This is the input file of two-component and single-phase flow. The injection
% point is at the front-left-down corner and the production point is at the front-left-up
% corner. At the initial time, there is C3H8 in the well. Then we impose pressure
% at the production point and inject CH4 into the well from the injection point. 
% The temperature is constant in the whole process. 

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on November 13rd, 2012

path('/Users/wuy/Research/MulticomponentFlow_matlab', path);

model.numDims = 3;  
model.Nc = 2;
model.phase = 'g'; 

Lx = 4; % unit: meter  
Ly = 4;
Lz = 4;
timeEnd = 0.1*365*24*3600;
nx = 40;  
ny = 40; 
nz = 40;
nt = 0.1*365*24;
model.nx = nx;   
model.ny = ny;   
model.nz = nz; 
model.nt = nt;
model.xs = (0:nx)*Lx/nx; 
model.ys = (0:ny)*Ly/ny; 
model.zs = (0:nz)*Lz/nz;  
model.ts = (0:nt)*timeEnd/nt; 

K_const = 9.869233*10^(-15); % unit: m^2 
model.Kxx = zeros(model.nx, model.ny, model.nz);  
model.Kxx(1:end, 1:end, 1:end) = K_const; 
model.Kyy = zeros(model.nx, model.ny, model.nz);  
model.Kyy(1:end, 1:end, 1:end) = K_const; 
model.Kzz = zeros(model.nx, model.ny, model.nz);  
model.Kzz(1:end, 1:end, 1:end) = K_const; 

model.poro = zeros(model.nx, model.ny, model.nz);  
model.poro(1:end, 1:end, 1:end) = 0.2;

model.T = 480; % unit: K

model.gravX = 0;
model.gravY = -9.807; % m/s^2
model.gravZ = 0;

model.ct = zeros(model.Nc, 1); % critical T
model.cp = zeros(model.Nc, 1); % critical pressure
model.cv = zeros(model.Nc, 1); % critical volume
model.af = zeros(model.Nc, 1); % acentric factor
model.mw = zeros(model.Nc, 1);
model.src = zeros(model.Nc, model.nx, model.ny, model.nz); 

% the first component is methane
model.ct(1) = 190; % unit: K
model.cp(1) = 46*10^5; % unit: Pa
model.cv(1) = 0.0062; % unit: m^3/kg
model.af(1) = 0.01;
model.mw(1) = 16*10^(-3); % unit: kg/mol

% the second component is propane
model.ct(2) = 370; % unit: K
model.cp(2) = 42*10^5; % unit: Pa
model.cv(2) = 0.0045; % unit: m^3/kg
model.af(2) = 0.15;
model.mw(2) = 44*10^(-3); 

% the binary interaction parameters
% reference: Page 155 in A. Firoozabadi's book
model.delta = zeros(model.Nc, model.Nc);
model.delta(1,2) = 0.0289 + 1.633*0.1*model.mw(2);
model.delta(2,1) = model.delta(1,2);

model.isDiriX = zeros(2, model.ny, model.nz);  
model.isDiriY = zeros(model.nx, 2, model.nz); 
model.isDiriZ = zeros(model.nx, model.ny, 2);
model.isDiriX(1, end, 1) = 1; 
model.isDiriY(1, 2, 1) = 1; 
model.isDiriZ(1, end, 1) = 1; 

model.xInit = zeros(model.Nc, model.nx, model.ny, model.nz);
model.xInit(1, 1:end, 1:end, 1:end) = 0.0;
model.xInit(2, 1:end, 1:end, 1:end) = 1.0;

model.xBdryX = zeros(model.Nc, 2, model.ny, model.nz); 
model.xBdryY = zeros(model.Nc, model.nx, 2, model.nz);  
model.xBdryZ = zeros(model.Nc, model.nx, model.ny, 2);
model.xBdryX(1, 1, 1, 1) = 1.0;
model.xBdryY(1, 1, 1, 1) = 1.0;
model.xBdryZ(1, 1, 1, 1) = 1.0;
model.xBdryX(2, 1, 1:end, 1:end) = 1.0;
model.xBdryX(2, 2, 1:end, 1:end) = 1.0;
model.xBdryY(2, 1:end, 1, 1:end) = 1.0;
model.xBdryY(2, 1:end, 2, 1:end) = 1.0;
model.xBdryZ(2, 1:end, 1:end, 1) = 1.0;
model.xBdryZ(2, 1:end, 1:end, 2) = 1.0;
model.xBdryX(2, 1, 1, 1) = 0.0; % supplement method
model.xBdryY(2, 1, 1, 1) = 0.0;
model.xBdryZ(2, 1, 1, 1) = 0.0;

uconst = 2*10^(-6); % unit: m/s
model.uBdryX = zeros(2, model.ny, model.nz); 
model.uBdryY = zeros(model.nx, 2, model.nz);  
model.uBdryZ = zeros(model.nx, model.ny, 2); 
model.uBdryX(1, 1, 1) = -uconst;
model.uBdryY(1, 1, 1) = -uconst;
model.uBdryZ(1, 1, 1) = -uconst;

% load initial pressures
initPres = zeros(1,ny);
temp = load('../../InitialReservoirData/P1.txt');
k = 1;
for i = 1 : ny
    initPres(i) = temp(k);
    k = k + 1;
end

model.pBdryX = zeros(2, model.ny, model.nz); 
model.pBdryY = zeros(model.nx, 2, model.nz); 
model.pBdryZ = zeros(model.nx, model.ny, 2); 
model.pBdryX(1, end, 1) = initPres(end);
model.pBdryY(1, 2, 1) = initPres(end);
model.pBdryZ(1, end, 1) = initPres(end);

model.pInit = zeros(model.nx, model.ny, model.nz); 
for i = 1 : nx
    for k = 1 : nz
        model.pInit(i, 1:end, k) = initPres;
    end
end

model.soludoc = 'case4';

RST_singlePhaseFlow( model ) 

rmpath('/Users/wuy/Research/MulticomponentFlow_matlab');

