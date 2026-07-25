% The function RST_PREOS is to use the PR-EOS to get the fluid molar
% density, the mass density and the derivation of xi to p of the cell

% Input parameters:
% model: the model
% x: the mole fraction of the cell, it is an array of size (Nc,1)
% p: the pressure of the cell, unit: Pa

% Return value:
% xi: the fluid molar density of the cell, unit: mol/m^3 
% rho: the mass density of the cell, unit: kg/m^3
% deri_xi_p: the derivation of xi to p

% Reference: Chap 7. in <<Reservoir Simulation: Mathematical Techniques in Oil Recovery>>

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on October 23rd, 2013      

function [ xi, rho, deri_xi_p ] = RST_PREOS( model, x, p )
    
    % simplify the notations of model
    Nc = model.Nc;
    T = model.T;
    phase = model.phase;
    ct = model.ct;
    cp = model.cp;
    mw = model.mw;
    af = model.af;
    delta = model.delta;
    
    R = 8.314; % gas constant, unit: J/(mol*K)
    
    lambda = zeros(Nc, 1);
    for m = 1 : Nc
        lambda(m) = 0.37464 + 1.5423*af(m) - 0.26992*af(m)^2;
    end
    
    alpha = zeros(Nc, 1);
    for m = 1 : Nc
        alpha(m) = (1+lambda(m)*(1-sqrt(T/ct(m))))^2;
    end
    
    am = zeros(Nc, 1);
    bm = zeros(Nc, 1);
    for m = 1 : Nc
        am(m) = 0.45724*alpha(m)*R^2*ct(m)^2/cp(m);
        bm(m) = 0.077796*R*ct(m)/cp(m);
    end
    
    a = 0;
    for m = 1 : Nc
        for v = 1 : Nc
            a = a + x(m)*x(v)*(1-delta(m,v))*sqrt(am(m)*am(v));
        end
    end
    
    b = 0;
    for m = 1 : Nc
        b = b + x(m)*bm(m);
    end
    
    A = a*p/(R*T)^2;
    B = b*p/(R*T);
    
    % compressibility factor
    Z = roots([1 -(1-B) (A-3*B^2-2*B) -(A*B-B^2-B^3)]);
    
    if(phase == 'l')
        ZRP = [];
        for i = 1 : 3
            if(isreal(Z(i)) && (Z(i)>0))
                ZRP = [ZRP Z(i)];   
            end
        end
        Z = min(ZRP);
        
        % the fluid molar density, before shifting
        xi_bs = p/(R*T*Z);
        
        % shift the liquid volume and correct the molar density in liquid phase
        % reference: <<GENERALIZED LIQUID VOLUME SHIFTS FOR THE PENGROBINSON
        % EQUATION OF STATE FOR C1 TO C8 HYDROCARBONS>>
        c = zeros(Nc,1); % volume correct item for each component
        C1toC2 = zeros(Nc,1);
        for m = 1 : Nc
            C1toC2(m) = 110.07*af(m)^4 - 83.807*af(m)^3 + 18.926*af(m)^2 - 1.6348*af(m) - 0.0066;
        end
        C2 = 2.013645*10^(-3); % unit: m^3/kg
        C3 = 0.89;
        for m = 1 : Nc
            c(m) = C1toC2(m)*C2 + C2*(T/ct(m)-C3)^2;
        end
        ctotal = 0;
        for m = 1 : Nc 
            ctotal = ctotal + x(m)*c(m)*mw(m); % unit: m^3/mol
        end
        xi = 1/(1/xi_bs + ctotal); % unit: m^3/mol
    
        % the mass density
        rho = 0;
        for m = 1 : Nc
            rho = rho + x(m)*mw(m);
        end
        rho = rho*xi;
    
        deri_A_p = a/(R*T)^2;
        deri_B_p = b/(R*T);
        deri_Z_p = -(deri_B_p*Z^2+(deri_A_p-2*(1+3*B)*deri_B_p)*Z-(deri_A_p*B+(A-2*B-3*B^2)*deri_B_p))/(3*Z^2-2*(1-B)*Z+(A-2*B-3*B^2));
        deri_xi_p_bs = 1/(R*T*Z) - p/(R*T*Z^2)*deri_Z_p;
        deri_xi_p = (1/(1+ctotal*xi_bs)^2)*deri_xi_p_bs;
   
    elseif(phase == 'g')
        ZR = [];
        for i = 1 : 3
            if isreal(Z(i))
                ZR = [ZR Z(i)];   
            end
        end   
        Z = max(ZR); 

        % the fluid molar density
        xi = p/(R*T*Z);
    
        % the mass density
        rho = 0;
        for m = 1 : Nc
            rho = rho + x(m)*mw(m);
        end
        rho = rho*xi;
    
        deri_A_p = a/(R*T)^2;
        deri_B_p = b/(R*T);
        deri_Z_p = -(deri_B_p*Z^2+(deri_A_p-2*(1+3*B)*deri_B_p)*Z-(deri_A_p*B+(A-2*B-3*B^2)*deri_B_p))/(3*Z^2-2*(1-B)*Z+(A-2*B-3*B^2));
        deri_xi_p = 1/(R*T*Z) - p/(R*T*Z^2)*deri_Z_p;
    end
    
end
