% This is the RST_singlePhaseFlow() function which can be called by the
% input files. It computes the pressures of the cells and then give out 
% the velocities. And finally, it tells you the mole fraction of each
% component. After that, it outputs the results to a series of solution
% files and then calls the RST_plot() function to draw the images of the 
% results.

% Input parameters:
% model: the model

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on October 23rd, 2012

function RST_singlePhaseFlow( model )

    % simplify the symbols of the model
    Nc = model.Nc;
    nx = model.nx;
    ny = model.ny;
    nz = model.nz;
    nt = model.nt;
    xs = model.xs;
    ys = model.ys;
    zs = model.zs;
    ts = model.ts;
    poro = model.poro;
    Kxx = model.Kxx;
    Kyy = model.Kyy;   
    Kzz = model.Kzz; 
    isDiriX = model.isDiriX;
    isDiriY = model.isDiriY;
    isDiriZ = model.isDiriZ;
    pBdryX = model.pBdryX;
    pBdryY = model.pBdryY;
    pBdryZ = model.pBdryZ;
    uBdryX = model.uBdryX;
    uBdryY = model.uBdryY;
    uBdryZ = model.uBdryZ;
    xBdryX = model.xBdryX;
    xBdryY = model.xBdryY; 
    xBdryZ = model.xBdryZ; 
    pInit = model.pInit;
    xInit = model.xInit;
    gravX = model.gravX;
    gravY = model.gravY;
    gravZ = model.gravZ;
    src = model.src;
    soludoc = model.soludoc;
    
    % create the new document to store the results
    if(exist(soludoc, 'dir') == 0)
         mkdir(soludoc);
    end
    
    % the file to store the mole fraction history of component 2
    fmhtxt = [soludoc, '/soln_1PhFlw_moleFractionHistory.txt'];
    fmhtxtid = fopen(fmhtxt, 'w');
    
    % the file to store the left mole ratio history
    fmrtxt = [soludoc, '/soln_1PhFlw_moleRatioHistory.txt'];
    fmrtxtid = fopen(fmrtxt, 'w');
  
    % rectify the directions of the velocities
    uBdryX(1, 1:end, 1:end) = -uBdryX(1, 1:end, 1:end);
    uBdryY(1:end, 1, 1:end) = -uBdryY(1:end, 1, 1:end);
    uBdryZ(1:end, 1:end, 1) = -uBdryZ(1:end, 1:end, 1);
    
    % define the pressures and initialize them
    p = zeros(nx+2, ny+2, nz+2);
    p(1, 2:ny+1, 2:nz+1) = pBdryX(1, 1:end, 1:end);
    p(nx+2, 2:ny+1, 2:nz+1) = pBdryX(2, 1:end, 1:end);
    p(2:nx+1, 1, 2:nz+1) = pBdryY(1:end, 1, 1:end);
    p(2:nx+1, ny+2, 2:nz+1) = pBdryY(1:end, 2, 1:end);
    p(2:nx+1, 2:ny+1, 1) = pBdryZ(1:end, 1:end, 1);
    p(2:nx+1, 2:ny+1, nz+2) = pBdryZ(1:end, 1:end, 2);
    
    for k = 2 : nz+1
        for j = 2 : ny+1
            for i = 2 : nx+1
                p(i,j,k) = pInit(i-1, j-1, k-1);
            end
        end
    end
    
    % define the fluid velocities and initialize them
    ux = zeros(nx+1, ny, nz);
    uy = zeros(nx, ny+1, nz);
    uz = zeros(nx, ny, nz+1);
    ux(1, 1:end, 1:end) = uBdryX(1, 1:end, 1:end);
    ux(nx+1, 1:end, 1:end) = uBdryX(2, 1:end, 1:end);
    uy(1:end, 1, 1:end) = uBdryY(1:end, 1, 1:end);
    uy(1:end, ny+1, 1:end) = uBdryY(1:end, 2, 1:end);
    uz(1:end, 1:end, 1) = uBdryZ(1:end, 1:end, 1);
    uz(1:end, 1:end, nz+1) = uBdryZ(1:end, 1:end, 2);
        
    % define the mole fraction of each component in the cell and initialize them
    x = zeros(Nc, nx, ny, nz);
    for m = 1 : Nc
        x(m, 1:end, 1:end, 1:end) = xInit(m, 1:end, 1:end, 1:end);
    end
    
    % define the mole fraction of each component on the faces and initialize them
    xfacex = zeros(Nc, nx+1, ny, nz);
    xfacey = zeros(Nc, nx, ny+1, nz);
    xfacez = zeros(Nc, nx, ny, nz+1);
    
    % define the mobility*(fluid molar density) on the faces
    lambdax = zeros(nx+1, ny, nz);
    lambday = zeros(nx, ny+1, nz);
    lambdaz = zeros(nx, ny, nz+1);
    
    % define the fluid density
    xi = zeros(nx, ny, nz);
    
    % define the mole density on the bars
    xifacex = zeros(nx+1, ny, nz);
    xifacey = zeros(nx, ny+1, nz);
    xifacez = zeros(nx, ny, nz+1);
    
    % define the mass density in the cell
    rho = zeros(nx, ny, nz);
    
    % define the mass density on the faces
    rhofacex = zeros(nx+1, ny, nz);
    rhofacey = zeros(nx, ny+1, nz);
    rhofacez = zeros(nx, ny, nz+1);
    
    % define the derivation of fluid density to pressure 
    deri_xi_p = zeros(nx, ny, nz);
    
    % define the viscosity
    mu = zeros(nx, ny, nz);    
    
    % define the permeabilities on the faces
    Kxxface = zeros(nx+1, ny, nz);
    Kyyface = zeros(nx, ny+1, nz);
    Kzzface = zeros(nx, ny, nz+1);
    
    % initialize Kxxface using harmonic weighting method
    for k = 1 : nz
        for j = 1 : ny
            for i = 1 : nx+1
                if(i == 1)
                    Kxxface(i,j,k) = Kxx(i,j,k);
                elseif(i == nx+1)
                    Kxxface(i,j,k) = Kxx(i-1, j, k);
                else
                    ltotal = xs(i+1) - xs(i-1);
                    lleft = xs(i) - xs(i-1);
                    lright = xs(i+1) - xs(i);
                    Kxxface(i,j,k) = ltotal / (lleft/Kxx(i-1, j, k)+lright/Kxx(i,j,k));
                end
            end
        end
    end
    
    % initialize Kyyface using harmonic weighting method
    for k = 1 : nz
        for j = 1 : ny+1  
            for i = 1 : nx
                if(j == 1)
                    Kyyface(i,j,k) = Kyy(i,j,k);
                elseif(j == ny+1)
                    Kyyface(i,j,k) = Kyy(i, j-1, k);
                else
                    ltotal = ys(j+1) - ys(j-1);
                    lup = ys(j+1) - ys(j);
                    ldown = ys(j) - ys(j-1);
                    Kyyface(i,j,k) = ltotal / (ldown/Kyy(i, j-1, k)+lup/Kyy(i,j,k));
                end
            end
        end
    end
    
    % initialize Kzzface using harmonic weighting method
    for k = 1 : nz+1
        for j = 1 : ny  
            for i = 1 : nx
                if(k == 1)
                    Kzzface(i,j,k) = Kzz(i,j,k);
                elseif(k == nz+1)
                    Kzzface(i,j,k) = Kzz(i, j, k-1);
                else
                    ltotal = zs(k+1) - zs(k-1);
                    lback = zs(k+1) - zs(k);
                    lfront = zs(k) - zs(k-1);
                    Kzzface(i,j,k) = ltotal / (lfront/Kzz(i,j, k-1)+lback/Kzz(i,j,k));
                end
            end
        end
    end
    
    firsttime = 1;
    totalmole = 0.0;
    % define the coefficient matrices A and b
    A = sparse(nx*ny*nz, nx*ny*nz);
    b = zeros(nx*ny*nz, 1);
    moleincell = zeros(Nc,nx,ny,nz);
    leftmole = zeros(Nc,1);
    % define the new molar density of each component
    c = zeros(Nc, nx, ny, nz);
    
    % time iteration 
    for t = 2 : nt+1
        
        t    
        
        % compute the PR EOS and viscosity
        for k = 1 : nz
            for j = 1 : ny
                for i = 1 : nx
                    [ xi(i,j,k), rho(i,j,k), deri_xi_p(i,j,k) ] = RST_PREOS( model, x(1:end, i, j, k), p(i+1, j+1, k+1) ); 
                    [ mu(i,j,k) ] = RST_viscosity( model, x(1:end, i, j, k), xi(i,j,k), p(i+1, j+1, k+1) ); 
                end
            end
        end
        
        for k = 1 : nz
            for j = 1 : ny
                for i = 1 : nx
                    for m = 1 : Nc
                        moleincell(m,i,j,k) = xi(i,j,k)*x(m,i,j,k)*(xs(i+1)-xs(i))*(ys(j+1)-ys(j))*(zs(k+1)-zs(k))*poro(i,j,k);
                    end
                end 
            end 
        end 

        leftmole(:,1) = 0;
        for k = 1 : nz
            for j = 1 : ny
                for i = 1 : nx
                    for m = 1 : Nc
                        leftmole(m) = leftmole(m) + moleincell(m,i,j,k);
                    end
                end 
            end 
        end 

        totaldesiredleftmole = 0.0;
        for m = 2 : Nc
            totaldesiredleftmole = totaldesiredleftmole + leftmole(m);
        end 

        if(firsttime == 1)
            firsttime = 0;
            totalmole = totaldesiredleftmole;
        end
        
        totaldesiredleftmole = 0.0;
        for m = 2 : Nc
            totaldesiredleftmole = totaldesiredleftmole + leftmole(m);
        end
        
        fprintf(fmhtxtid, '%10e\n', (totalmole-totaldesiredleftmole)/totalmole);
        
        % print the mole ratio of desired components to component 1 in
        % the well
        fprintf(fmrtxtid, '%10e\n', totaldesiredleftmole/leftmole(1));

		% compute the time step
		timestep = ts(2) - ts(1);
        
        % compute lambda = mobility*(fluid molar density) on the faces,
        % using single-point upstream weighting
        for k = 1 : nz
            for j = 1 : ny
                if(ux(1,j,k) > 0)
                    [ xifacex(1,j,k), rhofacex(1,j,k), ~ ] = RST_PREOS( model, xBdryX(1:Nc,1,j,k), p(2,j+1,k+1) ); 
                    [ mutemp ] = RST_viscosity( model, xBdryX(1:Nc,1,j,k), xifacex(1,j,k), p(2,j+1,k+1) );
                    lambdax(1,j,k) = Kxxface(1,j,k)*xifacex(1,j,k)/mutemp;
                    xfacex(1:Nc,1,j,k) = xBdryX(1:Nc,1,j,k);
                elseif(ux(1,j,k) < 0)
                    rhofacex(1,j,k) = rho(1,j,k);
                    xifacex(1,j,k) = xi(1,j,k);
                    lambdax(1,j,k) = Kxxface(1,j,k)*xifacex(1,j,k)/mu(1,j,k);  
                    xfacex(1:Nc,1,j,k) = x(1:Nc,1,j,k);
                else
                    if(isDiriX(1,j,k) == 1)
                        rhofacex(1,j,k) = rho(1,j,k);
                        xifacex(1,j,k) = xi(1,j,k);
                        lambdax(1,j,k) = Kxxface(1,j,k)*xifacex(1,j,k)/mu(1,j,k);  
                        xfacex(1:Nc,1,j,k) = x(1:Nc,1,j,k);
                    end
                end
            end
        end
        
        for k = 1 : nz
            for j = 1 : ny
                if(ux(nx+1,j,k) < 0)
                    [ xifacex(nx+1,j,k), rhofacex(nx+1,j,k), ~ ] = RST_PREOS( model, xBdryX(1:Nc,2,j,k), p(nx+1,j+1,k+1) ); 
                    [ mutemp ] = RST_viscosity( model, xBdryX(1:Nc,2,j,k), xifacex(nx+1,j,k), p(nx+1,j+1,k+1) );
                    lambdax(nx+1,j,k) = Kxxface(nx+1,j,k)*xifacex(nx+1,j,k)/mutemp;
                    xfacex(1:Nc,nx+1,j,k) = xBdryX(1:Nc,2,j,k);
                elseif(ux(nx+1,j,k) > 0)
                    rhofacex(nx+1,j,k) = rho(nx,j,k);
                    xifacex(nx+1,j,k) = xi(nx,j,k);
                    lambdax(nx+1,j,k) = Kxxface(nx+1,j,k)*xifacex(nx+1,j,k)/mu(nx,j,k);  
                    xfacex(1:Nc,nx+1,j,k) = x(1:Nc,nx,j,k);
                else
                    if(isDiriX(2,j,k) == 1)
                        rhofacex(nx+1,j,k) = rho(nx,j,k);
                        xifacex(nx+1,j,k) = xi(nx,j,k);
                        lambdax(nx+1,j,k) = Kxxface(nx+1,j,k)*xifacex(nx+1,j,k)/mu(nx,j,k);  
                        xfacex(1:Nc,nx+1,j,k) = x(1:Nc,nx,j,k);
                    end
                end
            end
        end
        
        for k = 1 : nz
            for j = 1 : ny
                for i = 2 : nx
                    if(ux(i,j,k) > 0)
                        rhofacex(i,j,k) = rho(i-1,j,k);
                        xifacex(i,j,k) = xi(i-1,j,k);
                        lambdax(i,j,k) = Kxxface(i,j,k)*xifacex(i,j,k)/mu(i-1,j,k);  
                        xfacex(1:Nc,i,j,k) = x(1:Nc,i-1,j,k);
                    else
                        rhofacex(i,j,k) = rho(i,j,k);
                        xifacex(i,j,k) = xi(i,j,k);
                        lambdax(i,j,k) = Kxxface(i,j,k)*xifacex(i,j,k)/mu(i,j,k);  
                        xfacex(1:Nc,i,j,k) = x(1:Nc,i,j,k);
                    end
                end
            end
        end
        
        for k = 1 : nz
            for i = 1 : nx
                if(uy(i,1,k) > 0)
                    [ xifacey(i,1,k), rhofacey(i,1,k), ~ ] = RST_PREOS( model, xBdryY(1:Nc,i,1,k), p(i+1,2,k+1) ); 
                    [ mutemp ] = RST_viscosity( model, xBdryY(1:Nc,i,1,k), xifacey(i,1,k), p(i+1,2,k+1) );
                    lambday(i,1,k) = Kyyface(i,1,k)*xifacey(i,1,k)/mutemp;
                    xfacey(1:Nc,i,1,k) = xBdryY(1:Nc,i,1,k);
                elseif(uy(i,1,k) < 0)
                    rhofacey(i,1,k) = rho(i,1,k);
                    xifacey(i,1,k) = xi(i,1,k);
                    lambday(i,1,k) = Kyyface(i,1,k)*xifacey(i,1,k)/mu(i,1,k);  
                    xfacey(1:Nc,i,1,k) = x(1:Nc,i,1,k);
                else
                    if(isDiriY(i,1,k) == 1)
                        rhofacey(i,1,k) = rho(i,1,k);
                        xifacey(i,1,k) = xi(i,1,k);
                        lambday(i,1,k) = Kyyface(i,1,k)*xifacey(i,1,k)/mu(i,1,k);  
                        xfacey(1:Nc,i,1,k) = x(1:Nc,i,1,k);
                    end
                end
            end
        end
        
        for k = 1 : nz
            for i = 1 : nx
                if(uy(i,ny+1,k) < 0)
                    [ xifacey(i,ny+1,k), rhofacey(i,ny+1,k), ~ ] = RST_PREOS( model, xBdryY(1:Nc,i,2,k), p(i+1,ny+1,k+1) ); 
                    [ mutemp ] = RST_viscosity( model, xBdryY(1:Nc,i,2,k), xifacey(i,ny+1,k), p(i+1,ny+1,k+1) );
                    lambday(i,ny+1,k) = Kyyface(i,ny+1,k)*xifacey(i,ny+1,k)/mutemp;
                    xfacey(1:Nc,i,ny+1,k) = xBdryY(1:Nc,i,2,k);
                elseif(uy(i,ny+1,k) > 0)
                    rhofacey(i,ny+1,k) = rho(i,ny,k);
                    xifacey(i,ny+1,k) = xi(i,ny,k);
                    lambday(i,ny+1,k) = Kyyface(i,ny+1,k)*xifacey(i,ny+1,k)/mu(i,ny,k);  
                    xfacey(1:Nc,i,ny+1,k) = x(1:Nc,i,ny,k);
                else
                    if(isDiriY(i,2,k) == 1)
                        rhofacey(i,ny+1,k) = rho(i,ny,k);
                        xifacey(i,ny+1,k) = xi(i,ny,k);
                        lambday(i,ny+1,k) = Kyyface(i,ny+1,k)*xifacey(i,ny+1,k)/mu(i,ny,k);  
                        xfacey(1:Nc,i,ny+1,k) = x(1:Nc,i,ny,k);
                    end
                end
            end
        end
        
        for k = 1 : nz
            for j = 2 : ny
                for i = 1 : nx
                    if(uy(i,j,k) > 0)
                        rhofacey(i,j,k) = rho(i,j-1,k);
                        xifacey(i,j,k) = xi(i,j-1,k);
                        lambday(i,j,k) = Kyyface(i,j,k)*xifacey(i,j,k)/mu(i,j-1,k);  
                        xfacey(1:Nc,i,j,k) = x(1:Nc,i,j-1,k);
                    else
                        rhofacey(i,j,k) = rho(i,j,k);
                        xifacey(i,j,k) = xi(i,j,k);
                        lambday(i,j,k) = Kyyface(i,j,k)*xifacey(i,j,k)/mu(i,j,k);  
                        xfacey(1:Nc,i,j,k) = x(1:Nc,i,j,k);
                    end
                end
            end
        end
        
        for j = 1 : ny
            for i = 1 : nx
                if(uz(i,j,1) > 0)
                    [ xifacez(i,j,1), rhofacez(i,j,1), ~ ] = RST_PREOS( model, xBdryZ(1:Nc,i,j,1), p(i+1,j+1,2) ); 
                    [ mutemp ] = RST_viscosity( model, xBdryZ(1:Nc,i,j,1), xifacez(i,j,1), p(i+1,j+1,2) );
                    lambdaz(i,j,1) = Kzzface(i,j,1)*xifacez(i,j,1)/mutemp;
                    xfacez(1:Nc,i,j,1) = xBdryZ(1:Nc,i,j,1);
                elseif(uz(i,j,1) < 0)
                    rhofacez(i,j,1) = rho(i,j,1);
                    xifacez(i,j,1) = xi(i,j,1);
                    lambdaz(i,j,1) = Kzzface(i,j,1)*xifacez(i,j,1)/mu(i,j,1);  
                    xfacez(1:Nc,i,j,1) = x(1:Nc,i,j,1);
                else
                    if(isDiriZ(i,j,1) == 1)
                        rhofacez(i,j,1) = rho(i,j,1);
                        xifacez(i,j,1) = xi(i,j,1);
                        lambdaz(i,j,1) = Kzzface(i,j,1)*xifacez(i,j,1)/mu(i,j,1);  
                        xfacez(1:Nc,i,j,1) = x(1:Nc,i,j,1);
                    end
                end
            end
        end
        
        for j = 1 : ny
            for i = 1 : nx
                if(uz(i,j,nz+1) < 0)
                    [ xifacez(i,j,nz+1), rhofacez(i,j,nz+1), ~ ] = RST_PREOS( model, xBdryZ(1:Nc,i,j,2), p(i+1,j+1,nz+1) ); 
                    [ mutemp ] = RST_viscosity( model, xBdryZ(1:Nc,i,j,2), xifacez(i,j,nz+1), p(i+1,j+1,nz+1) );
                    lambdaz(i,j,nz+1) = Kzzface(i,j,nz+1)*xifacez(i,j,nz+1)/mutemp;
                    xfacez(1:Nc,i,j,nz+1) = xBdryZ(1:Nc,i,j,2);
                elseif(uz(i,j,nz+1) > 0)
                    rhofacez(i,j,nz+1) = rho(i,j,nz);
                    xifacez(i,j,nz+1) = xi(i,j,nz);
                    lambdaz(i,j,nz+1) = Kzzface(i,j,nz+1)*xifacez(i,j,nz+1)/mu(i,j,nz);  
                    xfacez(1:Nc,i,j,nz+1) = x(1:Nc,i,j,nz);
                else
                    if(isDiriZ(i,j,2) == 1)
                        rhofacez(i,j,nz+1) = rho(i,j,nz);
                        xifacez(i,j,nz+1) = xi(i,j,nz);
                        lambdaz(i,j,nz+1) = Kzzface(i,j,nz+1)*xifacez(i,j,nz+1)/mu(i,j,nz);  
                        xfacez(1:Nc,i,j,nz+1) = x(1:Nc,i,j,nz);
                    end
                end
            end
        end
        
        for k = 2 : nz
            for j = 1 : ny
                for i = 1 : nx
                    if(uz(i,j,k) > 0)
                        rhofacez(i,j,k) = rho(i,j,k-1);
                        xifacez(i,j,k) = xi(i,j,k-1);
                        lambdaz(i,j,k) = Kzzface(i,j,k)*xifacez(i,j,k)/mu(i,j,k-1);  
                        xfacez(1:Nc,i,j,k) = x(1:Nc,i,j,k-1);
                    else
                        rhofacez(i,j,k) = rho(i,j,k);
                        xifacez(i,j,k) = xi(i,j,k);
                        lambdaz(i,j,k) = Kzzface(i,j,k)*xifacez(i,j,k)/mu(i,j,k);  
                        xfacez(1:Nc,i,j,k) = x(1:Nc,i,j,k);
                    end
                end
            end
        end
        
        % compute the coefficiencies of the pressure equations   
        for k = 2 : nz+1
            for j = 2 : ny+1
                for i = 2 : nx+1     

                    r = (k-2)*nx*ny + (j-2)*nx + (i-1); % r is the index of the current cell
                
                    xedge = xs(i) - xs(i-1); % x-edge of the current cell
                    yedge = ys(j) - ys(j-1); % y-edge of the current cell
                    zedge = zs(k) - zs(k-1); % z-edge of the current cell
            
                    % ledge means the x-edge of the left cell 
                    if(i ~= 2)
                        ledge = xs(i-1) - xs(i-2);
                    else
                        ledge = 0;
                    end
            
                    % redge means the x-edge of the right cell
                    if(i ~= nx+1)
                        redge = xs(i+1) - xs(i);    
                    else
                        redge = 0;
                    end
            
                    % uedge means the y-edge of the up cell
                    if(j ~= ny+1)
                        uedge = ys(j+1) - ys(j);
                    else
                        uedge = 0;
                    end
            
                    % dedge means the y-edge of the down cell
                    if(j ~= 2)
                        dedge = ys(j-1) - ys(j-2);
                    else
                        dedge = 0;
                    end
                    
                    % bedge means the z-edge of the back cell
                    if(k ~= nz+1)
                        bedge = zs(k+1) - zs(k);
                    else
                        bedge = 0;
                    end
            
                    % fedge means the z-edge of the front cell
                    if(k ~= 2)
                        fedge = zs(k-1) - zs(k-2);
                    else
                        fedge = 0;
                    end

                    % compute the coefficiencies of the pressure equations
                    left = -2*lambdax(i-1, j-1, k-1)/xedge/(xedge+ledge); % left means the coefficiency of the pressure of the left cell
                    right = -2*lambdax(i, j-1, k-1)/xedge/(xedge+redge); % right means the coefficiency of the pressure of the right cell
                    up = -2*lambday(i-1, j, k-1)/yedge/(yedge+uedge); % up means the coefficiency of the pressure of the up cell
                    down = -2*lambday(i-1, j-1, k-1)/yedge/(yedge+dedge); % down means the coefficiency of the pressure of the down cell   
                    front = -2*lambdaz(i-1, j-1, k-1)/zedge/(zedge+fedge); % front means the coefficiency of the pressure of the front cell
                    back = -2*lambdaz(i-1, j-1, k)/zedge/(zedge+bedge); % back means the coefficiency of the pressure of the back cell

                    cp = poro(i-1, j-1, k-1) * deri_xi_p(i-1, j-1, k-1); % suppose porosity was constant here
                
                    q = 0;
                    for m = 1 : Nc
                        q = q + src(m, i-1, j-1, k-1);
                    end
                
                    A(r,r) = cp/timestep;
                    b(r,1) = q + cp*p(i,j,k)/timestep;         
                
                    if((i == 2) && (isDiriX(1, j-1, k-1) == 0))
                        b(r,1) = b(r,1) + uBdryX(1, j-1, k-1)*xifacex(i-1, j-1, k-1)/xedge;
                    elseif((i == 2) && (isDiriX(1, j-1, k-1) == 1))
                        b(r,1) = b(r,1) - left*p(i-1, j, k) + lambdax(i-1, j-1, k-1)*rhofacex(i-1, j-1, k-1)*gravX/xedge;
                        A(r,r) = A(r,r) - left;
                    else
                        A(r, r-1) = left;
                        A(r,r) = A(r,r) - left;
                        b(r,1) = b(r,1) + lambdax(i-1, j-1, k-1)*rhofacex(i-1, j-1, k-1)*gravX/xedge;
                    end
                
                    if((i == nx+1) && (isDiriX(2, j-1, k-1) == 0))
                        b(r,1) = b(r,1) - uBdryX(2, j-1, k-1)*xifacex(i, j-1, k-1)/xedge;  
                    elseif((i == nx+1) && (isDiriX(2, j-1, k-1) == 1))
                        b(r,1) = b(r,1) - right*p(i+1, j, k) - lambdax(i, j-1, k-1)*rhofacex(i, j-1, k-1)*gravX/xedge;
                        A(r,r) = A(r,r) - right;
                    else
                        A(r, r+1) = right;
                        A(r,r) = A(r,r) - right;
                        b(r,1) = b(r,1) - lambdax(i, j-1, k-1)*rhofacex(i, j-1, k-1)*gravX/xedge;
                    end
                    
                    if((j == ny+1) && (isDiriY(i-1, 2, k-1) == 0)) % Neumann boundary
                        b(r,1) = b(r,1) - uBdryY(i-1, 2, k-1)*xifacey(i-1, j, k-1)/yedge;
                    elseif((j == ny+1) && (isDiriY(i-1, 2, k-1) == 1)) % Dirichlet boundary
                        b(r,1) = b(r,1) - up*p(i, j+1, k) - lambday(i-1, j, k-1)*rhofacey(i-1, j, k-1)*gravY/yedge;
                        A(r,r) = A(r,r) - up;
                    else
                        A(r, r+nx) = up;
                        A(r,r) = A(r,r) - up;
                        b(r,1) = b(r,1) - lambday(i-1, j, k-1)*rhofacey(i-1, j, k-1)*gravY/yedge;
                    end
                
                    if((j == 2) && (isDiriY(i-1, 1, k-1) == 0))
                        b(r,1) = b(r,1) + uBdryY(i-1, 1, k-1)*xifacey(i-1, j-1, k-1)/yedge;
                    elseif((j == 2) && (isDiriY(i-1, 1, k-1) == 1))
                        b(r,1) = b(r,1) - down*p(i, j-1, k) + lambday(i-1, j-1, k-1)*rhofacey(i-1, j-1, k-1)*gravY/yedge;
                        A(r,r) = A(r,r) - down;
                    else
                        A(r, r-nx) = down;
                        A(r,r) = A(r,r) - down;
                        b(r,1) = b(r,1) + lambday(i-1, j-1, k-1)*rhofacey(i-1, j-1, k-1)*gravY/yedge;
                    end
                    
                    if((k == nz+1) && (isDiriZ(i-1, j-1, 2) == 0)) % Neumann boundary
                        b(r,1) = b(r,1) - uBdryZ(i-1, j-1, 2)*xifacez(i-1, j-1, k)/zedge;
                    elseif((k == nz+1) && (isDiriZ(i-1, j-1, 2) == 1)) % Dirichlet boundary
                        b(r,1) = b(r,1) - back*p(i, j, k+1) - lambdaz(i-1, j-1, k)*rhofacez(i-1, j-1, k)*gravZ/zedge;
                        A(r,r) = A(r,r) - back;
                    else
                        A(r, r+nx*ny) = back;
                        A(r,r) = A(r,r) - back;
                        b(r,1) = b(r,1) - lambdaz(i-1, j-1, k)*rhofacez(i-1, j-1, k)*gravZ/zedge;
                    end
                
                    if((k == 2) && (isDiriZ(i-1, j-1, 1) == 0))
                        b(r,1) = b(r,1) + uBdryZ(i-1, j-1, 1)*xifacez(i-1, j-1, k-1)/zedge;
                    elseif((k == 2) && (isDiriZ(i-1, j-1, 1) == 1))
                        b(r,1) = b(r,1) - front*p(i, j, k-1) + lambdaz(i-1, j-1, k-1)*rhofacez(i-1, j-1, k-1)*gravZ/zedge;
                        A(r,r) = A(r,r) - front;
                    else
                        A(r, r-nx*ny) = front;
                        A(r,r) = A(r,r) - front;
                        b(r,1) = b(r,1) + lambdaz(i-1, j-1, k-1)*rhofacez(i-1, j-1, k-1)*gravZ/zedge;
                    end
                    
                end
            end
        end
    
        % compute the new pressures
        xx = A\b;
        for k = 2 : nz+1
            for j = 2 : ny+1
                for i = 2 : nx+1
                    p(i,j,k) = xx((k-2)*nx*ny + (j-2)*nx + (i-1), 1);
                end
            end
        end
        
        % compute the new velocities in x direction
        for k = 1 : nz
            for j = 1 : ny
                for i = 1 : nx+1
                    if((i == 1) && (isDiriX(1,j,k) == 1))
                        ux(i,j,k) = -lambdax(i,j,k)/xifacex(i,j,k) * ((p(i+1, j+1, k+1)-p(i, j+1, k+1))*2/(xs(i+1)-xs(i)) - rhofacex(i,j,k)*gravX);
                    elseif((i == 1) && (isDiriX(1,j,k) == 0))
                        ux(i,j,k) = uBdryX(1,j,k);
                    elseif((i == nx+1) && (isDiriX(2,j,k) == 1))
                        ux(i,j,k) = -lambdax(i,j,k)/xifacex(i,j,k) * ((p(i+1, j+1, k+1)-p(i, j+1, k+1))*2/(xs(i)-xs(i-1)) - rhofacex(i,j,k)*gravX);
                    elseif((i == nx+1) && (isDiriX(2,j,k) == 0))
                        ux(i,j,k) = uBdryX(2,j,k);
                    else
                        ux(i,j,k) = -lambdax(i,j,k)/xifacex(i,j,k) * ((p(i+1, j+1, k+1)-p(i, j+1, k+1))*2/(xs(i+1)-xs(i-1)) - rhofacex(i,j,k)*gravX);
                    end
                end
            end
        end
    
        % compute the new velocities in y direction
        for k = 1 : nz
            for j = 1 : ny+1
                for i = 1 : nx
                    if((j == 1) && (isDiriY(i,1,k) == 1))
                        uy(i,j,k) = -lambday(i,j,k)/xifacey(i,j,k) * ((p(i+1, j+1, k+1)-p(i+1, j, k+1))*2/(ys(j+1)-ys(j)) - rhofacey(i,j,k)*gravY);
                    elseif((j == 1) && (isDiriY(i,1,k) == 0))
                        uy(i,j,k) = uBdryY(i,1,k);
                    elseif((j == ny+1) && (isDiriY(i,2,k) == 1))
                        uy(i,j,k) = -lambday(i,j,k)/xifacey(i,j,k) * ((p(i+1, j+1, k+1)-p(i+1, j, k+1))*2/(ys(j)-ys(j-1)) - rhofacey(i,j,k)*gravY);
                    elseif((j == ny+1) && (isDiriY(i,2,k) == 0))
                        uy(i,j,k) = uBdryY(i,2,k);
                    else
                        uy(i,j,k) = -lambday(i,j,k)/xifacey(i,j,k) * ((p(i+1, j+1, k+1)-p(i+1, j, k+1))*2/(ys(j+1)-ys(j-1)) - rhofacey(i,j,k)*gravY);
                    end
                end
            end
        end
        
        % compute the new velocities in z direction
        for k = 1 : nz+1
            for j = 1 : ny
                for i = 1 : nx
                    if((k == 1) && (isDiriZ(i,j,1) == 1))
                        uz(i,j,k) = -lambdaz(i,j,k)/xifacez(i,j,k) * ((p(i+1, j+1, k+1)-p(i+1, j+1, k))*2/(zs(k+1)-zs(k)) - rhofacez(i,j,k)*gravZ);
                    elseif((k == 1) && (isDiriZ(i,j,1) == 0))
                        uz(i,j,k) = uBdryZ(i,j,1);
                    elseif((k == nz+1) && (isDiriZ(i,j,2) == 1))
                        uz(i,j,k) = -lambdaz(i,j,k)/xifacez(i,j,k) * ((p(i+1, j+1, k+1)-p(i+1, j+1, k))*2/(zs(k)-zs(k-1)) - rhofacez(i,j,k)*gravZ);
                    elseif((k == nz+1) && (isDiriZ(i,j,2) == 0))
                        uz(i,j,k) = uBdryZ(i,j,2);
                    else
                        uz(i,j,k) = -lambdaz(i,j,k)/xifacez(i,j,k) * ((p(i+1, j+1, k+1)-p(i+1, j+1, k))*2/(zs(k+1)-zs(k-1)) - rhofacez(i,j,k)*gravZ);
                    end
                end
            end
        end
        
        % compute the new mole fraction of each component
        for k = 1 : nz
            for j = 1 : ny
                for i = 1 : nx
                    for m = 1 : Nc
                        div = (xfacex(m, i+1, j, k)*ux(i+1, j, k)*xifacex(i+1, j, k) - xfacex(m,i,j,k)*ux(i,j,k)*xifacex(i, j, k))/(xs(i+1)-xs(i)) + ...
                            (xfacey(m, i, j+1, k)*uy(i, j+1, k)*xifacey(i, j+1, k) - xfacey(m,i,j,k)*uy(i,j,k)*xifacey(i, j, k))/(ys(j+1)-ys(j)) + ...
                            (xfacez(m, i, j, k+1)*uz(i, j, k+1)*xifacez(i, j, k+1) - xfacez(m,i,j,k)*uz(i,j,k)*xifacez(i, j, k))/(zs(k+1)-zs(k));
                        c(m,i,j,k) = (src(m,i,j,k)-div)*timestep/poro(i,j,k) + x(m,i,j,k)*xi(i,j,k);
                        if(c(m,i,j,k)<0)
                            c(m,i,j,k);
                        end
                    end
                    ctotal = 0;
                    for m = 1 : Nc 
                        ctotal = ctotal + c(m,i,j,k);
                    end
                    for m = 1 : Nc 
                        x(m,i,j,k) = c(m,i,j,k)/ctotal;
                    end
                end
            end
        end

    end
    
    % close the mole history file
    fclose(fmhtxtid);
    
    % close the mole rate history file
    fclose(fmrtxtid);
    
    % write the results to the solution files
    [ fptxt, fuxtxt, fuytxt, fuztxt, fmftxt ] = RST_writeFile( model, p, ux, uy, uz, x );
    
    % draw the result images
    RST_plot( model, fptxt, fuxtxt, fuytxt, fuztxt, fmftxt, fmhtxt, fmrtxt );
    
end