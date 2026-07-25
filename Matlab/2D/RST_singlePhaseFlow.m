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
    numDims = model.numDims; 
    Nc = model.Nc;
    nx = model.nx;
    ny = model.ny;
    nt = model.nt;
    xs = model.xs;
    ys = model.ys;
    ts = model.ts;
    poro = model.poro;
    Kxx = model.Kxx;
    Kyy = model.Kyy;   
    isDiriX = model.isDiriX;
    isDiriY = model.isDiriY;
    pBdryX = model.pBdryX;
    pBdryY = model.pBdryY;
    uBdryX = model.uBdryX;
    uBdryY = model.uBdryY; 
    xBdryX = model.xBdryX;
    xBdryY = model.xBdryY;
    pInit = model.pInit;
    xInit = model.xInit;
    gravX = model.gravX;
    gravY = model.gravY;
    src = model.src;
    soludoc = model.soludoc;
    
    if(numDims ~= 2)
        error('Such dimension simulation is not implemented yet!\n');
    end
    
    % create the new document to store the results
    if(exist(soludoc, 'dir') == 0)
        mkdir(soludoc);
    end
       
    % the file to store the issue mole fraction history
    fmhtxt = [soludoc, '/soln_1PhFlw_moleHistory.txt'];
    fmhtxtid = fopen(fmhtxt, 'w');
    [message, errnum] = ferror(fmhtxtid);
    if(errnum ~= 0)
        error(message, errnum);
    end
    
    % the file to store the left mole ratio history
    fmrtxt = [soludoc, '/soln_1PhFlw_moleRatioHistory.txt'];
    fmrtxtid = fopen(fmrtxt, 'w');
    [message, errnum] = ferror(fmrtxtid);
    if(errnum ~= 0)
        error(message, errnum);
    end
  
    % rectify the directions of the velocities
    uBdryX(1, 1:end) = -uBdryX(1, 1:end);
    uBdryY(1:end ,1) = -uBdryY(1:end, 1);
    
    % define the pressures and initialize them
    p = zeros(ny+2, nx+2);
    p(1, 2:nx+1) = pBdryY(1:end, 1);
    p(ny+2, 2:nx+1) = pBdryY(1:end, 2);
    p(2:ny+1, 1) = pBdryX(1, 1:end);
    p(2:ny+1, nx+2) = pBdryX(2, 1:end);
    for i = 2 : ny+1
        for j = 2 : nx+1
            p(i,j) = pInit(j-1,i-1);
        end
    end
    
    % define the fluid velocities and initialize them
    ux = zeros(ny, nx+1);
    uy = zeros(ny+1, nx);
    ux(1:end, 1) = uBdryX(1, 1:end);
    ux(1:end, nx+1) = uBdryX(2, 1:end);
    uy(1, 1:end) = uBdryY(1:end, 1);
    uy(ny+1, 1:end) = uBdryY(1:end, 2);
        
    % define the mole fraction of each component in the cell and initialize them
    x = zeros(Nc, ny, nx);
    for m = 1 : Nc
        for i = 1 : ny
            for j = 1 : nx
                x(m, i, j) = xInit(m, j, i);
            end
        end
    end
    
    % define the mole fraction of each component on the bars and initialize them
    xbarx = zeros(Nc, ny, nx+1);
    xbary = zeros(Nc, ny+1, nx);
    
    % define the mobility*(fluid molar density) on the bars
    lambdax = zeros(ny, nx+1);
    lambday = zeros(ny+1, nx);
    
    % define the mole density
    xi = zeros(ny, nx);
    
    % define the mole density on the bars
    xibarx = zeros(ny, nx+1);
    xibary = zeros(ny+1, nx);
    
    % define the mass density in the cell
    rho = zeros(ny, nx);
    
    % define the mass density on the bars
    rhobarx = zeros(ny, nx+1);
    rhobary = zeros(ny+1, nx);
    
    % define the derivation of fluid density to pressure 
    deri_xi_p = zeros(ny, nx);
    
    % define the viscosity
    mu = zeros(ny, nx);   
    
    % define the permeabilities on the bars
    Kxxbar = zeros(ny, nx+1);
    Kyybar = zeros(ny+1, nx);
    
    % initialize Kxxbar using harmonic weighting method
    Kxxbar(1:ny,1) = Kxx(1,1:ny);
    Kxxbar(1:ny,nx+1) = Kxx(nx, 1:ny);
    for i = 1 : ny
        for j = 2 : nx  
            ltotal = xs(j+1) - xs(j-1);
            lleft = xs(j) - xs(j-1);
            lright = xs(j+1) - xs(j);
            Kxxbar(i,j) = ltotal / (lleft/Kxx(j-1,i)+lright/Kxx(j,i));
        end
    end
    
    % initialize Kyybar using harmonic weighting method
    Kyybar(1,1:nx) = Kyy(1:nx,1);
    Kyybar(ny+1,1:nx) = Kyy(1:nx, ny);
    for i = 2 : ny
        for j = 1 : nx  
            ltotal = ys(i+1) - ys(i-1);
            lup = ys(i+1) - ys(i);
            ldown = ys(i) - ys(i-1);
            Kyybar(i,j) = ltotal / (ldown/Kyy(j, i-1)+lup/Kyy(j,i));
        end
    end
    
    % define the mole of each component in each cell
    moleincell = zeros(Nc,ny,nx);
    
    % define the left mole amount of each component
    leftmole = zeros(Nc, 1);
    
    % define the coefficient matrices A and b
    row = [];
    col = [];
    value = [];
    A = sparse(row, col, value, nx*ny, nx*ny);
    b = zeros(nx*ny, 1);
    % define the new molar density of each component
    c = zeros(Nc, ny, nx);
    firsttime = 1;
    totalmole = 0.0;
    
    % time iteration 
    for t = 2 : nt+1
        
        t    
        % compute the PR EOS and viscosity
        for j = 1 : nx
            for i = 1 : ny
                [ xi(i,j), rho(i,j), deri_xi_p(i,j) ] = RST_PREOS( model, x(1:end,i,j), p(i+1,j+1) ); 
                [ mu(i,j) ] = RST_viscosity( model, x(1:end,i,j), xi(i,j), p(i+1,j+1) ); 
            end
        end

        for j = 1 : nx
            for i = 1 : ny
                for m = 1 : Nc
                    moleincell(m,i,j) = xi(i,j)*x(m,i,j)*(xs(j+1)-xs(j))*(ys(i+1)-ys(i))*poro(j,i);
                end 
            end 
        end 

        leftmole(:,1) = 0;
        for j = 1 : nx
            for i = 1 : ny
                for m = 1 : Nc
                    leftmole(m) = leftmole(m) + moleincell(m,i,j);
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

        fprintf(fmhtxtid, '%10e\n', (totalmole-totaldesiredleftmole)/totalmole);
        
        % print the mole ratio of desired components to component 1 in
        % the well
        fprintf(fmrtxtid, '%10e\n', totaldesiredleftmole/leftmole(1));
        
        % compute the time step
        %maxux = max(max(abs(ux)));
        %maxuy = max(max(abs(uy)));
        %deltax = xs(2)-xs(1);
        %deltay = ys(2)-ys(1);
        %timestep = floor(poro(1,1)*deltax*deltay/(maxux*deltay+maxuy*deltax));
        %ts(t+1) = ts(t) + timestep;
        timestep = ts(2) - ts(1);
        
        % compute lambda = mobility*(fluid molar density) on the bars,
        % using single-point upstream weighting
        for i = 1 : ny
            if(ux(i,1) > 0)
                [ xibarx(i,1), rhobarx(i,1), ~ ] = RST_PREOS( model, xBdryX(1:end,1,i), p(i+1,2) ); 
                [ mutemp ] = RST_viscosity( model, xBdryX(1:end,1,i), xibarx(i,1), p(i+1,2) );           
                lambdax(i,1) = Kxxbar(i,1)*xibarx(i,1)/mutemp;
                xbarx(1:Nc,i,1) = xBdryX(1:Nc,1,i);
            elseif(ux(i,1) < 0)
                rhobarx(i,1) = rho(i,1);
                xibarx(i,1) = xi(i,1);
                lambdax(i,1) = Kxxbar(i,1)*xibarx(i,1)/mu(i,1);  
                xbarx(1:Nc,i,1) = x(1:Nc,i,1);
            else
                if(isDiriX(1,i) == 1)
                    rhobarx(i,1) = rho(i,1);
                    xibarx(i,1) = xi(i,1);
                    lambdax(i,1) = Kxxbar(i,1)*xibarx(i,1)/mu(i,1); 
                    xbarx(1:Nc,i,1) = x(1:Nc,i,1);
                end
            end
        end
        
        for i = 1 : ny
            if(ux(i,nx+1) < 0)
                [ xibarx(i,nx+1), rhobarx(i,nx+1), ~ ] = RST_PREOS( model, xBdryX(1:end,2,i), p(i+1,nx+1) ); 
                [ mutemp ] = RST_viscosity( model, xBdryX(1:end,2,i), xibarx(i,nx+1), p(i+1,nx+1) );  
                lambdax(i,nx+1) = Kxxbar(i,nx+1)*xibarx(i,nx+1)/mutemp;
                xbarx(1:Nc,i,nx+1) = xBdryX(1:Nc,2,i);
            elseif(ux(i,nx+1) > 0)
                rhobarx(i,nx+1) = rho(i,nx);
                xibarx(i,nx+1) = xi(i,nx);  
                lambdax(i,nx+1) = Kxxbar(i,nx+1)*xibarx(i,nx+1)/mu(i,nx); 
                xbarx(1:Nc,i,nx+1) = x(1:Nc,i,nx);
            else
                if(isDiriX(2,i) == 1)
                    rhobarx(i,nx+1) = rho(i,nx);
                    xibarx(i,nx+1) = xi(i,nx);  
                    lambdax(i,nx+1) = Kxxbar(i,nx+1)*xibarx(i,nx+1)/mu(i,nx);
                    xbarx(1:Nc,i,nx+1) = x(1:Nc,i,nx);
                end
            end
        end
                    
        for j = 2 : nx  
            for i = 1 : ny             
                if(ux(i,j) > 0) 
                    rhobarx(i,j) = rho(i,j-1);
                    xibarx(i,j) = xi(i,j-1); 
                    lambdax(i,j) = Kxxbar(i,j)*xibarx(i,j)/mu(i,j-1); 
                    xbarx(1:Nc,i,j) = x(1:Nc,i,j-1);
                else
                    rhobarx(i,j) = rho(i,j);
                    xibarx(i,j) = xi(i,j);
                    lambdax(i,j) = Kxxbar(i,j)*xibarx(i,j)/mu(i,j);  
                    xbarx(1:Nc,i,j) = x(1:Nc,i,j);
                end
            end
        end
        
        for j = 1 : nx
            if(uy(1,j) > 0)
                [ xibary(1,j), rhobary(1,j), ~ ] = RST_PREOS( model, xBdryY(1:end,j,1), p(2,j+1) ); 
                [ mutemp ] = RST_viscosity( model, xBdryY(1:end,j,1), xibary(1,j), p(2,j+1) );    
                lambday(1,j) = Kyybar(1,j)*xibary(1,j)/mutemp;
                xbary(1:Nc,1,j) = xBdryY(1:Nc,j,1);
            elseif(uy(1,j) < 0)
                rhobary(1,j) = rho(1,j);
                xibary(1,j) = xi(1,j);
                lambday(1,j) = Kyybar(1,j)*xibary(1,j)/mu(1,j); 
                xbary(1:Nc,1,j) = x(1:Nc,1,j);
            else
                if(isDiriY(j,1) == 1)
                    rhobary(1,j) = rho(1,j);
                    xibary(1,j) = xi(1,j);
                    lambday(1,j) = Kyybar(1,j)*xibary(1,j)/mu(1,j); 
                    xbary(1:Nc,1,j) = x(1:Nc,1,j);
                end
            end
        end

        for j = 1 : nx
            if(uy(ny+1,j) < 0)
                [ xibary(ny+1,j), rhobary(ny+1,j), ~ ] = RST_PREOS( model, xBdryY(1:end,j,2), p(ny+1,j+1) ); 
                [ mutemp ] = RST_viscosity( model, xBdryY(1:end,j,2), xibary(ny+1,j), p(ny+1,j+1) );
                lambday(ny+1,j) = Kyybar(ny+1,j)*xibary(ny+1,j)/mutemp;
                xbary(1:Nc,ny+1,j) = xBdryY(1:Nc,j,2);
            elseif(uy(ny+1,j) > 0)
                rhobary(ny+1,j) = rho(ny,j);
                xibary(ny+1,j) = xi(ny,j);
                lambday(ny+1,j) = Kyybar(ny+1,j)*xibary(ny+1,j)/mu(ny,j);
                xbary(1:Nc,ny+1,j) = x(1:Nc,ny,j);
            else
                if(isDiriY(j,2) == 1)
                    rhobary(ny+1,j) = rho(ny,j);
                    xibary(ny+1,j) = xi(ny,j);
                    lambday(ny+1,j) = Kyybar(ny+1,j)*xibary(ny+1,j)/mu(ny,j);
                    xbary(1:Nc,ny+1,j) = x(1:Nc,ny,j);
                end
            end
        end
                    
        for j = 1 : nx
            for i = 2 : ny  
                if(uy(i,j) > 0)
                    rhobary(i,j) = rho(i-1,j);
                    xibary(i,j) = xi(i-1,j);
                    lambday(i,j) = Kyybar(i,j)*xibary(i,j)/mu(i-1,j);   
                    xbary(1:Nc,i,j) = x(1:Nc,i-1,j);
                else
                    rhobary(i,j) = rho(i,j);
                    xibary(i,j) = xi(i,j); 
                    lambday(i,j) = Kyybar(i,j)*xibary(i,j)/mu(i,j); 
                    xbary(1:Nc,i,j) = x(1:Nc,i,j);
                end
            end
        end
        
        % compute the coefficiencies of the pressure equations   
        for j = 2 : nx+1
            for i = 2 : ny+1

                r = (i-2)*nx + (j-1); % r is the index of the current cell
                
                xedge = xs(j) - xs(j-1); % x-edge of the current cell
                yedge = ys(i) - ys(i-1); % y-edge of the current cell
            
                % ledge means the x-edge of the left cell 
                if(j ~= 2)
                    ledge = xs(j-1) - xs(j-2);
                else
                    ledge = 0;
                end
            
                % redge means the x-edge of the right cell
                if(j ~= nx+1)
                    redge = xs(j+1) - xs(j);    
                else
                    redge = 0;
                end
            
                % uedge means the y-edge of the up cell
                if(i ~= ny+1)
                    uedge = ys(i+1) - ys(i);
                else
                    uedge = 0;
                end
            
                % dedge means the y-edge of the down cell
                if(i ~= 2)
                    dedge = ys(i-1) - ys(i-2);
                else
                    dedge = 0;
                end

                % compute the coefficiencies of the pressure equations
                up = -2*lambday(i,j-1)/yedge/(yedge+uedge); % up means the coefficiency of the pressure of the up cell
                down = -2*lambday(i-1,j-1)/yedge/(yedge+dedge); % down means the coefficiency of the pressure of the down cell
                left = -2*lambdax(i-1,j-1)/xedge/(xedge+ledge); % left means the coefficiency of the pressure of the left cell
                right = -2*lambdax(i-1,j)/xedge/(xedge+redge); % right means the coefficiency of the pressure of the right cell

                cp = poro(j-1, i-1)*deri_xi_p(i-1, j-1); % suppose porosity was constant here
                
                q = 0;
                for m = 1 : Nc
                    q = q + src(m, j-1, i-1);
                end
                
                A(r,r) = cp/timestep;
                b(r,1) = q + cp*p(i,j)/timestep;
                
                if((i == ny+1) && (isDiriY(j-1,2) == 0)) % Neumann boundary
                    b(r,1) = b(r,1) - uBdryY(j-1,2)*xibary(i,j-1)/yedge;
                elseif((i == ny+1) && (isDiriY(j-1,2) == 1)) % Dirichlet boundary
                    b(r,1) = b(r,1) - up*p(i+1,j) - lambday(i,j-1)*rhobary(i,j-1)*gravY/yedge;
                    A(r,r) = A(r,r) - up;
                else
                    A(r,r+nx) = up;
                    A(r,r) = A(r,r) - up;
                    b(r,1) = b(r,1) - lambday(i,j-1)*rhobary(i,j-1)*gravY/yedge;
                end
                
                if((i == 2) && (isDiriY(j-1,1) == 0))
                    b(r,1) = b(r,1) + uBdryY(j-1,1)*xibary(i-1,j-1)/yedge;
                elseif((i == 2) && (isDiriY(j-1,1) == 1))
                    b(r,1) = b(r,1) - down*p(i-1,j) + lambday(i-1,j-1)*rhobary(i-1,j-1)*gravY/yedge;
                    A(r,r) = A(r,r) - down;
                else
                    A(r,r-nx) = down;
                    A(r,r) = A(r,r) - down;
                    b(r,1) = b(r,1) + lambday(i-1,j-1)*rhobary(i-1,j-1)*gravY/yedge;
                end
                
                if((j == 2) && (isDiriX(1,i-1) == 0))
                    b(r,1) = b(r,1) + uBdryX(1,i-1)*xibarx(i-1,j-1)/xedge;
                elseif((j == 2) && (isDiriX(1,i-1) == 1))
                    b(r,1) = b(r,1) - left*p(i,j-1) + lambdax(i-1,j-1)*rhobarx(i-1,j-1)*gravX/xedge;
                    A(r,r) = A(r,r) - left;
                else
                    A(r,r-1) = left;
                    A(r,r) = A(r,r) - left;
                    b(r,1) = b(r,1) + lambdax(i-1,j-1)*rhobarx(i-1,j-1)*gravX/xedge;
                end
                
                if((j == nx+1) && (isDiriX(2,i-1) == 0))
                    b(r,1) = b(r,1) - uBdryX(2,i-1)*xibarx(i-1,j)/xedge;  
                elseif((j == nx+1) && (isDiriX(2,i-1) == 1))
                    b(r,1) = b(r,1) - right*p(i,j+1) - lambdax(i-1,j)*rhobarx(i-1,j)*gravX/xedge;
                    A(r,r) = A(r,r) - right;
                else
                    A(r,r+1) = right;
                    A(r,r) = A(r,r) - right;
                    b(r,1) = b(r,1) - lambdax(i-1,j)*rhobarx(i-1,j)*gravX/xedge;
                end
            end
        end
    
        % compute the new pressures
        xx = A\b;
        for j = 2 : nx+1
            for i = 2 : ny+1 
                p(i,j) = xx((i-2)*nx+(j-1),1);
            end
        end
        
         % compute the new velocities in x direction  
        for i = 1 : ny
            if(isDiriX(1,i) == 1)
                ux(i,1) = -lambdax(i,1)/xibarx(i,1) * ((p(i+1,2)-p(i+1,1))*2/(xs(2)-xs(1)) - rhobarx(i,1)*gravX);
            else
                ux(i,1) = uBdryX(1,i);
            end
        end
        
        for i = 1 : ny
            if(isDiriX(2,i) == 1)
                ux(i,nx+1) = -lambdax(i,nx+1)/xibarx(i,nx+1) * ((p(i+1,nx+2)-p(i+1,nx+1))*2/(xs(nx+1)-xs(nx)) - rhobarx(i,nx+1)*gravX);
            else
                ux(i,nx+1) = uBdryX(2,i);
            end
        end
        
        for i = 1 : ny
            for j = 2 : nx
                ux(i,j) = -lambdax(i,j)/xibarx(i,j) * ((p(i+1,j+1)-p(i+1,j))*2/(xs(j+1)-xs(j-1)) - rhobarx(i,j)*gravX);
            end
        end
    
        % compute the new velocities in y direction 
        for j = 1 : nx
            if(isDiriY(j,1) == 1)
                uy(1,j) = -lambday(1,j)/xibary(1,j) * ((p(2,j+1)-p(1,j+1))*2/(ys(2)-ys(1)) - rhobary(1,j)*gravY);
            else
                uy(1,j) = uBdryY(j,1);
            end
        end
        
        for j = 1 : nx
            if(isDiriY(j,2) == 1)
                uy(ny+1,j) = -lambday(ny+1,j)/xibary(ny+1,j) * ((p(ny+2,j+1)-p(ny+1,j+1))*2/(ys(ny+1)-ys(ny)) - rhobary(ny+1,j)*gravY);
            else
                uy(ny+1,j) = uBdryY(j,2);
            end
        end
                    
        for j = 1 : nx
            for i = 2 : ny 
                uy(i,j) = -lambday(i,j)/xibary(i,j) * ((p(i+1,j+1)-p(i,j+1))*2/(ys(i+1)-ys(i-1)) - rhobary(i,j)*gravY);
            end
        end
        
        % compute the new mole fraction of each component
        for j = 1 : nx
            for i = 1 : ny  
                for m = 1 : Nc           
                    div = (xbarx(m,i,j+1)*ux(i,j+1)*xibarx(i,j+1) - xbarx(m,i,j)*ux(i,j)*xibarx(i,j))/(xs(j+1)-xs(j)) + ...
                        (xbary(m,i+1,j)*uy(i+1,j)*xibary(i+1,j) - xbary(m,i,j)*uy(i,j)*xibary(i,j))/(ys(i+1)-ys(i));
                    c(m,i,j) = (src(m,j,i)-div)*timestep/poro(j,i) + x(m,i,j)*xi(i,j);
                    if(c(m,i,j)<0)
                        c(m,i,j)
                        error('Please tune the time step.');
                    end
                end
                ctotal = 0;
                for m = 1 : Nc 
                    ctotal = ctotal + c(m,i,j);
                end
                for m = 1 : Nc 
                    x(m,i,j) = c(m,i,j)/ctotal;
                end
            end
        end  

    end
    
    % close the mole history file
    fclose(fmhtxtid);
    
    % close the mole rate history file
    fclose(fmrtxtid);
    
    % write the results to the solution files
    [ fptxt, fuxtxt, fuytxt, fmftxt ] = RST_writeFile( model, p, ux, uy, x );
    
    % draw the result images
    RST_plot( model, fptxt, fuxtxt, fuytxt, fmftxt, fmhtxt, fmrtxt );
    
end