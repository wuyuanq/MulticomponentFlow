% This is the RST_writeFile() function which outputs the results of 
% RST_singlePhaseFlow() to a series of solution files.

% Input parameters:
% model: the model
% p: the pressure
% ux: the velocities in the x dirction
% uy: the velocities in the y dirction
% uz: the velocities in the z dirction
% x: the mole fraction

% Return value:
% fptxt: the pressure file in the txt form
% fuxtxt: the x-dirction velocity file in the txt form
% fuytxt: the y-dirction velocity file in the txt form
% fuztxt: the z-dirction velocity file in the txt form
% fmftxt: the mole fraction files in the txt form

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on October 23rd, 2013

function [ fptxt, fuxtxt, fuytxt, fuztxt, fmftxt ] = RST_writeFile( model, p, ux, uy, uz, x )

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
    soludoc = model.soludoc;
    
    % define the names of the solution files
    fp = [soludoc, '/soln_1PhFlw_P_RSTo.m'];
    fux = [soludoc, '/soln_1PhFlw_Ux_RSTo.m'];
    fuy = [soludoc, '/soln_1PhFlw_Uy_RSTo.m'];
    fuz = [soludoc, '/soln_1PhFlw_Uz_RSTo.m'];
    fptxt = [soludoc, '/soln_1PhFlw_P_raw.txt'];
    fuxtxt = [soludoc, '/soln_1PhFlw_Ux_raw.txt'];
    fuytxt = [soludoc, '/soln_1PhFlw_Uy_raw.txt'];
    fuztxt = [soludoc, '/soln_1PhFlw_Uz_raw.txt'];
    fmf = [];
    for m = 1 : Nc
        fk = [soludoc, '/soln_1PhFlw_X', num2str(m), '_RSTo.m'];
        fmf = [fmf; fk];
    end
    cdfmf = cellstr(fmf);
    fmftxt = [];
    for m = 1 : Nc
        fk = [soludoc, '/soln_1PhFlw_X', num2str(m), '_raw.txt'];
        fmftxt = [fmftxt; fk];
    end
    cdfmftxt = cellstr(fmftxt);
    
    fid = fopen(fptxt, 'w');
    for k = 2 : nz+1
        for j = 2 : ny+1
             for i = 2 : nx+1
                fprintf(fid, '%10e\n', p(i,j,k));
            end
        end
    end
    fclose(fid);
    
    fid = fopen(fuxtxt, 'w');
    for k = 1 : nz 
        for j = 1 : ny
            for i = 1 : nx+1
                fprintf(fid, '%10e\n', ux(i,j,k));
            end
        end
    end
    fclose(fid);
    
    fid = fopen(fuytxt, 'w');
    for k = 1 : nz 
        for j = 1 : ny+1
            for i = 1 : nx
                fprintf(fid, '%10e\n', uy(i,j,k));
            end
        end
    end
    fclose(fid); 
    
    fid = fopen(fuztxt, 'w');
    for k = 1 : nz+1
        for j = 1 : ny
            for i = 1 : nx
                fprintf(fid, '%10e\n', uz(i,j,k));
            end
        end
    end
    fclose(fid); 
    
    for m = 1 : Nc
        fid = fopen(char(cdfmftxt(m)), 'w');
        for k = 1 : nz  
            for j = 1 : ny
                for i = 1 : nx
                    fprintf(fid, '%10e\n', x(m,i,j,k));
                end
            end
        end
        fclose(fid); 
    end
    
    fid = fopen(fp, 'w');
    fprintf(fid, 'p.type = ''cell-centered_data'';\n');
    fprintf(fid, 'p.simTime =      %.4f    ;\n', ts(nt+1));
    fprintf(fid, 'p.mesh.type = ''rectangular_mesh'';\n');
    fprintf(fid, 'p.mesh.xs = [');
    for i = 1 : nx
        fprintf(fid, '%.4f,\t', xs(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', xs(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'p.mesh.ys = [');
    for i = 1 : ny
        fprintf(fid, '%.4f,\t', ys(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', ys(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'p.mesh.zs = [');
    for i = 1 : nz
        fprintf(fid, '%.4f,\t', zs(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', zs(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'p.array = zeros(%d, %d, %d);\n', nx, ny, nz); 
    fprintf(fid, 'p.array(1:%d, 1:%d, 1:%d) = reshape( load(''%s''), [%d, %d, %d] );', nx, ny, nz, fptxt, nx, ny, nz); 
    fclose(fid);

    fid = fopen(fux, 'w');
    fprintf(fid, 'ux.type = ''face-centered_data'';\n');
    fprintf(fid, 'ux.simTime =      %.4f    ;\n', ts(nt+1));
    fprintf(fid, 'ux.mesh.type = ''rectangular_mesh'';\n');
    fprintf(fid, 'ux.mesh.xs = [');
    for i = 1 : nx
        fprintf(fid,'%.4f,\t', xs(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', xs(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'ux.mesh.ys = [');
    for i = 1 : ny
        fprintf(fid, '%.4f,\t', ys(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', ys(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'p.mesh.zs = [');
    for i = 1 : nz
        fprintf(fid, '%.4f,\t', zs(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', zs(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'ux.array = zeros(%d,%d,%d);\n', nx+1, ny, nz);
    fprintf(fid, 'ux.array(1:%d, 1:%d, 1:%d) = reshape( load(''%s''), [%d, %d, %d] );\n', nx+1, ny, nz, fuxtxt, nx+1, ny, nz); 
    fclose(fid);
    
    fid = fopen(fuy, 'w');
    fprintf(fid, 'uy.type = ''face-centered_data'';\n');
    fprintf(fid, 'uy.simTime =      %.4f    ;\n', ts(nt+1));
    fprintf(fid, 'uy.mesh.type = ''rectangular_mesh'';\n');
    fprintf(fid, 'uy.mesh.xs = [');
    for i = 1 : nx
        fprintf(fid, '%.4f,\t', xs(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', xs(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'uy.mesh.ys = [');
    for i = 1 : ny
        fprintf(fid, '%.4f,\t', ys(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', ys(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'p.mesh.zs = [');
    for i = 1 : nz
        fprintf(fid, '%.4f,\t', zs(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', zs(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'uy.array = zeros(%d,%d,%d);\n', nx, ny+1, nz);
    fprintf(fid, 'uy.array(1:%d, 1:%d, 1:%d) = reshape( load(''%s''), [%d, %d, %d] );\n', nx, ny+1, nz, fuytxt, nx, ny+1, nz); 
    fclose(fid);
    
    fid = fopen(fuz, 'w');
    fprintf(fid, 'uz.type = ''face-centered_data'';\n');
    fprintf(fid, 'uz.simTime =      %.4f    ;\n', ts(nt+1));
    fprintf(fid, 'uz.mesh.type = ''rectangular_mesh'';\n');
    fprintf(fid, 'uz.mesh.xs = [');
    for i = 1 : nx
        fprintf(fid, '%.4f,\t', xs(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', xs(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'uy.mesh.ys = [');
    for i = 1 : ny
        fprintf(fid, '%.4f,\t', ys(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', ys(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'p.mesh.zs = [');
    for i = 1 : nz
        fprintf(fid, '%.4f,\t', zs(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', zs(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'uz.array = zeros(%d,%d,%d);\n', nx, ny, nz+1);
    fprintf(fid, 'uz.array(1:%d, 1:%d, 1:%d) = reshape( load(''%s''), [%d, %d, %d] );\n', nx, ny, nz+1, fuztxt, nx, ny, nz+1); 
    fclose(fid);

    for m = 1 : Nc
        fid = fopen(char(cdfmf(m)), 'w');
        fprintf(fid, 'x%d.type = ''cell-centered_data'';\n',m);
        fprintf(fid, 'x%d.simTime =      %.4f    ;\n', m, ts(nt+1));
        fprintf(fid, 'x%d.mesh.type = ''rectangular_mesh'';\n', m);
        fprintf(fid, 'x%d.mesh.xs = [', m);
        for i = 1 : nx
            fprintf(fid, '%.4f,\t', xs(i));
        end
        i = i + 1;
        fprintf(fid, '%.4f\t', xs(i));
        fprintf(fid, '];\n');
        fprintf(fid, 'x%d.mesh.ys = [', m);
        for i = 1 : ny
            fprintf(fid, '%.4f,\t', ys(i));
        end
        i = i + 1;
        fprintf(fid, '%.4f\t', ys(i));
        fprintf(fid, '];\n');
        fprintf(fid, 'x%d.mesh.zs = [', m);
        for i = 1 : nz
            fprintf(fid, '%.4f,\t', zs(i));
        end
        i = i + 1;
        fprintf(fid, '%.4f\t', zs(i));
        fprintf(fid, '];\n');
        fprintf(fid, 'x%d.array = zeros(%d,%d,%d);\n', m, nx, ny, nz);
        fprintf(fid, 'x%d.array(1:%d, 1:%d, 1:%d) = reshape( load(''%s''), [%d, %d, %d] );\n', m, nx, ny, nz, char(cdfmftxt(m)), nx, ny, nz); 
        fclose(fid);
    end

end

