% This is the RST_writeFile() function which outputs the results of 
% RST_singlePhaseFlow() to a series of solution files.

% Input parameters:
% model: the model
% p: the pressure
% ux: the velocities in the x dirction
% uy: the velocities in the y dirction
% x: the mole fraction

% Return value:
% fptxt: the pressure file in the txt form
% fuxtxt: the x-dirction velocity file in the txt form
% fuytxt: the y-dirction velocity file in the txt form
% fmftxt: the mole fraction files in the txt form

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on October 23rd, 2012

function [ fptxt, fuxtxt, fuytxt, fmftxt ] = RST_writeFile( model, p, ux, uy, x )

    % simplify the symbols of the model
    Nc = model.Nc;
    nx = model.nx;
    ny = model.ny;
    nt = model.nt;
    xs = model.xs;
    ys = model.ys;
    ts = model.ts;
    soludoc = model.soludoc;
    
    % define the names of the solution files
    fp = [soludoc, '/soln_1PhFlw_P_RSTo.m'];
    fux = [soludoc, '/soln_1PhFlw_Ux_RSTo.m'];
    fuy = [soludoc, '/soln_1PhFlw_Uy_RSTo.m'];
    fptxt = [soludoc, '/soln_1PhFlw_P_raw.txt'];
    fuxtxt = [soludoc, '/soln_1PhFlw_Ux_raw.txt'];
    fuytxt = [soludoc, '/soln_1PhFlw_Uy_raw.txt'];
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
    [message, errnum] = ferror(fid);
    if(errnum ~= 0)
        error(message, errnum);
    end
    for j = 2 : nx+1
        for i = 2 : ny+1
            fprintf(fid, '%10e\n', p(i,j));
        end
    end
    fclose(fid);
    
    fid = fopen(fuxtxt, 'w');
    [message, errnum] = ferror(fid);
    if(errnum ~= 0)
        error(message, errnum);
    end
    for j = 1 : nx+1
        for i = 1 : ny
            fprintf(fid, '%10e\n', ux(i,j));
        end
    end
    fclose(fid);
    
    fid = fopen(fuytxt, 'w');
    [message, errnum] = ferror(fid);
    if(errnum ~= 0)
        error(message, errnum);
    end
    for j = 1 : nx
        for i = 1 : ny+1
            fprintf(fid, '%10e\n', uy(i,j));
        end
    end
    fclose(fid); 
    
    for m = 1 : Nc
        fid = fopen(char(cdfmftxt(m)), 'w');
        [message, errnum] = ferror(fid);
        if(errnum ~= 0)
            error(message, errnum);
        end
        for j = 1 : nx
            for i = 1 : ny
                fprintf(fid, '%10e\n', x(m,i,j));
            end
        end
        fclose(fid); 
    end
    
    fid = fopen(fp, 'w');
    [message, errnum] = ferror(fid);
    if(errnum ~= 0)
        error(message, errnum);
    end
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
    fprintf(fid, 'p.array = zeros(%d, %d);\n', nx, ny); 
    fprintf(fid, 'Pw.array(1:%d, 1:%d) = reshape( load(''%s''), [%d, %d] );', nx, ny, fptxt, nx, ny); 
    fclose(fid);

    fid = fopen(fux, 'w');
    [message, errnum] = ferror(fid);
    if(errnum ~= 0)
        error(message, errnum);
    end
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
    fprintf(fid, 'ux.array = zeros(%d,%d);\n', nx+1, ny);
    fprintf(fid, 'ux.array(1:%d, 1:%d) = reshape( load(''%s''), [%d, %d] );\n', nx+1, ny, fuxtxt, nx+1, ny); 
    fclose(fid);
    
    fid = fopen(fuy, 'w');
    [message, errnum] = ferror(fid);
    if(errnum ~= 0)
        error(message, errnum);
    end
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
    fprintf(fid, 'uy.array = zeros(%d,%d);\n', nx, ny+1);
    fprintf(fid, 'uy.array(1:%d, 1:%d) = reshape( load(''%s''), [%d, %d] );\n', nx, ny+1, fuytxt, nx, ny+1); 
    fclose(fid);

    for m = 1 : Nc
        fid = fopen(char(cdfmf(m)), 'w');
        [message, errnum] = ferror(fid);
        if(errnum ~= 0)
            error(message, errnum);
        end
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
        fprintf(fid, 'x%d.array = zeros(%d,%d);\n', m, nx, ny);
        fprintf(fid, 'x%d.array(1:%d, 1:%d) = reshape( load(''%s''), [%d, %d] );\n', m, nx, ny, char(cdfmftxt(m)), nx, ny); 
        fclose(fid);
    end
    
end

