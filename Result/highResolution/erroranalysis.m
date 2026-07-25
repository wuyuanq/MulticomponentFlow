
nx = 40;
ny = 80;
location = '/Users/wuy/Dropbox/Research/MulticomponentFlow_results/highResolution/';

ftxt{1} = [location,'case4-40*80/soln_cmp2PhFlw_x1_p0_raw.txt'];
ftxt{2} = [location,'case4-80*160/soln_cmp2PhFlw_x1_p0_raw.txt'];
ftxt{3} = [location,'case4-160*320/soln_cmp2PhFlw_x1_p0_raw.txt'];

% load fields
for n = 1 : 2
    field{n} = zeros(ny*2^(n-1), nx*2^(n-1));
    temp = load(ftxt{n});
    c = 0;
    for j = 1 : ny*2^(n-1)
        for i = 1 : nx*2^(n-1) 
            c = c + 1;
            field{n}(j,i) = temp(c);
        end
    end
end

field{3} = zeros(ny*2^2, nx*2^2);
temp = load(ftxt{3});
c = 0;
for j = 1 : ny*2^2
    for i = 1 : nx*2^2 
        c = c + 1;
        field{3}(j,i) = temp(c);
    end
end

norm = 0;
for j = 1 : ny*2^2    
    for i = 1 : nx*2^2 
        norm = norm + field{3}(j,i)^2;
    end
end
norm = sqrt(norm/(ny*2^2*nx*2^2));

% average and error
for n = 1 : 2
    convfield{n} = zeros(ny*2^(n-1), nx*2^(n-1));
    error{n} = zeros(ny*2^(n-1), nx*2^(n-1));
    aveerror = 0;
    for j = 1 : ny*2^(n-1)
        for i = 1 : nx*2^(n-1) 
            convfield{n}(j,i) = 0;
            for nj = 1 : 4/2^(n-1)
                for ni = 1 : 4/2^(n-1)
                    convfield{n}(j,i) = convfield{n}(j,i) + field{3}((j-1)*4/2^(n-1)+nj,(i-1)*4/2^(n-1)+ni);
                end
            end
            convfield{n}(j,i) = convfield{n}(j,i)/(4/2^(n-1))^2;
            error{n}(j,i) = abs(convfield{n}(j,i)-field{n}(j,i));
            aveerror = aveerror + error{n}(j,i)^2;
        end
    end
    aveerror = sqrt(aveerror/(nx*2^(n-1)*ny*2^(n-1)))/norm;
    disp(['The L2 norm error of grid ', num2str(n), ' is ', num2str(aveerror)]);
end


