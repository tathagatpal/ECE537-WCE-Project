function[hdk_array,theta,Hrk_array,w_array]=channel_generation(k,m,n,r1,r2)
    hdk_array=cell(k,k); %k x k (each element m x 1)
    
    for i = 1:k
        for j = 1:k
            for x = 1:m
                hdk_temp = zeros(m,1);
                R = (r2-r1).*rand(1,1) + r1;
                lambda = 3*10^8/(10*10^6);
                hdk = fspl(R,lambda);
                hdk_temp(x,1) = hdk; %LoS coeff (Tx to RIS)
            end
            hdk_array{i,j} = hdk_temp;
        end
    end
    
    w_array=cell(k,k); %k x k (each element m x 1)
    
    for i = 1:k
        for j = 1:k
            for x = 1:m
                w_temp = zeros(m,1);
                R = (r2-r1).*rand(1,1) + r1;
                lambda = 3*10^8/(10*10^6);
                w = fspl(R,lambda);
                w_temp(x,1) = hdk; %LoS coeff (Tx to RIS)
            end
            w_array{i,j} = w_temp;
        end
    end
    
    a=randn(1,n);
    b=randn(1,n);
    theta=complex(a,b)';
    theta = bsxfun(@rdivide, theta, sqrt(sum(theta.*conj(theta), 2)));
    
    Hrk_array=cell(k,k);   %k x 1 (each element n x m)
    for i = 1:k
        for j = 1:k
            for x = 1:n
                for y = 1:m
                    hrk_temp = zeros(n,m);
                    R = abs((r2-r1).*rand(1,1) + r1);
                    lambda = 3*10^8/(10*10^6);
                    hrk = fspl(R,lambda);   %RIS to Rx
                    hrk_temp(x,y) = hrk;
                end
            end
            Hrk_array{i,j} = (hrk_temp);
        end
    end
    
end