function [net_sumrate]=sumrate_calc(hdk_array,theta,Hrk_array,k,w_array)%,sumrate1,sumrate2,sumrate3,sumrate4)
sumrate1=[];
sumrate2=[];
sumrate3=[];
sumrate4=[];

for i=1:k
            sum3=0;
            for j=1:k
                if(j~=i)
%                     hdk_array(i,j)
% size(hdk_array{i,j}')
% size(theta')
% size(Hrk_array{j})
% size(w_array{i,j})
                    sum3=sum3+(abs(hdk_array{i,j}'+theta'*Hrk_array{j})*w_array{i,j})^2;
                end
            end
            gamma=(abs(hdk_array{i,i}'+theta'*Hrk_array{i})* w_array{i,i})^2/sum3;
            if(i==1)
            sumrate1=[sumrate1 log2(1+gamma)];
            end
            if(i==2)
            sumrate2=[sumrate2 log2(1+gamma)];
            end
            if(i==3)
            sumrate3=[sumrate3 log2(1+gamma)];
            end
            if(i==4)
            sumrate4=[sumrate4 log2(1+gamma)];
            end

    end

    net_sumrate=norm((sumrate1+sumrate2+sumrate3+sumrate4));
end