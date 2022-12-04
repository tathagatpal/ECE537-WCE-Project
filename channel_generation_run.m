k = 5;
m = 3; n = 4;
a = 50; %50 m
b = 10000;  %10000 m

%hdk_array,theta,Hrk_array,G_array --> in dB
[hdk_array,theta,Hrk_array]=channel_generation(k,m,n,a,b);


% hdk_array = db2mag(hdk_array);



% r1 = 50; %50 m
% r2 = 10000;  %10000 m
% 
% k = 3;
% m = 2;
% hdk_array=cell(k,k);
% 
% for r = 1:k
%     for c = 1:k
%         for i = 1:m
%                 hdk = zeros(m,1);
%                 R = (r2-r1).*rand(1,1) + r1;
%                 lambda = 3*10^8/(10*10^6);
%                 hdk(i,1) = fspl(R,lambda);
%         end
%         
%         hdk_array{r,c} = hdk;
%     end
% end
%     
