clc
clear all
run=1;
scaling=10^6;
final=[];

for runs=1:run
    k=4;
    n=100;
    m=30;
    r1=50;
    r2=10000;
   
    [hdk_array,theta,Hrk_array,w_array]=channel_generation(k,m,n,r1,r2);
    net_sumrate=sumrate_calc(hdk_array,theta,Hrk_array,k,w_array);
    final=[final net_sumrate];
    betas=cell(k,1);
    alphas=zeros(k,1);
    
    for iterates=1:50
            for i=1:k
                sum4=0;
                for j=1:k
        %             if(j~=i)
                    hk_=abs((hdk_array{i,j}'+ theta'*Hrk_array{i})*w_array{i,j})^2;
                    sum4=sum4+hk_;
        %             end
                end
                betai=(sqrt((1+alphas(i))) * (hdk_array{i,i}'+ theta'*Hrk_array{i}))*w_array{i,i}/sum4;
                betas{i}=betai;
                hi=hdk_array{i,i}+Hrk_array{i}'*theta;
                kaii=real(conj(betas{i})*hi'*w_array{i,i});
                alphai=(kaii^2+kaii*sqrt(kaii^2+4))/2;
                alphas(i)=alphai;
            end
           %U and V calculation
           sum_aik=0;   
            for i=1:k
                aik=Hrk_array{i}*w_array{i,i};
                aik=aik*aik';
                sum_aik=sum_aik+aik;
            end
            
            U=0;
            for i=1:k
                U=U+(norm(betas{i})^2)*sum_aik;
            end
            
            v=0;
            sum_ab=0;
            for i=1:k
                aik=Hrk_array{i}*w_array{i,i};
                bik=hdk_array{i,i}'*w_array{i,i};
                sum_ab=sum_ab+conj(bik)*aik;
            end
            v=0;
            for i=1:k
                aik=Hrk_array{i}*w_array{i,i};
                v=v+(sqrt((1+alphas(i)))*conj(betas{i})*aik)-(abs(betas{i}))^2*sum_ab;
            end
            
            phi=angle(theta);
            oof=U*exp(j.*phi)-v;
            oof2=-j.*exp(-j.*phi);
            grad=2*real(oof.*oof2);
            
            diagtheta=diag(theta');
              Ltheta=SCA_phi_step_para(U,v,n,diagtheta)*scaling;
%               Ltheta=length(theta);
              rho0=1/(Ltheta)*100;
              rhom=0.95;
              m=0;
              sig=0.4;
              while(1)
    %             step=1/(start);
                rho=rho0*rhom^m;
                phi1=phi-grad*rho;
                f1=exp(j.*phi)'*U*exp(j.*phi)-2*real(v'*exp(j.*phi));
                f2=exp(j.*phi1)'*U*exp(j.*phi1)-2*real(v'*exp(j.*phi1));
                if(f1-f2>(-grad'*grad)*sig*rho)
                    break
                end
                if(rho<1/Ltheta/10)
                    break
                end
                m=m+1;
    %             start=start-0.2*start;
              end
    %         step=100000;
%             phi2=phi-grad*step;
            phi2=phi-grad*rho;
%                     grad
%                 phi2=phi-grad;
            theta2=exp(j.*phi2);

            for i=1:k
                sum3=0;
                for j=1:k
                    if(j~=i)
                        sum3=sum3+(abs(hdk_array{i,j}'+theta'*Hrk_array{j})*w_array{i,j})^2;
                    end
                end
                gamma=(abs(hdk_array{i,i}'+theta'*Hrk_array{i})* w_array{i,i})^2/sum3;
                if(i==1)
                sumrate1=log2(1+gamma);
                end
                if(i==2)
                sumrate2=log2(1+gamma);
                end
                if(i==3)
                sumrate3=log2(1+gamma);
                end
                if(i==4)
                sumrate4=log2(1+gamma);
                end

            end

            new_net=(sumrate1+sumrate2+sumrate3+sumrate4);

            if(new_net>net_sumrate)
                net_sumrate=new_net;
            end
%             net_sumrate=new_net;
            phi=phi2;
            theta=exp(j.*phi);
            final=[final net_sumrate];   
        
            
            

    end
end

hold on    
plot((final)/run,"linewidth",1.5)
xlabel("Iterations")
ylabel("Sum Rate")
