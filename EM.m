%EM
M=3;          % M����˹�ֲ����
N=600;        % ������
th=0.000001;  % ������ֵ
K=2;          % ����ά��
% ���������ݵĲ���
a_real =[2/3;1/6;1/6];%���ģ���л�ģ�͸�˹�ܶȺ�����Ȩ��
mu_real=[3 4 6;
         5 3 7];%��ֵ
cov_real(:,:,1)=[5 0;0 0.2];%Э����
cov_real(:,:,2)=[0.1 0;0 0.1];
cov_real(:,:,3)=[0.1 0;0 0.1];                    
%���ɷ��ϱ�׼����������(ÿһ��Ϊһ������)
x=[ mvnrnd( mu_real(:,1) , cov_real(:,:,1) , round(N*a_real(1)) )' ,...
    mvnrnd( mu_real(:,2) , cov_real(:,:,2) , round(N*a_real(2)) )' ,...
    mvnrnd( mu_real(:,3) , cov_real(:,:,3) , round(N*a_real(3)) )' ];
%��ʼ������
a=[1/3;1/3;1/3];
mu=[1 2 3;2 1 4];
cov(:,:,1)=[1 0;0 1];
cov(:,:,2)=[1 0;0 1];
cov(:,:,3)=[1 0;0 1];
t=inf;
while t>=th
    a_old  = a;
    mu_old = mu;
    cov_old= cov;     
    rznk_temp=zeros(M,N);
    for k=1:M
        for n=1:N
            %����P(x|mu_cm,cov_cm)
            rznk_temp(k,n)=exp(-1/2*(x(:,n)-mu(:,k))'*inv(cov(:,:,k))*(x(:,n)-mu(:,k)));
        end
        rznk_temp(k,:)=rznk_temp(k,:)/sqrt(det(cov(:,:,k)));
    end
    rznk_temp=rznk_temp*(2*pi)^(-K/2);
%E step
    %��rznk
    rznk=zeros(M,N);
    for n=1:N
        for k=1:M
            rznk(k,n)=a(k)*rznk_temp(k,n);
        end
        rznk(:,n)=rznk(:,n)/sum(rznk(:,n));
    end
% M step
    %��Nk
    nk=zeros(1,M);
    nk=sum(rznk');
   
    % ��a
    a=nk/N;
       
    % ��MU
    for k=1:M
        mu_k_sum=0;
        for n=1:N
            mu_k_sum=mu_k_sum+rznk(k,n)*x(:,n);
        end
        mu(:,k)=mu_k_sum/nk(k);
    end
   
    % ��COV  
    for k=1:M
        cov_k_sum=0;
        for n=1:N
            cov_k_sum=cov_k_sum+rznk(k,n)*(x(:,n)-mu(:,k))*(x(:,n)-mu(:,k))';
        end
        cov(:,:,k)=cov_k_sum/nk(k);
    end
    %��ֹ��������Ȩֵ����������ֵ��������Э������������ߵķ��������С��th(��ֵ)  
    t=max([norm(a_old(:)-a(:))/norm(a_old(:));norm(mu_old(:)-mu(:))/norm(mu_old(:));norm(cov_old(:)-cov(:))/norm(cov_old(:))]); 
end 
%���������Ƚ�
a_real
a
mu_real
mu
cov_real
cov