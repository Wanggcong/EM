%EM
M=3;          % M个高斯分布混合
N=600;        % 样本数
th=0.000001;  % 收敛阈值
K=2;          % 样本维数
% 待生成数据的参数
a_real =[2/3;1/6;1/6];%混合模型中基模型高斯密度函数的权重
mu_real=[3 4 6;
         5 3 7];%均值
cov_real(:,:,1)=[5 0;0 0.2];%协方差
cov_real(:,:,2)=[0.1 0;0 0.1];
cov_real(:,:,3)=[0.1 0;0 0.1];                    
%生成符合标准的样本数据(每一列为一个样本)
x=[ mvnrnd( mu_real(:,1) , cov_real(:,:,1) , round(N*a_real(1)) )' ,...
    mvnrnd( mu_real(:,2) , cov_real(:,:,2) , round(N*a_real(2)) )' ,...
    mvnrnd( mu_real(:,3) , cov_real(:,:,3) , round(N*a_real(3)) )' ];
%初始化参数
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
            %计算P(x|mu_cm,cov_cm)
            rznk_temp(k,n)=exp(-1/2*(x(:,n)-mu(:,k))'*inv(cov(:,:,k))*(x(:,n)-mu(:,k)));
        end
        rznk_temp(k,:)=rznk_temp(k,:)/sqrt(det(cov(:,:,k)));
    end
    rznk_temp=rznk_temp*(2*pi)^(-K/2);
%E step
    %求rznk
    rznk=zeros(M,N);
    for n=1:N
        for k=1:M
            rznk(k,n)=a(k)*rznk_temp(k,n);
        end
        rznk(:,n)=rznk(:,n)/sum(rznk(:,n));
    end
% M step
    %求Nk
    nk=zeros(1,M);
    nk=sum(rznk');
   
    % 求a
    a=nk/N;
       
    % 求MU
    for k=1:M
        mu_k_sum=0;
        for n=1:N
            mu_k_sum=mu_k_sum+rznk(k,n)*x(:,n);
        end
        mu(:,k)=mu_k_sum/nk(k);
    end
   
    % 求COV  
    for k=1:M
        cov_k_sum=0;
        for n=1:N
            cov_k_sum=cov_k_sum+rznk(k,n)*(x(:,n)-mu(:,k))*(x(:,n)-mu(:,k))';
        end
        cov(:,:,k)=cov_k_sum/nk(k);
    end
    %终止条件，让权值的增量，均值的增量，协方差的增量三者的范数最大者小于th(阀值)  
    t=max([norm(a_old(:)-a(:))/norm(a_old(:));norm(mu_old(:)-mu(:))/norm(mu_old(:));norm(cov_old(:)-cov(:))/norm(cov_old(:))]); 
end 
%输出结果并比较
a_real
a
mu_real
mu
cov_real
cov