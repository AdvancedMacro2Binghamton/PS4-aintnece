% Parameter setups
alpha=1/3;
beta=0.99;
sigma=2;
delta=0.025;
a_bar=0;
l_bar=1;
rho=0.5;
sigma_e=0.2;
sigma_z=sigma_e/sqrt(1-rho^2);

% Discretization using Tauchen's method
[lz,lzprob]=TAUCHEN(5,rho,sigma_e,3);
z=exp(lz');
P=lzprob^1000;
N_s=sum(P(1,:).*z);

a_min=a_bar;
a_max=80;
num_a=500;
pas=(a_min+a_max)/num_a;

d=1;
% Discretization of assets

a=a_min:pas:a_max;

% VFI
k_min=20;
k_max=40;

dis=1;

while abs(dis)>=0.01
    if abs(k_max-k_min)<1e-04
        break;
    end
    k_guess=(k_min+k_max)/2;
    r=alpha*(N_s/k_guess)^(1-alpha)+(1-delta);
    w=(1-alpha)*(k_guess/N_s)^alpha;
     % CURRENT RETURN (UTILITY) FUNCTION
    cons = bsxfun(@minus,r*a,a');
    cons = bsxfun(@plus, cons, permute(z*w*l_bar, [1 3 2]));
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
    ret(cons<0)=-Inf;
    
    % INITIAL VALUE FUNCTION GUESS
    v_guess = zeros(numel(z), numel(a));
    
    % VALUE FUNCTION ITERATION
    v_tol = 1;
    while v_tol >1e-06
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        for i=1:numel(z)
        E(:,:,i)=repmat(lzprob(i,:)*v_guess, [numel(a) 1]);
        end
        value_mat=ret+beta*E;
        %value_mat=ret+beta*repmat(permute(lzprob*v_guess,[3 2 1]),[numel(a) 1 1]);
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        [vfn,pol_indx]=max(value_mat,[],2);
        vfn=permute(vfn,[3 1 2]);
        v_tol = abs(max(v_guess(:) - vfn(:)));
        v_guess=vfn;
  
    end
    
    % KEEP DECSISION RULE
    pol_indx = permute(pol_indx,[3 1 2]);
    pol_fn=a(pol_indx);

    % SET UP INITITAL DISTRIBUTION

    Mu=zeros(size(pol_fn));
    Mu(end,end)=1;
    % ITERATE OVER DISTRIBUTIONS
   
    MuNew = zeros(size(Mu));
    mu_tol=1;
    while mu_tol>1e-07
        [emp_ind, a_ind] = find(Mu);
        MuNew = zeros(size(Mu));
    for ii = 1:length(emp_ind)
        apr_ind = pol_indx(emp_ind(ii), a_ind(ii)); 
        MuNew(:, apr_ind) = MuNew(:, apr_ind) + ...
            (lzprob(emp_ind(ii), :)*Mu(emp_ind(ii),a_ind(ii)))';
    end
    mu_tol=max(abs(MuNew(:)-Mu(:)));
    Mu=MuNew;
    end
    s=sum(Mu*a')
    dis=k_guess-s;
    if dis>=0
        k_max=k_guess;
    else
        k_min=k_guess;
    end
end
% Policy functions graph

figure;
plot(a,pol_fn);
title('Policy Function with different Z');
legendcell=cellstr(num2str(z'));
legend(legendcell,'Location','northwest');

% Wealth Distribution
figure;
plot(a,Mu);
title('Wealth Distribution with Different Z');
legend(legendcell,'Location','northwest');

% Lorenz Curve & Gini Coefficient
w=repmat(reshape(a,[numel(a) 1]),numel(z),1);
mu=reshape(Mu,[numel(pol_fn) 1]);
c=cumsum(sortrows([mu mu.*w w],3));
L=bsxfun(@rdivide,c,c(end,:))*100;
wi=repmat(w,1,numel(pol_fn));
wj=repmat(w',numel(pol_fn),1);
g=sum(sum(abs(wi-wj)))/sum(w)/numel(pol_fn)/2;
figure;
area(L(:,1),L(:,2));
axis([0 100 0 100]);
axis('auto y')
hold on;
fplot(@(x) x,':');
title(['Gini Coefficient=',num2str(g)]);

while abs(k_guess-aggsav)>=0.01;
    k_guess=1/2*(k_min+k_max);
    %r_guess=1/2*(r_min+r_max);
    %k_guess=((r_guess-1+delta)/alpha)^(1/(alpha-1))*L;
r=alpha*(k_guess/L)^(alpha-1)+1-delta;
w=(1-alpha)*(k_guess/L)^(alpha);
y_s=Z'*w*l;
    % CURRENT RETURN (UTILITY) FUNCTION
    cons = bsxfun(@minus, r*a', a);
    cons = bsxfun(@plus, cons, permute(y_s, [1 3 2]));
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
    ret(cons<0)=-Inf;
    % INITIAL VALUE FUNCTION GUESS
    v_guess = zeros(nz, num_a);
    % VALUE FUNCTION ITERATION
    v_tol = 1;
    while v_tol >.000001;
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        for i=1:nz;
        E(:,:,i)=repmat(PI(i,:)*v_guess, [num_a 1]);
        end
        v=ret+beta*E;
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        [vmax, pol_indx]=max(v, [], 2);
        vmax=permute(vmax, [3 1 2]);
        v_tol = max(max(abs(vmax-v_guess)));
        %v_tol = abs(max(max(vmax-v_guess)));
        v_guess=vmax;
    end;
     % KEEP DECSISION RULE
    pol_indx=permute(pol_indx, [3 1 2]);
    pol_fn = a(pol_indx);
    
   % SET UP INITITAL DISTRIBUTION
    %Mu=ones(nz, num_a)/(nz*num_a);
    Mu = zeros(nz,num_a);
  Mu(1, 4) = 1; % initial guess: everyone employed, 0 assets
 % Mu = ones(nz, num_a); %alternative initial guess: same mass in all states
%Mu = Mu / sum(Mu(:)); % normalize total mass to 1

% ITERATE OVER DISTRIBUTIONS
% way 1: loop over all non-zeros states
mu_tol = 1; 
while mu_tol > 1e-07
    [emp_ind, a_ind,mass] = find(Mu ); % find non-zero indices
    
    MuNew = zeros(size(Mu));
    for ii = 1:length(emp_ind)
        apr_ind = pol_indx(emp_ind(ii), a_ind(ii)); 
        MuNew(:, apr_ind) = MuNew(:, apr_ind) + ...
            (PI(emp_ind(ii), :)*Mu(emp_ind(ii), a_ind(ii)) )';
    end

    mu_tol = max(max(abs(MuNew-Mu)));
    
    Mu = MuNew ;
end
    aggsav=sum(Mu*a');
    %aggsav=sum(Mu*a');
    if aggsav>=k_guess;
       k_min=k_guess;
    else
       k_max=k_guess;
    end
    [k_guess,aggsav,k_max,k_min,r]
    if abs(k_max-k_min)<10^(-5);
       break;
   end
end
r_cm=1/beta
figure(1)
plot(a,pol_fn)
legend('State 1','State 2','State 3','State 4','State 5','location','northwest')
title(['Policy Function'])
Mu1=Mu';
pop=[Mu1(:,1);Mu1(:,2);Mu1(:,3);Mu1(:,4);Mu1(:,5)];
wealth=repmat(a',[nz 1]);
figure(2)
c_w=gini(pop, wealth,true);
title(['Wealth, gini=',num2str(c_w)])
M=sum(Mu)
figure(3)
plot(a,M)
title(['Distribution of wealth'])
Y=repmat(y_s',[1,num_a]);
A=repmat(a,[nz,1])
c=Y+r*A-pol_fn;
cf=c(:,pol_indx');
cf1=reshape(cf,[nz num_a nz]);
i=1;
while i < nz+1
c1(i,:)=PI(i,:)*cf1(:,:,i);
i=i+1;
end
Eulererror=sum(sum(abs(c.^(-sigma)-beta*c1.^(-sigma)*r).*Mu))