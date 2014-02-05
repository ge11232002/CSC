function PWM = PWMrandomizeBayes(PCM,Dprior,N,W)

% generates N (default 1) random PWM of width drawn from the posterior distribution
% of PWMs. The posterior is propertional to the Dirichlet prior distribution Dprior (which
% might be a mixture) times the mulitnomial likelihood with count matrix PCM.

[A,Win] = size( PCM ) ;

if nargin < 3
    N = 1 ;
end
if nargin < 4
    W = Win ;
end
 
% parameters of prior
alpha0 = Dprior.alpha0 ; % alpha0 is the estimated Dirichlet parameters A x K
pmix = Dprior.pmix ; % pmix mixing proportions 1 x K 

K = length( pmix ) ;

% accumulative distribution for mixing proportions
apmix = pmix ;
for k=2:K % generate cumulative distribution
    apmix(k) = apmix(k) + apmix(k-1) ;
end

for n = 1:N
 
    if W == Win 
    
        for w=1:W
       
            % draw from the mixture component
            k = find( rand < apmix , 1) ;
        
            % draw from component k of dirichlet posterior
            PWM(n,:,w) = dirichlet_sample(alpha0(:,k)+PCM(:,w)) ;
        
        end
        
    else
    
        for w=1:W
            
            % draw from the mixture component
            k = find( rand < apmix , 1) ;
        
            % draw from component k of dirichlet posterior and use random
            % column in PCM
            PWM(n,:,w) = dirichlet_sample(alpha0(:,k)+PCM(:,ceil(Win*rand))) ;
            
        end
        
    end
            
end

if N==1
    PWM = squeeze(PWM)' ;
end


