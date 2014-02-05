function [alpha0,pmix,ll] = Dirichlet_mixture_EM_estimation(c,K,alpha0,pmix)
% [alpha0,pmix] = Dirichlet_mixture_EM_estimation(c,K)
% returns the maximum Liklihood estimate of the parameters of
% Dirichlet mixture model.
% c is the samples summarized as counts for each of A letters, N x A
% K is the number of sought component.
% alpha0 is the estimated Dirichlet parameters A x K
% pmix mixing proportions 1 x K 
% ll log likelihood

% Ole Winther, binf.ku.dk, June 2006

[ N , A ] = size( c ) ;
C = sum( c , 2 ) ; % N x 1  
oN = ones( N , 1 ) ;
oK = ones( K , 1 ) ;
oA = ones( A , 1 ) ;

gam = zeros( N , K ) ;  
contrib_n = zeros(A,K) ;
contrib_d = zeros(1,K) ;

% random initialization
if nargin < 3 
    epsilon = 0.1 ;
    alpha0 = epsilon * ( 2 * rand( A , K ) - 1 ) + 2; % A x K
end
Alpha0 = sum( alpha0 ) ; % 1 x K
if nargin < 4
    pmix = ones( 1 , K ) / K ; % 1 x K 
end

% minimum alpha0 value
alpha0_min = 1e-8;
ftol = 1e-10 ;
iteouter_max = 10000 ;
iteinner_max = 1;
dll_min = 1e-6 ;

ite=0; dll = Inf; 
while ite < iteouter_max && ( dll > dll_min || ite == 1 )  
    
    ite = ite + 1 ;
    
    % E-step
    contrib = log( pmix ) + gammaln( Alpha0 ) - sum( gammaln( alpha0 ) ) ; 
    for k=1:K
        alpha0_k = alpha0(:,k) ; % A x 1
        gam(:,k) = sum( gammaln( c + alpha0_k( : , oN )' ) , 2 ) ...
            - gammaln( C + Alpha0(oN,k) ) ; % alpha0sum( c + alpha0_k( : , oN )' , 2 ) ) ; % N x K
    end
    gam = gam + contrib( oN , : )  ; % N x K 
    maxgam = max( gam , [] , 2 ) ; % N x 1 
    gam = exp( gam - maxgam( : , oK ) ) ;
    
    ll(ite) = sum( maxgam + log( sum( gam , 2 ) ) ) ;
    
    dll = ll(ite)-ll(max(ite-1,1)) ;
    
    sumgam = 1 ./ sum( gam , 2 ) ; 
    gam = gam .* sumgam( : , oK ) ; % responsibility N x K
    
    % M-step
    pmix = sum( gam ) / N ;  % 1 x K 
   
    dalpha0 = Inf ; iteinner = 0 ;
    while sum(sum(abs(dalpha0))) > ftol & iteinner < iteinner_max
        iteinner = iteinner + 1 ;
        for k=1:K
            alpha0_k = alpha0(:,k) ;
            contrib_n(:,k) = psi( c' + alpha0_k( : , oN ) ) * gam(:,k) ; % A x 1
            contrib_d(k) = psi( C' + Alpha0( k ) ) * gam(:,k) - pmix(k) * N * psi( Alpha0(k) ) ; %
        end
        contrib_n = contrib_n - N * psi( alpha0 ) .* pmix( oA , : ) ;
        %sum(sum(abs(dalpha0)))
       
    
        % ( contrib_n ./ contrib_d( oA , : ) )
        
        dalpha0 = max( alpha0 .* ( contrib_n ./ contrib_d( oA , : ) ) , alpha0_min ) - alpha0 ;
        alpha0 = alpha0 + dalpha0 ;
        
        Alpha0 = sum( alpha0 ) ;
        
    end
    if ~mod(ite,iteouter_max/10)
        dll
        sum(sum(abs(dalpha0)))
        pmix
        alpha0
    end
    
end
