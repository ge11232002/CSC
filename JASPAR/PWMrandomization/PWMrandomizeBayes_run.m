function PWMrandomizeBayes(infile,nmat,width)
% this script makes an toy example call of PWMrandomizeBayes 

% motif model
motifmodel = 'transfac' ;% 'jaspar'; % default transfac
switch motifmodel
    case 'jaspar'
        load Dirichlet_mixture_EM_estimation_jaspar_all_K6_ll-70020.6.mat  
    case 'transfac'
        load Dirichlet_mixture_EM_estimation_all_K6_ll-459717.mat 
    case 'ones'
        alpha1 = ones( A, W(1) );
        pmix = 1 ;
    case 'countmatrix'
        % user defined count matrix for biased searches
end
% end motifs model

Dprior.alpha0 = alpha1 ;
Dprior.pmix = pmix ;

PCM = load(infile,'ASCII');

%PCM = [ [13 1 1 1]' [ 4 4 4 4]' [16 0 0 0]' [0 0 8 8]' [1 0 3 12]' ] ; 

%nmat = 2;
%width = 10;
for n = 1:1

    PWMsample = PWMrandomizeBayes(PCM,Dprior,nmat,width) 

    for matno = 1:nmat
    	fd = fopen(strcat('matrix',num2str(matno),'.txt'),'w');
    	for i=1:4
		for j=1:width
		    fprintf(fd,'3.2%f\t',10*PWMsample(matno,i,j));
		end
                fprintf(fd,'\n');
    	end	
	fclose(fd);
    end	
    %squeeze(PWMsample)'
end



