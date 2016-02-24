function [noise] = ssmEstimateNoiseParameter( data, params, variational)
%_____________________________
% ssmEstimateNoiseParameter.m
%
% estimates the Sigma_eps and inverse Sigma_eps parameters
    
% $Id: ssmEstimateNoiseParameter.m v0.01 2012-06-23 14:34:01 fj $
% $Id: ssmEstimateNoiseParameter.m v0.02 2012-07-05 11:14:07 fj $ using progress_line
    
%%%% propaganda
myLogo						= cafe_logo( mfilename, 'messg', [ 'subject : ' data.subject_id ' * ' datestr( now, 31)]) ;
noise.logo					= myLogo.tmp ;
    
% $$$     display( 'computing noise parameters ...') ;
    fprintf( '\ncomputing noise parameters ...\n') ;
    
    D						= data.D ;
    L						= data.L ;
    T						= data.T ;
    Y						= data.y ;
    s						= data.s ;
    TR						= data.TR ;
    K						= data.K ;
    tol						= data.tol ;
    H						= params.H ;
    
    Sigma_eps					= zeros(D,D) ;
    
    %%%% PISTI : 2012-04-03
% $$$     h = waitbar(0, 'Estimating noise', 'Name', 'ssmEstimateNoiseParameter') ;
% $$$     fprintf('\n* noise - t=') ;
% $$$     fprintf([ '\n* noise * L=%d * t(%d) = '], L, T) ;
    progress_line( 'init') ;
    
    for t = 1 : T
	
	%%%% PISTI : 2012-04-03
% $$$         waitbar( t/T, h) ;
% $$$ 	fprintf( '%d ',t) ;
	progress_line( 1, T, t) ;
	
        y_t					= Y(t,:)' ;
        Sigma_eps				= Sigma_eps + y_t*y_t' ;        
	
        for l_idx = 1 : L   %indexing variable
	    
% $$$ 	    fprintf('.') ;
            l					= l_idx - 1 ;  %0-indexed value 
            if t-l < 1
                continue ;
            end
            [ nu_qt_l, Lambda_qt_l]		= variationalMomentsOfZ( variational, t-l, K) ;
            tmp					= y_t*(H(:,:,l_idx)*nu_qt_l)' ;
            Sigma_eps				= Sigma_eps - tmp - tmp' + H(:,:,l_idx)*Lambda_qt_l*H(:,:,l_idx) ;
	    
            for m_idx = 1 : L   %indexing variable
                m				= m_idx - 1 ;  %0-indexed value
                if t-m < 1 || m == l
                    continue ;
                end 
                [ nu_qt_m, Lambda_qt_m]		= variationalMomentsOfZ( variational, t-m, K) ;
                tmp				= H(:,:,l_idx)*nu_qt_l*nu_qt_m'*H(:,:,l_idx) ;
                Sigma_eps			= Sigma_eps + tmp + tmp' ;
            end %m            
        end %l
    end %t
    
    %%%% PISTI : 2012-04-03
% $$$     close (h) ;
% $$$     fprintf('\n\n') ;
    progress_line( 'term') ;

    Sigma_eps					= Sigma_eps / T ;
    
    if Sigma_eps ~= Sigma_eps'
        error( 'ssmEstimateNoiseParameter : Sigma_eps not symmetric') ;
    end
    
    noise.invSigma_eps				= CholeskyInverse(Sigma_eps) ;
    noise.invSigma_eps				= Sparsify(noise.invSigma_eps) ;
    noise.Sigma_eps				= CholeskyInverse(noise.invSigma_eps) ;

    noise.rel_change.Sigma_eps			= norm(noise.Sigma_eps(:) - params.Sigma_eps(:)) / norm(params.Sigma_eps(:)) ;
    noise.rel_change.invSigma_eps		= norm(noise.invSigma_eps(:) - params.invSigma_eps(:))/norm(params.invSigma_eps(:)) ;
    
% $$$ end
