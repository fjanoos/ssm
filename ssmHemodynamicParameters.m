function [H_estimate] = ssmHemodynamicParameters(data, params, variational)
%____________________________
% ssmHemodynamicParameters.m
%
% function to estimate hemodynamic parameters h
% deal with each dimension of Y independently

% $Id: ssmHemodynamicParameters.m v0.01 2012-06-23 14:20:07 fj $
% $Id: ssmHemodynamicParameters.m v0.02 2012-07-05 12:01:54 fj $ using progress_line
        
    %%%% propaganda
    myLogo					= cafe_logo( mfilename, 'messg', [ 'subject : ' data.subject_id ' * ' datestr( now, 31)]) ;
    H_estimate.logo				= myLogo.tmp ;

% $$$     display('computing hemodynamic parameters ...') ;
    fprintf('\ncomputing hemodynamic parameters ...\n') ;
    
    D						= data.D ;
    L						= data.L ;
    T						= data.T ;
    Y						= data.y ;
    s						= data.s ;
    TR						= data.TR ;
    K						= data.K ;
    tol						= data.tol ;
    
    invSigma_eps				= params.invSigma_eps ;
    Sigma_h					= kron( speye(D), data.Sigma_h ) ;
    mu_h					= kron( ones(D,1), data.mu_h ) ; 
        
    % initialize
    H						= params.H ;

    %%%% PISTI : 2012-04-03
% $$$     wb_h = waitbar(0,'... building Ax=b') ;
    
    % we'll solve Ah = b using conjugate gradients where h is the h[d] vectors stacked from d= 1...L
    A						= zeros( L*D, L*D) ; 
    b						= zeros( L*D, 1) ;
    
% $$$     fprintf('\n* H_estimate - t=')
% $$$     fprintf([ '\n* H_estimate * D=%d * t(%d) = '], D, T) ;
    progress_line( 'init') ;

    for t = 1 : T % --------------------------------------------------------------------------------
	
	%%%% PISTI : 2012-04-03
% $$$         waitbar(t/T) ;
% $$$ 	fprintf('%d',t)
	progress_line( 1, T, t) ;
	fprintf(' hrf') ;
        [nu_t, Lambda_t]			= buildMomentMatrices( variational, t, data) ;
% $$$ 	fprintf('e')
	
        for d = 1 : D             % the column index of A            

	    %%%% PISTI : 2012-04-03
% $$$             waitbar(((t-1)*D+d)/D/T) ;
% $$$ 	    fprintf('\n\td=%d ',d)
% $$$ 	    fprintf('d=%d ',d)
% $$$ 	    fprintf('%d ',d) ; pause( .01)
	    
            d_idx				= [ (d-1) * L+1 : d * L ] ;
            for dp = 1 : D
		
		%%%% PISTI : 2012-04-03
% $$$                 waitbar((((t-1)*D+d-1)*D+dp)/D/D/T) ;
% $$$ 		fprintf('.')
		
		dp_idx				= [(dp-1) * L+1 : dp * L ] ;
                %accumulate \sum_t \Sigma_\eps\inv[d,d']\Lambda_t[d,d']
                A(d_idx, dp_idx)		= A( d_idx, dp_idx) + invSigma_eps( dp, d) * Lambda_t{ d, dp} ;
                b(d_idx)			= b( d_idx) + invSigma_eps( dp, d) * Y( t, dp) * nu_t{ d} ;
            end %dp
        end % d    
    end % t
    
    %%%% PISTI : 2012-04-03
% $$$     close(wb_h)  ;
% $$$     fprintf('\n\n')
    progress_line( 'term') ;
    
    % add in the 't' independent terms
    A						= Sigma_h * A + speye( L* D) ; % LxL identity matrices along the diagonal
    b						= Sigma_h * b + mu_h ;
    
    % we need to make the thing sparse by ensuring that invSigma_eps(dp,d) is sparse 
    
% $$$     display( '... solving Ax=b') ;
    fprintf( '\n... solving Ax=b') ;
    [ h, flag, relres, iter, resvec]		= pcg( A, b, 1e-3, 50)  ;
    display ( [ ' : d = ',			num2str( d), ...
                ' , cdg flag = ',		num2str( flag),...
                ' , relative resdiual norm = ',	num2str( relres), ...
                ' , iterations = ',		num2str( iter)]) ;
    % reorder h as the H matrices
    for d = 1:D
        H_estimate.H(d,d,:)			= h( (d-1)*L+1 : d*L )  ;    
    end %
   
    H_estimate.rel_change			= norm( H_estimate.H(:) - params.H(:))/norm( params.H(:)) ;
    
% $$$     display('computing hemodynamic parameters ... DONE') ;         
    fprintf('\ncomputing hemodynamic parameters ... DONE\n') ;

end


%%%% function --------------------------------------------------------------------------------------
function [ nu_t, Lambda_t] = buildMomentMatrices( variational, t, data)
% construct the nu_t vector and Lambda_t matrix, for all 'd'

    D						= data.D ;
    L						= data.L ;
    K						= data.K ;
    
    nu_t					= cell( D) ;
    Lambda_t					= cell( D, D) ;
    for d = 1 : D
        nu_t{d}					= sparse( L , 1) ;
        for dp = 1 : D
            Lambda_t{d,dp}			= sparse( L , L) ;
        end
    end

    for l_idx = 1 : L   %indexing variable
	fprintf('.')
        l					= l_idx - 1 ;  %0-indexed value 
        if t-l < 1
            continue ;
        end
        [ nu_qt_l, Lambda_qt_l]			= variationalMomentsOfZ( variational, t-l, K) ;
        for d = 1 : D
            nu_t{d}(l_idx)			= nu_qt_l( d) ;
            for dp = 1 : D
                Lambda_t{d,dp}(l_idx,l_idx)	=  Lambda_qt_l( d, dp) ;   
            end
        end
        for m_idx = 1 : L   %indexing variable
            m					= m_idx - 1 ;  %0-indexed value
            if t-m < 1 || m == l
                continue ;
            end 
            nu_qt_m				= variationalMomentsOfZ( variational, t-m, K) ;
            for d = 1 : D
                for dp = 1 : D
                    Lambda_t{d,dp}(l_idx,m_idx)	= nu_qt_l( d) * nu_qt_m( dp) ;
                end
            end
        end
    end %l_idx
end
