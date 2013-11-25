function vScales = computescale ( cellDFB, dRatio, nStart, nEnd, coefMode )
% COMPUTESCALE   Comupute display scale for PDFB coefficients
%
%       computescale(cellDFB, [dRatio, nStart, nEnd, coefMode])
%
% Input:
%	cellDFB:	a cell vector, one for each layer of 
%		subband images from DFB.
%   dRatio:
%       display ratio. It ranges from 1.2 to 10.
%
%   nStart:
%       starting index of the cell vector cellDFB for the computation. 
%       Its default value is 1.
%   nEnd:
%       ending index of the cell vector cellDFB for the computation. 
%       Its default value is the length of cellDFB.
%   coefMode: 
%       coefficients mode (a string): 
%           'real' ----  Highpass filters use the real coefficients. 
%           'abs' ------ Highpass filters use the absolute coefficients. 
%                        It's the default value
% Output:
%	vScales ---- 1 X 2 vectors for two scales.
%
% History: 
%   10/03/2003  Creation.
%   04/01/2004  Limit the display scale into the range of 
%               [min(celldfb), max(celldfb)] or [min(abs(celldfb)), max(abs(celldfb))]
%
% See also:     SHOWPDFB

if ~iscell(cellDFB)
    error ('Error in computescale.m! The first input must be a cell vector!');    
end

% Display ratio
if ~exist('dRatio', 'var')
    dRatio = 2 ;
elseif dRatio < 1
    display ('Warning! the display ratio must be larger than 1!Its defualt value is 2!');
end

% Starting index for the cell vector cellDFB
if ~exist('nStart', 'var')
    nStart = 1 ;
elseif nStart < 1 | nStart > length(cellDFB)
    display ('Warning! The starting index from 1 to length(cellDFB)! Its defualt value is 1!');
    nStart = 1 ;
end

% Starting index for the cell vector cellDFB
if ~exist('nEnd', 'var')
    nEnd = length(cellDFB) ;
elseif nEnd < 1 | nEnd > length(cellDFB)
    display ('Warning! The ending index from 1 to length(cellDFB)! Its defualt value is length(cellDFB)!');
    nEnd = length( cellDFB ) ;
end

% Coefficient mode
if ~exist('coefMode', 'var')
    coefMode = 'abs' ;
elseif ~strcmp(coefMode,'real') & ~strcmp(coefMode, 'abs')
    display ('Warning! There are only two coefficients mode: real, abs! Its defualt value is "abs"!');
    coefMode = 'abs' ;
end

% Initialization
dSum = 0 ;
dMean = 0 ;
% Added on 04/01/04 by jpzhou
dMin = 1.0e14 ;
dMax = -1.0e14 ;
dAbsMin = 1.0e14 ;
dAbsMax = -1.0e14 ;
dAbsSum = 0 ;
nCount = 0 ;
vScales = zeros (1, 2);

if strcmp( coefMode, 'real') %Use the real coefficients
    % Compute the mean of all coefficients
    for i= nStart : nEnd
        if iscell( cellDFB{i} ) % Check whether it is a cell vector
            m = length(cellDFB{i});
            for j=1 : m
                
                % Added on 04/01/04 by jpzhou
                dUnitMin = min( min( cellDFB{i}{j} ) ) ;
                if dUnitMin < dMin
                    dMin = dUnitMin ;
                end
                dUnitMax = max( max( cellDFB{i}{j} ) ) ;
                if dUnitMax > dMax
                    dMax = dUnitMax ;
                end
                
                dSum = dSum + sum( sum( cellDFB{i}{j} ));
                nCount = nCount + prod( size( cellDFB{i}{j} ) );
            end
        else
            
            % Added on 04/01/04 by jpzhou
            dUnitMin = min( min( cellDFB{i} ) ) ;
            if dUnitMin < dMin
                dMin = dUnitMin ;
            end
            dUnitMax = max( max( cellDFB{i} ) ) ;
            if dUnitMax > dMax
                dMax = dUnitMax ;
            end
                
            dSum = dSum + sum( sum( cellDFB{i} ));
            nCount = nCount + prod( size( cellDFB{i} ) ); 
        end
    end
	if nCount < 2 | abs(dSum) < 1e-10
        error('Error in computescale.m! No data in this unit!' );
	else
        dMean = dSum / nCount ;
    end
    
    % Compute the STD.
    dSum = 0 ;
    for i= nStart : nEnd
        if iscell( cellDFB{i} ) %Check whether it is a cell vector
            m = length(cellDFB{i});
            for j=1 : m                   
                dSum = dSum + sum( sum( (cellDFB{i}{j}-dMean).^2 ));
                %nCount = nCount + prod( size( cellDFB{i}{j} ) );
            end
        else
            dSum = dSum + sum( sum( (cellDFB{i}-dMean).^2 ));
            %nCount = nCount + prod( size( cellDFB{i} ) ); 
        end
    end
	dStd = sqrt( dSum / (nCount-1) );
    
    % Modified on 04/01/04
    %dMin = -1.0e10 ;
    %dMax = 1.0e10 ;
    vScales( 1 ) = max( dMean - dRatio * dStd, dMin ) ;
    vScales( 2 ) = min( dMean + dRatio * dStd, dMax ) ;
    
else %Use the absolute coefficients
    % Compute the mean of absolute values
    for i= nStart : nEnd
        if iscell( cellDFB{i} ) % Check whether it is a cell vector
            m = length(cellDFB{i});
            for j=1 : m
                
                % Added on 04/01/04 by jpzhou
                dUnitMin = min( min( abs(cellDFB{i}{j}) ) ) ;
                if dUnitMin < dAbsMin
                    dAbsMin = dUnitMin ;
                end
                dUnitMax = max( max( abs(cellDFB{i}{j}) ) ) ;
                if dUnitMax > dAbsMax
                    dAbsMax = dUnitMax ;
                end
                    
                dAbsSum = dAbsSum + sum( sum( abs(cellDFB{i}{j}) ));
                nCount = nCount + prod( size( cellDFB{i}{j} ) );
            end
        else
            
            % Added on 04/01/04 by jpzhou
            dUnitMin = min( min( abs(cellDFB{i}) ) ) ;
            if dUnitMin < dAbsMin
                dAbsMin = dUnitMin ;
            end
            dUnitMax = max( max( abs(cellDFB{i}) ) ) ;
            if dUnitMax > dAbsMax
                dAbsMax = dUnitMax ;
            end
            
            dAbsSum = dAbsSum + sum( sum( abs(cellDFB{i}) ));
            nCount = nCount + prod( size( cellDFB{i} ) ); 
        end
    end
	if nCount < 2 | dAbsSum < 1e-10
        error('Error in computescale! No data in this unit!');
	else
        dAbsMean = dAbsSum / nCount ;
    end
    
    % Compute the std of absolute values
    dSum = 0 ;
    for i= nStart : nEnd
        if iscell( cellDFB{i} ) %Check whether it is a cell vector
            m = length( cellDFB{i} );
            for j = 1 : m
                dSum = dSum + sum( sum( (abs(cellDFB{i}{j})-dAbsMean).^2 ));
            end
        else
            dSum = dSum + sum( sum( (abs(cellDFB{i})-dAbsMean).^2 ));  
        end
    end
    dStd = sqrt( dSum / (nCount-1) ); 
    
    % Modified on 04/01/04   
    % Compute the scale values
    vScales( 1 ) = max( dAbsMean - dRatio * dStd, dAbsMin ) ;
    vScales( 2 ) = min( dAbsMean + dRatio * dStd, dAbsMax ) ;	
end
