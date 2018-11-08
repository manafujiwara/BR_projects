function v = boxcar(data,s)

ss = s/2;
nSamp = size(data,2);

for isamp =1 :nSamp
    if 1 <= isamp && isamp <= ss-1
        %         y = sum( data( :, 1 : isamp + ss ), 2 )/isamp+ss ;
        y = nanmean(data( :, 1 : isamp + ss ), 2) ;
        
    elseif ss <= isamp && isamp <= nSamp -ss-1
        %         y = sum( data( :, isamp-ss+1 : isamp+ss ), 2 )/s ;
        y = nanmean( data( :, isamp-ss+1 : isamp+ss ), 2 ) ;
        
    elseif nSamp-ss <= isamp && isamp<= nSamp
        %         y = sum( data( :, isamp-ss+1 : nSamp), 2 )/(nSamp-isamp+ss);
        y = nanmean( data( :, isamp-ss+1 : nSamp ), 2 );
        
    end
    v(:, isamp) = y;
    
    if sum(isinf(y))~=0
        
    end
    
%     if mod(isamp,1000)==0
%         disp([ 'Boxcar: ' num2str(isamp) '/' num2str(nSamp) ' samples'])
%     end
    
end

