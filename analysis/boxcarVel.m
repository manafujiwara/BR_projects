function smoothVel = boxcarVel(data,s,p)
%%
% s: n sample in moving window
% p: 

%%
for is = 1:length(s)
    ss = s(is);
    for ip = 1:length(p)
        pp = p(ip);
        
        for isamp = ss+pp+1 : length(data)
            
            lpres = sum( data( isamp-ss : isamp, : ), 1 )/(ss+1) ;
            lpast = sum( data( (isamp-pp)-ss : isamp-pp  , : ) , 1)/(ss +1);
            
            v(isamp,:) = (lpres - lpast)/(pp+1);
            smoothVel(is,ip,isamp,:) = squeeze(v(isamp,:));
            
        end
        
    end
end

end
