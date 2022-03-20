function [next] = gen_rand(len)
    next(1) = 10;
    for i=2:len
        next(i) = mod( int64( floor( (next(i-1) * 1103515245 + 12345)/65536 ) ), 32768 );
%         next(i) = mod(uint64(next(i)/65536), 32768);
    end
end