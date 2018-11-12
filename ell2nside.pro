FUNCTION ELL2NSIDE, l
result=2^(ceil(alog(l)/alog(2))-1) > 16
return,result
end
