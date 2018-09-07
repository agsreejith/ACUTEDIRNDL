function cute_random,nx,ny,spread,center

array=dblarr(nx,ny)
;column=dblarr(nx)
for i=0,ny-1 do begin
  rand=fix(20*randomu(seed,1))
  ran_val=rand[0]
  for rd=0,ran_val do value=spread*randomn(seed,nx)+center
  array[*,i]=value
endfor
return,array
end