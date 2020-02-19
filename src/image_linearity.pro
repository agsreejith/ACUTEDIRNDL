pro image_linearity,in_image,out_image,line_file
if not keyword_defined(line_file) then line_file = 0

if line_file eq 0 then begin
   limit=60000
;  Linear model Poly3:
;  f(x) = p1*x^3 + p2*x^2 + p3*x + p4
;  Coefficients (with 95% confidence bounds):
;  p1 =      -26.67  (-154.7, 101.4)
;  p2 =       85.71  (-337, 508.4)
;  p3 =       -90.5  (-554.5, 373.5)
;  p4 =       32.46  (-137, 201.9)
;
;  Goodness of fit:
;  SSE: 2.286e-05
;  R-square: 0.9992
;  Adjusted R-square: 0.9967
;  RMSE: 0.004781
 

  p1 =  -26.67  
  p2 =   85.71 
  p3 =  -90.5  
  p4 =   32.46

  s=size(in_image)

  row=s[2]
  column=s[1]
  for i=0,column-1 do begin
    for j=0,row-1 do begin
      value=in_image[i,j]
      if value gt 60000 then begin
        if value gt 72000 then y = 72000 else begin
          x=value/60000
          y=p1*x^3 + p2*x^2 + p3*x + p4
          y=y*60000
        endelse  
      endif else begin
        y=value
      endelse
  out_image[i,j]=y
  endfor
endfor
endif else begin
     readcol,line_file,input,output,'D,D'
     s=size(in_image)
     row=s[2]
     column=s[1]
     for i=0,column-1 do begin
        for j=0,row-1 do begin
           value=in_image[i,j]
           y=interpol(output,input,value)
           out_image[i,j]=y
      endfor
    endfor
endelse
end
