pro cute_scatter,nx,ny,area,z,reflection



;assume  z as 1 photons/s/cm2/sr
;slit_area=width*length/3600^2
photons_slit=1.z*area/3282.810874
;asuume 10-10 reduction
counts=photons_slit*reflection
im=make_array(nx,ny,value=counts)
return,im
end