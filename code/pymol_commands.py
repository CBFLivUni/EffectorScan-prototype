hide lines, all 
show cartoon

set_color n0, [0.051, 0.341, 0.827]
set_color n1, [0.416, 0.796, 0.945]
set_color n2, [0.996, 0.851, 0.212]
set_color n3, [0.992, 0.490, 0.302]
color n0, b < 100; color n1, b < 90
color n2, b < 70;  color n3, b < 50

set ray_trace_mode, 1; bg_color white; set antialias,3
set ray_trace_gain, 0.001 
set cartoon_fancy_helices, 1

ray 900, 900
png FOXG_07142.png, dpi=300
#png 1sj1.png, dpi=1000, ray=1