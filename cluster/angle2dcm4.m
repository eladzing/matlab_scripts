function mat = angle2dcm4(a1,a2,a3)
mat = angle2dcm(a1,a2,a3);
mat(4,4)=1;