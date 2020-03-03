# Codigos para gerar imagem 
set terminal png font 'Helvetica' 20.5 size 800,600

#set title 'FHN DDM Droplet 7.5um' 
#set title 'DDM ' 
#set title 'QCM ' 
set title 'DM ' 

unset key 
set tic scale 0
#set palette rgbformula -7,2,-7 
set cbrange [0:1.2] 
#set cbrange [0:0.12] 
set cblabel 'Concentration u '
#set cblabel 'Concentration v '
set xlabel 'X'
set ylabel 'Y'
set xrange [0:50] 
set yrange [0:50] 

#do for [i=0:399]{
#do for [i=0:10000:1000]{
#do for [i=0:50000:1000]{
#do for [i=0:200000:1000]{
#do for [i=0:300000:1000]{
do for [i=0:1000000:1000]{


	outfile = sprintf('Droplet_2d_%d.png',i)
	#outfile = sprintf('Droplet_2d_%d_VV.png',i)
	
	set output outfile
	
 
    plot 'UU_'.i.'.dat' matrix with image 
    #plot 'VV_'.i.'.dat' matrix with image 
    
}

