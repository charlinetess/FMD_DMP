# This file loads the date produced from the file DMPCoordSearchPref and generates a plot of the coordinate estimates x and y 

# 	
# 	
# 	`7MMF'        .g8""8q.      db      `7MM"""Yb.       `7MM"""Yb.      db   MMP""MM""YMM   db
# 	  MM        .dP'    `YM.   ;MM:       MM    `Yb.       MM    `Yb.   ;MM:  P'   MM   `7  ;MM:
# 	  MM        dM'      `MM  ,V^MM.      MM     `Mb       MM     `Mb  ,V^MM.      MM      ,V^MM.
# 	  MM        MM        MM ,M  `MM      MM      MM       MM      MM ,M  `MM      MM     ,M  `MM
# 	  MM      , MM.      ,MP AbmmmqMA     MM     ,MP       MM     ,MP AbmmmqMA     MM     AbmmmqMA
# 	  MM     ,M `Mb.    ,dP'A'     VML    MM    ,dP'       MM    ,dP'A'     VML    MM    A'     VML
# 	.JMMmmmmMMM   `"bmmd"'.AMA.   .AMMA..JMMmmmdP'       .JMMmmmdP'.AMA.   .AMMA..JMML..AMA.   .AMMA.
# 	
# 	
rats=load("experimentDMP.jld2"); # load documents 
include("DiverseFunctions.jl") # load fucntions 

parameters=rats["parameters"]
featuresexperiment=rats["features"]; # [numberofrats, numberofdays, numberoftrials,randomstartingpos,centres];
data=rats["data"]


# 	
# 	              ,,                                                  ,,                    ,,    ,,
# 	  .g8"""bgd `7MM                                                  db                  `7MM    db
# 	.dP'     `M   MM                                                                        MM
# 	dM'       `   MMpMMMb.  ,pW"Wq.   ,pW"Wq.  ,pP"Ybd  .gP"Ya      `7MM  `7MMpMMMb.   ,M""bMM  `7MM  ,p6"bo   .gP"Ya  ,pP"Ybd
# 	MM            MM    MM 6W'   `Wb 6W'   `Wb 8I   `" ,M'   Yb       MM    MM    MM ,AP    MM    MM 6M'  OO  ,M'   Yb 8I   `"
# 	MM.           MM    MM 8M     M8 8M     M8 `YMMMa. 8M""""""       MM    MM    MM 8MI    MM    MM 8M       8M"""""" `YMMMa.
# 	`Mb.     ,'   MM    MM YA.   ,A9 YA.   ,A9 L.   I8 YM.    ,       MM    MM    MM `Mb    MM    MM YM.    , YM.    , L.   I8
# 	  `"bmmmd'  .JMML  JMML.`Ybmd9'   `Ybmd9'  M9mmmP'  `Mbmmd'     .JMML..JMML  JMML.`Wbmd"MML..JMML.YMbmd'   `Mbmmd' M9mmmP'
# 	
# 	

# chose rat :
indexrat=2;
# Chose day :
indexday1=1;
indexday2=2;

indextrial1=2;
indextrial2=4;


# 	
# 	                                                                                                    ,,
# 	  .g8"""bgd                                                mm          `7MM"""YMM            mm     db                             mm
# 	.dP'     `M                                                MM            MM    `7            MM                                    MM
# 	dM'       ` ,pW"Wq.`7MMpMMMb.pMMMb. `7MMpdMAo.`7MM  `7MM mmMMmm .gP"Ya   MM   d    ,pP"Ybd mmMMmm `7MM  `7MMpMMMb.pMMMb.   ,6"Yb.mmMMmm .gP"Ya
# 	MM         6W'   `Wb MM    MM    MM   MM   `Wb  MM    MM   MM  ,M'   Yb  MMmmMM    8I   `"   MM     MM    MM    MM    MM  8)   MM  MM  ,M'   Yb
# 	MM.        8M     M8 MM    MM    MM   MM    M8  MM    MM   MM  8M""""""  MM   Y  , `YMMMa.   MM     MM    MM    MM    MM   ,pm9MM  MM  8M""""""
# 	`Mb.     ,'YA.   ,A9 MM    MM    MM   MM   ,AP  MM    MM   MM  YM.    ,  MM     ,M L.   I8   MM     MM    MM    MM    MM  8M   MM  MM  YM.    ,
# 	  `"bmmmd'  `Ybmd9'.JMML  JMML  JMML. MMbmmd'   `Mbod"YML. `Mbmo`Mbmmd'.JMMmmmmMMM M9mmmP'   `Mbmo.JMML..JMML  JMML  JMML.`Moo9^Yo.`Mbmo`Mbmmd'
# 	                                      MM
# 	                                    .JMML.


# establish the grid of points in the pool
a=[-parameters[:R]+2(k-1) for k=1:floor(2*parameters[:R]/2)];

# initalize the xestimate variable
global xestimatetrial1
xestimatetrial1 = zeros(length(a),length(a));
global yestimatetrial1 
yestimatetrial1 = zeros(length(a),length(a));
global xestimatetrial2
xestimatetrial2= zeros(length(a),length(a));
global yestimatetrial2
yestimatetrial2= zeros(length(a),length(a));

global xweighttrialbegin
xweighttrialbegin=data[indexrat][indexday1].day[indextrial1].Xcoord;
global yweighttrialbegin
yweighttrialbegin=data[indexrat][indexday1].day[indextrial1].Ycoord;
global xweighttrialend
xweighttrialend=data[indexrat][indexday2].day[indextrial2].Xcoord;
global yweighttrialend
yweighttrialend=data[indexrat][indexday2].day[indextrial2].Ycoord;

# for each place point in the grid, calculate the x and y values 
for i = 1:length(a)
    for j = 1:length(a)

        # make sure the point is in the pool
        if sqrt((a[i]^2+a[j]^2)) < parameters[:R]
        
            # determine the place cell activity at this point
            F = placecells([a[i],a[j]],parameters[:centres],parameters[:σPC])           
            # determine the actor activity
            global xestimatetrial1
            xestimatetrial1[i,j] = dot(xweighttrialbegin,F);
            global xestimatetrial2
            xestimatetrial2[i,j] = dot(xweighttrialend,F);
            global yestimatetrial1
            yestimatetrial1[i,j] = dot(yweighttrialbegin,F);
            global yestimatetrial2
            yestimatetrial2[i,j] = dot(yweighttrialend,F);
                
        else
        	global xestimatetrial1
            xestimatetrial1[i,j] = NaN;
            global xestimatetrial2
            xestimatetrial2[i,j] = NaN;
            global yestimatetrial1
            yestimatetrial1[i,j] = NaN;
            global yestimatetrial2
            yestimatetrial2[i,j] = NaN;
        end
    end
end

# 	
# 	                                                                   ,,
# 	  .g8"""bgd                          mm               `7MM"""Mq.   db
# 	.dP'     `M                          MM                 MM   `MM.
# 	dM'       ``7Mb,od8 .gP"Ya   ,6"Yb.mmMMmm .gP"Ya        MM   ,M9 `7MM  ,p6"bo
# 	MM           MM' "',M'   Yb 8)   MM  MM  ,M'   Yb       MMmmdM9    MM 6M'  OO
# 	MM.          MM    8M""""""  ,pm9MM  MM  8M""""""       MM         MM 8M
# 	`Mb.     ,'  MM    YM.    , 8M   MM  MM  YM.    ,       MM         MM YM.    ,
# 	  `"bmmmd' .JMML.   `Mbmmd' `Moo9^Yo.`Mbmo`Mbmmd'     .JMML.     .JMML.YMbmd'
# 	
# 	
using PyPlot

fig=figure(figsize=(7,7))
theta=0:pi/50:2*pi; # to draw circle 
 

# plot heatmap 

# plot estimated x coordinate for the first day chosen  :
ax1 = subplot2grid((2, 2), (0, 0),title="X estimate before learning ")

s1 = pcolormesh(a,a,xestimatetrial1)#,facecolors=get_cmap(homemadecoloragain).o(xestimatetrial1./maximum(xestimatetrial1[findall(x->!==(x, NaN), xestimatetrial1)])), linewidth=0, antialiased=true, shade=false,alpha=1 ,rstride=1, cstride=1,zorder=3)
colorbar()
plot(parameters[:R]*cos.(theta),parameters[:R]*sin.(theta),ls="--",color=[169/255,169/255,169/255],zorder=1)

ax=gca() 
ax[:set_axis_off]()
gca()[:grid](false);

ax2 = subplot2grid((2,2), (0,1),title="X estimate after learning ")
# plt heatmap 
s2 = pcolormesh(a,a,xestimatetrial2)
colorbar()
plot(parameters[:R]*cos.(theta),parameters[:R]*sin.(theta),ls="--",color=[169/255,169/255,169/255],zorder=1)

ax=gca() 
ax[:set_axis_off]()
gca()[:grid](false);

# plot heatmap
ax3 = subplot2grid((2,2), (1,0),title="Y estimate before learning ")

s3 = pcolormesh(a,a,yestimatetrial1)#,facecolors=get_cmap(homemadecoloragain).o(yestimatetrial1./maximum(yestimatetrial1[findall(x->!==(x, NaN), yestimatetrial1)])), linewidth=0, antialiased=true, shade=false,alpha=1 ,rstride=1, cstride=1,zorder=3)
colorbar()
plot(parameters[:R]*cos.(theta),parameters[:R]*sin.(theta),ls="--",color=[169/255,169/255,169/255],zorder=1)

ax=gca() 
ax[:set_axis_off]()
gca()[:grid](false);

# plot heatmap
ax4 = subplot2grid((2,2), (1,1),title="Y estimate at after learning ")

s4 = pcolormesh(a,a,yestimatetrial2)
colorbar()
plot(parameters[:R]*cos.(theta),parameters[:R]*sin.(theta),ls="--",color=[169/255,169/255,169/255],zorder=1)

ax4=gca() 
ax4[:set_axis_off]()
gca()[:grid](false);

# 	
# 	                                                     ,,
# 	 .M"""bgd                               `7MM"""Mq.   db
# 	,MI    "Y                                 MM   `MM.
# 	`MMb.      ,6"Yb.`7M'   `MF'.gP"Ya        MM   ,M9 `7MM  ,p6"bo
# 	  `YMMNq. 8)   MM  VA   ,V ,M'   Yb       MMmmdM9    MM 6M'  OO
# 	.     `MM  ,pm9MM   VA ,V  8M""""""       MM         MM 8M
# 	Mb     dM 8M   MM    VVV   YM.    ,       MM         MM YM.    ,
# 	P"Ybmmd"  `Moo9^Yo.   W     `Mbmmd'     .JMML.     .JMML.YMbmd'
# 	
# 	
savefig("Coordinates_days$(indexday1)_$(indexday2)_trials$(indextrial1)_$(indextrial2).png")

show()
