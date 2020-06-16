# 	
# 	                                     ,,             ,,
# 	`7MMF'                             `7MM           `7MM           mm
# 	  MM                                 MM             MM           MM
# 	  MM         ,pW"Wq.   ,6"Yb.   ,M""bMM        ,M""bMM   ,6"Yb.mmMMmm  ,6"Yb.
# 	  MM        6W'   `Wb 8)   MM ,AP    MM      ,AP    MM  8)   MM  MM   8)   MM
# 	  MM      , 8M     M8  ,pm9MM 8MI    MM      8MI    MM   ,pm9MM  MM    ,pm9MM
# 	  MM     ,M YA.   ,A9 8M   MM `Mb    MM      `Mb    MM  8M   MM  MM   8M   MM
# 	.JMMmmmmMMM  `Ybmd9'  `Moo9^Yo.`Wbmd"MML.     `Wbmd"MML.`Moo9^Yo.`Mbmo`Moo9^Yo.
# 	
# 	
using PyPlot
using LinearAlgebra
using Statistics
using JLD2
using FileIO


rats=load("experimentDMP.jld2"); # load documents 
include("DiverseFunctions.jl") # load fucntions 

parameters=rats["parameters"]
featuresexperiment=rats["features"]; # [numberofrats, numberofdays, numberoftrials];
data=rats["data"];

# 	
# 	            ,,
# 	`7MM"""Mq.`7MM           mm             db      `7MMF'      `7MMF'          `7MMF'            db   MMP""MM""YMM
# 	  MM   `MM. MM           MM            ;MM:       MM          MM              MM             ;MM:  P'   MM   `7
# 	  MM   ,M9  MM  ,pW"Wq.mmMMmm         ,V^MM.      MM          MM              MM            ,V^MM.      MM
# 	  MMmmdM9   MM 6W'   `Wb MM          ,M  `MM      MM          MM              MM           ,M  `MM      MM
# 	  MM        MM 8M     M8 MM          AbmmmqMA     MM      ,   MM      ,       MM      ,    AbmmmqMA     MM
# 	  MM        MM YA.   ,A9 MM         A'     VML    MM     ,M   MM     ,M       MM     ,M   A'     VML    MM
# 	.JMML.    .JMML.`Ybmd9'  `Mbmo    .AMA.   .AMMA..JMMmmmmMMM .JMMmmmmMMM     .JMMmmmmMMM .AMA.   .AMMA..JMML.
# 	
# 	
ioff()
fig = figure("Test plot latencies",figsize=(9,9))
ax = fig[:add_subplot](1,1,1)

xlabel("trials")
ylabel("latencies")         

indexdays=collect(1:featuresexperiment[:numberofdays])
filter!(x->!(x in parameters[:indexprobedays]),indexdays)

indexdays=indexdays.-1; # just for iteration purposes

global labels_code=[] 

let iter=0; # because we miss probe days to plot the latencies 
	for k in indexdays
		# Calculate the lower value for the error bar : 
		uppererror = [std([data[n][k+1].day[i].latency for n in 1:numberofrats]; corrected=false) for i in 1:numberoftrials]./sqrt(numberofrats) ;
		lowererror = [std([data[n][k+1].day[i].latency for n in 1:numberofrats]; corrected=false) for i in 1:numberoftrials]./sqrt(numberofrats) ;
		errs=[lowererror,uppererror];
		PyPlot.plot(iter*numberoftrials.+(0:numberoftrials-1), [mean([data[n][k+1].day[i].latency for n in 1:numberofrats]) for i in 1:numberoftrials ], marker="None",linestyle="-",color="darkgreen",label="Base Plot")
		PyPlot.errorbar(iter*numberoftrials.+(0:numberoftrials-1),[mean([data[n][k+1].day[i].latency for n in 1:numberofrats]) for i in 1:numberoftrials ],yerr=errs,fmt="o",color="k")
		global labels_code=vcat(labels_code,collect(1:1:numberoftrials))
		iter+=1;
	end # end iteration on all days except probe days 
end # end let iter 

mx = matplotlib.ticker.MultipleLocator(1) # Define interval of minor ticks
ax.xaxis.set_major_locator(mx) # Set interval of minor ticks

ax.spines["top"].set_color("none")
ax.spines["right"].set_color("none")
ax.spines["bottom"].set_visible("False")
ax.spines["left"].set_visible("False")
  
ax[:set_ylim]([0,ymax])
xmin, xmax = ax.get_xlim() 
ymin, ymax = ax.get_ylim()
# get width and height of axes object to compute 
# matching arrowhead length and width
dps = fig.dpi_scale_trans.inverted()
bbox = ax.get_window_extent().transformed(dps)
width, height = bbox.width, bbox.height
# manual arrowhead width and length
hw = 1/20*(ymax-ymin) 
hl = 1/20*(xmax-xmin)
lw = 1 # axis line width
ohg = 0.3 # arrow overhang
# compute matching arrowhead length and width
yhw = hw/(ymax-ymin)*(xmax-xmin)* height/width 
yhl = hl/(xmax-xmin)*(ymax-ymin)* width/height
ax.arrow(xmin, ymin, xmax-xmin, 0.,length_includes_head= "True", fc="k", ec="k", lw = lw,head_width=hw, head_length=hl, overhang = ohg,  clip_on = "False") 

ax.arrow(xmin, ymin, 0., ymax-ymin,length_includes_head= "True", fc="k", ec="k", lw = lw, head_width=yhw, head_length=yhl, overhang = ohg,  clip_on = "False")
#ax.xaxis.set_major_locator(MaxNLocator(integer=True))

#rc("font", family="serif",size=16)
#title("One-shot learning in artificial watermaze")
xlabel("Trials ")#, fontsize=18);
ylabel("Latencies")#, fontsize=18)

#labels = [item.get_text() for item in ax.get_xticklabels()]
#labels[1] = labels_code
labels_code=vcat("","",labels_code)
ax.set_xticklabels(labels_code)

SMALL_SIZE = 10
MEDIUM_SIZE = 20
BIGGER_SIZE = 30

plt.rc("font", size=SMALL_SIZE)          # controls default text sizes
plt.rc("axes", titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc("axes", labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc("xtick", labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc("ytick", labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)    # legend fontsize
plt.rc("figure", titlesize=BIGGER_SIZE)  # fontsize of the figure title
savefig("Latencies_full.png")
#show()





# 	
# 	
# 	`7MM"""Mq.`7MMF'        .g8""8q.  MMP""MM""YMM                    MMP""MM""YMM `7MM"""Mq.  `7MMF'      db      `7MMF'       .M"""bgd
# 	  MM   `MM. MM        .dP'    `YM.P'   MM   `7                    P'   MM   `7   MM   `MM.   MM       ;MM:       MM        ,MI    "Y
# 	  MM   ,M9  MM        dM'      `MM     MM               ,AM            MM        MM   ,M9    MM      ,V^MM.      MM        `MMb.
# 	  MMmmdM9   MM        MM        MM     MM              AVMM            MM        MMmmdM9     MM     ,M  `MM      MM          `YMMNq.
# 	  MM        MM      , MM.      ,MP     MM            ,W' MM            MM        MM  YM.     MM     AbmmmqMA     MM      , .     `MM
# 	  MM        MM     ,M `Mb.    ,dP'     MM          ,W'   MM            MM        MM   `Mb.   MM    A'     VML    MM     ,M Mb     dM
# 	.JMML.    .JMMmmmmMMM   `"bmmd"'     .JMML.        AmmmmmMMmm        .JMML.    .JMML. .JMM..JMML..AMA.   .AMMA..JMMmmmmMMM P"Ybmmd"
# 	                                                         MM
# 	                                                         MM

indexdays=collect(1:numberofdays)
filter!(x->!(x in indexprobedays),indexdays)

# plot means latencies
using PyPlot

# create mean latencies per rats : 
Mean_Lat=[ [mean([data[indexrat][indexday].day[indextrial].latency for indexday in indexdays] ) for indextrial=1:numberoftrials ] for indexrat=1:numberofrats ];
# Calculate the error bar : 
uppererror = [std([Mean_Lat[indexrat][indextrial] for indexrat in 1:numberofrats]; corrected=false)./sqrt(numberofrats) for indextrial in 1:numberoftrials] ;
lowererror = [std([Mean_Lat[indexrat][indextrial] for indexrat in 1:numberofrats]; corrected=false)./sqrt(numberofrats) for indextrial in 1:numberoftrials] ;
errs=[lowererror,uppererror]; # gather 


using PyPlot

clf()
ioff()
fig = figure("Test plot latencies",figsize=(9,9))
ax=gca() 
#for indextrial=1:numberoftrialstest
PyPlot.plot((1:1:numberoftrials), [mean([Mean_Lat[indexrat][indextrial] for indexrat in 1:numberofrats]) for indextrial in 1:numberoftrials], marker="None",linestyle="-",color="darkgreen",label="Base Plot")
PyPlot.errorbar((1:1:numberoftrials),[mean([Mean_Lat[indexrat][indextrial] for indexrat in 1:numberofrats]) for indextrial in 1:numberoftrials],yerr=errs,fmt="o",color="k")

mx = matplotlib.ticker.MultipleLocator(1) # Define interval of minor ticks
ax.xaxis.set_major_locator(mx) # Set interval of minor ticks

ax.spines["top"].set_color("none")
ax.spines["right"].set_color("none")
ax.spines["bottom"].set_visible("False")
ax.spines["left"].set_visible("False")
  
ax[:set_ylim]([0,ymax])
xmin, xmax = ax.get_xlim() 
ymin, ymax = ax.get_ylim()
# get width and height of axes object to compute 
# matching arrowhead length and width
dps = fig.dpi_scale_trans.inverted()
bbox = ax.get_window_extent().transformed(dps)
width, height = bbox.width, bbox.height
# manual arrowhead width and length
hw = 1/20*(ymax-ymin) 
hl = 1/20*(xmax-xmin)
lw = 1 # axis line width
ohg = 0.3 # arrow overhang
# compute matching arrowhead length and width
yhw = hw/(ymax-ymin)*(xmax-xmin)* height/width 
yhl = hl/(xmax-xmin)*(ymax-ymin)* width/height
ax.arrow(xmin, ymin, xmax-xmin, 0.,length_includes_head= "True", fc="k", ec="k", lw = lw,head_width=hw, head_length=hl, overhang = ohg,  clip_on = "False") 

ax.arrow(xmin, ymin, 0., ymax-ymin,length_includes_head= "True", fc="k", ec="k", lw = lw, head_width=yhw, head_length=yhl, overhang = ohg,  clip_on = "False")
#ax.xaxis.set_major_locator(MaxNLocator(integer=True))
#rc("font", family="serif",size=16)
#title("One-shot learning in artificial watermaze")
xlabel("Trials ")#, fontsize=18);
ylabel("Latencies")#, fontsize=18)

savefig("Latencies_4trials.png")
#show()



