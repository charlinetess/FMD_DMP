# This file loads the date produced from the file DMPCoordSearchPref and generates a plot of the search preference on all search rpef days and a plot of all platform positions 

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

rats=load("experimentDMP.jld2"); # load documents 
include("DiverseFunctions.jl") # load fucntions 

parameters=rats["parameters"]
featuresexperiment=rats["features"]; # [numberofrats, numberofdays, numberoftrials];
data=rats["data"];


argument=-pi/50:pi/50:2*pi+pi/50;




# 	
# 	           ,,
# 	         `7MM           mm
# 	           MM           MM
# 	`7MMpdMAo. MM  ,pW"Wq.mmMMmm     M"""MMV ,pW"Wq.`7MMpMMMb.  .gP"Ya  ,pP"Ybd
# 	  MM   `Wb MM 6W'   `Wb MM       '  AMV 6W'   `Wb MM    MM ,M'   Yb 8I   `"
# 	  MM    M8 MM 8M     M8 MM         AMV  8M     M8 MM    MM 8M"""""" `YMMMa.
# 	  MM   ,AP MM YA.   ,A9 MM        AMV  ,YA.   ,A9 MM    MM YM.    , L.   I8
# 	  MMbmmd'.JMML.`Ybmd9'  `Mbmo    AMMmmmM `Ybmd9'.JMML  JMML.`Mbmmd' M9mmmP'
# 	  MM
# 	.JMML.
using PyPlot
for k=1:length(parameters[:Xplatform])
plot(parameters[:Xplatform][k].+parameters[:r].*cos.(argument),parameters[:Yplatform][k].+parameters[:r].*sin.(argument),"m-") # plot platform
plot(parameters[:Xplatform][k].+parameters[:radiussearchpref].*cos.(argument),parameters[:Yplatform][k].+parameters[:radiussearchpref].*sin.(argument),"b-")
# Plot circle
plot(parameters[:R]*cos.(argument),parameters[:R]*sin.(argument),"k-")
ax7=gca() 
ax7[:set_axis_off]()
end
savefig("Platforms_Zones.png")


# using PyPlot
# scatter(parameters[:centres][1,:],parameters[:centres][2,:])
# # Plot circle
# plot(parameters[:R]*cos.(argument),parameters[:R]*sin.(argument),"k-")
# ax7=gca() 
# ax7[:set_axis_off]()
# show()

# 	
# 	           ,,                                                                 ,,                                           ,...
# 	         `7MM           mm        .M"""bgd                                  `7MM            `7MM"""Mq.                   .d' ""
# 	           MM           MM       ,MI    "Y                                    MM              MM   `MM.                  dM`
# 	`7MMpdMAo. MM  ,pW"Wq.mmMMmm     `MMb.      .gP"Ya   ,6"Yb.  `7Mb,od8 ,p6"bo  MMpMMMb.        MM   ,M9 `7Mb,od8 .gP"Ya  mMMmm
# 	  MM   `Wb MM 6W'   `Wb MM         `YMMNq. ,M'   Yb 8)   MM    MM' "'6M'  OO  MM    MM        MMmmdM9    MM' "',M'   Yb  MM
# 	  MM    M8 MM 8M     M8 MM       .     `MM 8M""""""  ,pm9MM    MM    8M       MM    MM        MM         MM    8M""""""  MM
# 	  MM   ,AP MM YA.   ,A9 MM       Mb     dM YM.    , 8M   MM    MM    YM.    , MM    MM        MM         MM    YM.    ,  MM
# 	  MMbmmd'.JMML.`Ybmd9'  `Mbmo    P"Ybmmd"   `Mbmmd' `Moo9^Yo..JMML.   YMbmd'.JMML  JMML.    .JMML.     .JMML.   `Mbmmd'.JMML.
# 	  MM
# 	.JMML.

searchprefs=[mean([data[n][k].day[2].SearchPref for n in 1:featuresexperiment[:numberofrats]]) for k in parameters[:indexprobedays] ]; # compute mean, the probe trial is always trial2
# plot bar plot 
using PyPlot
ioff()
fig = figure("Test plot search Preference",figsize=(4,9))
ax = gca()

SMALL_SIZE = 15
MEDIUM_SIZE = 30
BIGGER_SIZE = 30

plt.rc("font", size=SMALL_SIZE)          # controls default text sizes
plt.rc("axes", titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc("axes", labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc("xtick", labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc("ytick", labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)    # legend fontsize
plt.rc("figure", titlesize=BIGGER_SIZE)  # fontsize of the figure title

uppererror = [std([data[n][k].day[2].SearchPref*100 for n in 1:featuresexperiment[:numberofrats]]) for k in parameters[:indexprobedays] ]./sqrt(numberofrats); 
lowererror = [std([data[n][k].day[2].SearchPref*100 for n in 1:featuresexperiment[:numberofrats]]) for k in parameters[:indexprobedays] ]./sqrt(numberofrats); 
errs=[lowererror,uppererror];

bar(1:1:length(parameters[:indexprobedays]),[mean([data[n][k].day[2].SearchPref*100 for n in 1:featuresexperiment[:numberofrats]]) for k in parameters[:indexprobedays] ],width=0.5,yerr=errs,color=[60/255,179/255,113/255],align="center",alpha=0.4)
ax[:axes][:get_xaxis]()[:set_ticks]([])
#xlabel=["day $(indexprobedays[1])","day $(indexprobedays[2])","day $(indexprobedays[3])"]
axis("tight")
PyPlot.plot(0:length(parameters[:indexprobedays])+2*0.2/2+1,12.5*ones(size(0:length(parameters[:indexprobedays])+2*0.2/2+1,1),size(0:length(parameters[:indexprobedays])+2*0.2/2+1,2)),color="green",label="Chance",linestyle="--")

locs, labels = xticks()  
xticks(1:1:length(parameters[:indexprobedays]), ["day $(parameters[:indexprobedays][1])","day $(parameters[:indexprobedays][2])","day $(parameters[:indexprobedays][3])"])

ax[:grid](false);
ax.spines["right"].set_visible(false)
ax.spines["top"].set_visible(false)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
legend()
title("%time spent in the correct zone")

# show()

# 	
# 	                                     ,...,,
# 	                                   .d' ""db
# 	                                   dM`
# 	,pP"Ybd  ,6"Yb.`7M'   `MF'.gP"Ya  mMMmm`7MM  .P"Ybmmm
# 	8I   `" 8)   MM  VA   ,V ,M'   Yb  MM    MM :MI  I8
# 	`YMMMa.  ,pm9MM   VA ,V  8M""""""  MM    MM  WmmmP"
# 	L.   I8 8M   MM    VVV   YM.    ,  MM    MM 8M
# 	M9mmmP' `Moo9^Yo.   W     `Mbmmd'.JMML..JMML.YMMMMMb
# 	                                            6'     dP
# 	                                            Ybmmmd'

savefig("SearchPref.png")

