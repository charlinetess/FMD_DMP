# This file contains the code to implement the algorithm from Foster, Morris and Dayan (Hippocampus, 2000)  
using FileIO
using JLD2
#using Polynomials
using LinearAlgebra
using Statistics


# 	
# 	  ,,                          ,,          ,...                                      ,,
# 	`7MM                        `7MM        .d' ""                               mm     db
# 	  MM                          MM        dM`                                  MM
# 	  MM  ,pW"Wq.   ,6"Yb.   ,M""bMM       mMMmm`7MM  `7MM  `7MMpMMMb.  ,p6"bo mmMMmm `7MM  ,pW"Wq.`7MMpMMMb.  ,pP"Ybd
# 	  MM 6W'   `Wb 8)   MM ,AP    MM        MM    MM    MM    MM    MM 6M'  OO   MM     MM 6W'   `Wb MM    MM  8I   `"
# 	  MM 8M     M8  ,pm9MM 8MI    MM        MM    MM    MM    MM    MM 8M        MM     MM 8M     M8 MM    MM  `YMMMa.
# 	  MM YA.   ,A9 8M   MM `Mb    MM        MM    MM    MM    MM    MM YM.    ,  MM     MM YA.   ,A9 MM    MM  L.   I8
# 	.JMML.`Ybmd9'  `Moo9^Yo.`Wbmd"MML.    .JMML.  `Mbod"YML..JMML  JMML.YMbmd'   `Mbmo.JMML.`Ybmd9'.JMML  JMML.M9mmmP'
# 	
# 	

include("DiverseFunctions.jl")

# 	
# 	
# 	`7MM"""Mq.   db      `7MM"""Mq.        db      `7MMM.     ,MMF'`7MM"""YMM MMP""MM""YMM `7MM"""YMM  `7MM"""Mq.   .M"""bgd
# 	  MM   `MM. ;MM:       MM   `MM.      ;MM:       MMMb    dPMM    MM    `7 P'   MM   `7   MM    `7    MM   `MM. ,MI    "Y
# 	  MM   ,M9 ,V^MM.      MM   ,M9      ,V^MM.      M YM   ,M MM    MM   d        MM        MM   d      MM   ,M9  `MMb.
# 	  MMmmdM9 ,M  `MM      MMmmdM9      ,M  `MM      M  Mb  M' MM    MMmmMM        MM        MMmmMM      MMmmdM9     `YMMNq.
# 	  MM      AbmmmqMA     MM  YM.      AbmmmqMA     M  YM.P'  MM    MM   Y  ,     MM        MM   Y  ,   MM  YM.   .     `MM
# 	  MM     A'     VML    MM   `Mb.   A'     VML    M  `YM'   MM    MM     ,M     MM        MM     ,M   MM   `Mb. Mb     dM
# 	.JMML. .AMA.   .AMMA..JMML. .JMM..AMA.   .AMMA..JML. `'  .JMML..JMMmmmmMMM   .JMML.    .JMMmmmmMMM .JMML. .JMM.P"Ybmmd"
# 	
# 	


# Creating the circle and the place cells:
center=[0,0];
R= 100; # Radius of the circle in cm
r=6;# Radius of the platform  in cm
radiussearchpref=20; # radius of the area in which we calculate searchpreference 

# Motion characteristic 
dt=0.1; # timestep in s 
speed=30; # speed of the rat in cm.s-1
angles=[-3*pi/4, -2*pi/4, -pi/4, 0, pi/4, 2*pi/4, 3*pi/4, pi];# Different possible directions 

# Trial characteristic :
T=120; # maximal duration of a trial in seconds
Tprobetrials=60; # in probe trials the reward is removed and we study the behavior of the agent during 60 seconds 


# Number of Place cells 
NPC=493; # number of place cells 
#Xplacecell=sunflower(N,R,2)[:,1]; # absciss place cells  
#Yplacecell=sunflower(N,R,2)[:,2]; # y place cells 

# Place cell # initialize the centres of the place cells by random unifrom sampling across the pool
arguments= rand(1,NPC)*2*pi; # random angles covering the whole circle 
radii= sqrt.(rand(1,NPC))*R; # initialize the centres of the place cells by random unifrom sampling across the pool
radii=sort(radii,dims=2,rev=true); # sort the radius in increasing orders 
centres= [cos.(arguments).*radii; sin.(arguments).*radii]; 

σPC=0.30*100; # variability of place cell activity, in centimeters

# Action cells : 
NA=9; # number of action cells 


# Potential Starting positions of the rat :
Xstart=[0.95,0,-0.95,0].*R; # East, North, West, South
Ystart=[0,0.95,0,-0.95].*R;

times=collect(0:dt:T+dt);

Xplacecell=centres[1,:];
Yplacecell=centres[2,:];

# Potential positions of the platform : 
Xplatform=[0.3,0,-0.3,0,0.5,-0.5,0.5,-0.5].*R; # in cm
Yplatform=[0,0.3,0,-0.3,0.5,0.5,-0.5,-0.5].*R;# in cm

# Potential Starting positions of the rat :
Xstart=[0.95,0,-0.95,0].*R; # East, North, West, South
Ystart=[0,0.95,0,-0.95].*R;


# Parameter that regulate the choice between former angle and new angle 
momentum=1.0;

# Learning variables : 
γ=0.99; # Discount factor.  they dont precise the value  
actorLR=0.01; # actor learning rate for normal action cells 
actor2LR=0.0009; # actor learning rate for additional action cells 
criticLR=0.01; # critic learning rate
# learning rate for position:
LRxcoord=0.001; # learning rate for x coordinate 
LRycoord=0.001;  # learning rate for y coordinate 

# eligibility parameter for postion estimation 
λ=0.9;

# Probe days : 
indexprobedays=[1, 7, 10];

temperature=2;

randomstartingpos=0; # 1 if random starting position at every trial, 0 if every start loc are used every day


# Define number of rats, number of days and numbers of trials per day
numberofdays=10; # every day goal location changes 
numberofrats=20; # number of independant runs 
numberoftrials=4; # number of trials per goal location 



parameters=Dict(:momentum=>momentum,:γ=>γ,:actorLR=>actorLR,:criticLR=>criticLR,:centres=>centres,:R=>R,:LRxcoord=>LRxcoord,:LRycoord=>LRycoord,:λ=>λ,:actor2LR=>actor2LR,:r=>r,:speed=>speed,:angles=>angles,:NPC=>NPC,:NA=>NA,:σPC=>σPC,:Xstart=>Xstart,:Ystart=>Ystart,:dt=>dt,:T=>T,:times=>times,:Xplatform=>Xplatform,:Yplatform=>Yplatform,:temperature=>temperature,:randomstartingpos=>randomstartingpos,:indexprobedays=>indexprobedays,:Tprobetrials=>Tprobetrials,:radiussearchpref=>radiussearchpref);
featuresexperiment=Dict(:numberofrats=>numberofrats, :numberofdays=>numberofdays, :numberoftrials=>numberoftrials);

NameOfFile="experimentDMP.jld2";


# 	
# 	                    ,,    ,,      ,...                                      ,,
# 	                  `7MM  `7MM    .d' ""                               mm     db
# 	                    MM    MM    dM`                                  MM
# 	 ,p6"bo   ,6"Yb.    MM    MM   mMMmm`7MM  `7MM  `7MMpMMMb.  ,p6"bo mmMMmm `7MM  ,pW"Wq.`7MMpMMMb.
# 	6M'  OO  8)   MM    MM    MM    MM    MM    MM    MM    MM 6M'  OO   MM     MM 6W'   `Wb MM    MM
# 	8M        ,pm9MM    MM    MM    MM    MM    MM    MM    MM 8M        MM     MM 8M     M8 MM    MM
# 	YM.    , 8M   MM    MM    MM    MM    MM    MM    MM    MM YM.    ,  MM     MM YA.   ,A9 MM    MM
# 	 YMbmd'  `Moo9^Yo..JMML..JMML..JMML.  `Mbod"YML..JMML  JMML.YMbmd'   `Mbmo.JMML.`Ybmd9'.JMML  JMML.
# 	
# 	
@time begin # get the time it takes to run it 
DMP(parameters,featuresexperiment,NameOfFile)
end # end time 






