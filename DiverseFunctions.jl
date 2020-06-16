# This file contains all functions that are shared with both data-generation computations and plotting codes related to Foster, Morris and Dayan (Hippocampus, 2000) paper modelling an artificial agent wihin a morris watermaze spatial navigation task 


#
#
#     `7MM"""Mq.`7MMF'            db       .g8"""bgd `7MM"""YMM        .g8"""bgd `7MM"""YMM  `7MMF'      `7MMF'       .M"""bgd
#       MM   `MM. MM             ;MM:    .dP'     `M   MM    `7      .dP'     `M   MM    `7    MM          MM        ,MI    "Y
#       MM   ,M9  MM            ,V^MM.   dM'       `   MM   d        dM'       `   MM   d      MM          MM        `MMb.
#       MMmmdM9   MM           ,M  `MM   MM            MMmmMM        MM            MMmmMM      MM          MM          `YMMNq.
#       MM        MM      ,    AbmmmqMA  MM.           MM   Y  ,     MM.           MM   Y  ,   MM      ,   MM      , .     `MM
#       MM        MM     ,M   A'     VML `Mb.     ,'   MM     ,M     `Mb.     ,'   MM     ,M   MM     ,M   MM     ,M Mb     dM
#     .JMML.    .JMMmmmmMMM .AMA.   .AMMA. `"bmmmd'  .JMMmmmmMMM       `"bmmmd'  .JMMmmmmMMM .JMMmmmmMMM .JMMmmmmMMM P"Ybmmd"
#
#

function  placecells(pos,cent,width)
	# pos vector size 1,2
	# cent vector size (Number of place cells, 2), centers of place cells 
	# width is a positive scalar that refers to the width of place cells activity profile
#
# PLACECELLS(POSITION,CENTRES,WIDTH) calculates the activity of the place cells
#in the simulation. The returned vector F is of length N, where N is the number of place
#cells, and it contains the activity of each place cell given the simulated rat's current
#POSITION (a 2 element column vector). The activity of the place cells is modelled as a
#rate-of-fire (i.e. a scalar value) determined by a gaussian function. The CENTRES of the
#gaussian functions are an argument, and must be a 2 x N matrix containing each place
#cell's preferred location in 2D space. The WIDTH of the place cell fields must
#also be provided as a scalar value (all place cells are assumed to have the same
#width).
#
#The returned vector, F, must be a N element column vector.
    # calculate place cell activity

F = exp.(-sum((repeat(pos,1,size(cent,2))-cent).^2,dims=1)./(2*width.^2));
Fbis=zeros(length(F),1)
transpose!(Fbis,F)
    return Fbis
end

#
#                                                                 ,,
#     `7MM"""Mq.                                                `7MM
#       MM   `MM.                                                 MM
#       MM   ,M9  .gP"Ya `7M'    ,A    `MF',6"Yb.  `7Mb,od8  ,M""bMM
#       MMmmdM9  ,M'   Yb  VA   ,VAA   ,V 8)   MM    MM' "',AP    MM
#       MM  YM.  8M""""""   VA ,V  VA ,V   ,pm9MM    MM    8MI    MM
#       MM   `Mb.YM.    ,    VVV    VVV   8M   MM    MM    `Mb    MM
#     .JMML. .JMM.`Mbmmd'     W      W    `Moo9^Yo..JMML.   `Wbmd"MML.
#
#

# Calculate reward as a function of position, reward is given at the position of the goal 
function reward(x,y,xp,yp,r) # x,y position of the rat and xp,yp position of the platform, r radius of the platform
    if (x-xp)^2+(y-yp)^2<= r^2 # if the rat is in the platform
        R=1;
    else # else 
        R=0;
    end 
end

#
#       ,,                    ,,    ,,
#       db                  `7MM    db
#                             MM
#     `7MM  `7MMpMMMb.   ,M""bMM  `7MM  ,p6"bo   .gP"Ya
#       MM    MM    MM ,AP    MM    MM 6M'  OO  ,M'   Yb
#       MM    MM    MM 8MI    MM    MM 8M       8M""""""
#       MM    MM    MM `Mb    MM    MM YM.    , YM.    ,
#     .JMML..JMML  JMML.`Wbmd"MML..JMML.YMbmd'   `Mbmmd'
#
#

# This function takes as input a vector that defines a probability distribution and draws a random number according to this probability distribution. 

function indice(A)
    Acum=[sum(A[1:k]) for k=1:length(A)] # Compute summed probability distribution:
    x=rand(1)[1]
    for i=1:length(Acum)
       if i==1
           if x<Acum[i] # if the random number generated is before the first 
                return i
            end
        else
            if Acum[i-1]<x<=Acum[i]
                return i
            end
        end
    end   
end

#
#                                                          ,,
#     `7MM"""YMM                                           db                                       mm
#       MM    `7                                                                                    MM
#       MM   d    `7M'   `MF'`7MMpdMAo.  .gP"Ya `7Mb,od8 `7MM  `7MMpMMMb.pMMMb.  .gP"Ya `7MMpMMMb.mmMMmm
#       MMmmMM      `VA ,V'    MM   `Wb ,M'   Yb  MM' "'   MM    MM    MM    MM ,M'   Yb  MM    MM  MM
#       MM   Y  ,     XMX      MM    M8 8M""""""  MM       MM    MM    MM    MM 8M""""""  MM    MM  MM
#       MM     ,M   ,V' VA.    MM   ,AP YM.    ,  MM       MM    MM    MM    MM YM.    ,  MM    MM  MM
#     .JMMmmmmMMM .AM.   .MA.  MMbmmd'   `Mbmmd'.JMML.   .JMML..JMML  JMML  JMML.`Mbmmd'.JMML  JMML.`Mbmo
#                              MM
#                            .JMML.


function DMP(parameters,featuresexperiment,NameOfFile)
    experiment=[];
    for indexrat=1:featuresexperiment[:numberofrats]
        let criticweights=zeros(parameters[:NPC],1), actorweights=zeros(parameters[:NPC],parameters[:NA]), weightsxcoord=zeros(parameters[:NPC],1), weightsycoord=zeros(parameters[:NPC],1), currentexperiment=[];# Initialisation variables   
                    ##########  ##########  ##########  ##########   ########## 
                ##########  ##########  START EXPERIMENT  ##########  ##########  
                    ##########  ##########  ##########  ##########   ########## 
                for indexday=1:featuresexperiment[:numberofdays]
                    let indexplatform=rand(1:length(parameters[:Xplatform])),xp=parameters[:Xplatform][indexplatform],yp=parameters[:Yplatform][indexplatform], currentday=[],platform=0,Xplatformestimate=0,Yplatformestimate=0; # Everyday the location of the platofrm changes, chosen randomly. platform=0 indicates that the coordinate action defines now random motion, the estimates of the platform location are reset. 
                                ##########  ##########  ##########  ##########  
                            ##########  ##########  START DAY ##########  ##########  
                                ##########  ##########  ##########  ##########  
                            for indextrial=1:featuresexperiment[:numberoftrials]  
                                let positionstart,indexstart,currentposition, re,k,t,historyX,historyY,TDerrors, searchpref, searchinzones,arg,timeout, prevdir, indexaction, probaction9, X, Y, C,actplacecell, formeractplacecell, Xestimate, Yestimate, actactioncell, Pactioncell, probaction9, SumPactioncell, formerposition, Xf, Yf, err, eligibilitytrace
                                    if randomstartingpos==1  ## Chose starting position :
                                        # just to try if it learns better
                                        indexstart=rand(1:length(parameters[:Xstart])); # chose randomnly between 4 possibilities 1 East 2 North 3 West 4 South
                                    elseif randomstartingpos==0
                                        indexstart=mod(indextrial+indexrat+indexday,4)+1; # use all the 4 starting position everyday 
                                    end
                                    positionstart=[parameters[:Xstart][indexstart] parameters[:Ystart][indexstart]];# take indexstart-th starting position 
                                    currentposition=positionstart;
                                    re=0; # Initialize reward 
                                    k=1; # Initialise index to save the trajectory and the values 
                                    t=parameters[:times][k]; # initialise time
                                    historyX=Float64[]; # to store the trajectory
                                    historyY=Float64[];
                                    TDerrors=Float64[];
                                    searchpref=0; # initialise search preference measure (amount of time spent in an area centetred around the platform)
                                    searchinzones=0;
                                    arg=0;     
                                    timeout=0;   
                                    prevdir=[0 0];
                                    indexaction=0; 
                                    probaction9=[];
                                    if (indexday in parameters[:indexprobedays])&&(indextrial==2) # if we are on probe days, the second trial is different than the others 
                                        while t<=parameters[:Tprobetrials] # we put no learning on actorweights and criticweights on probe trials
                                                if t==parameters[:Tprobetrials]
                                                    X=xp;
                                                    Y=yp;
                                                    currentposition=[X Y];
                                                    timeout=1; # if we have to put the rat on the platform then we dont reinforce the actor but only the critic
                                                    platform=1;
                                                    Xplatformestimate=dot(weightsxcoord,placecells([currentposition[1],currentposition[2]],parameters[:centres],parameters[:σPC])); # we register our estimate of the position of the paltform
                                                    Yplatformestimate=dot(weightsycoord,placecells([currentposition[1],currentposition[2]],parameters[:centres],parameters[:σPC]));
                                                end                     
                                            push!(historyX,currentposition[1]); 
                                            push!(historyY,currentposition[2]);
                                            # compute new activity of pace cells :
                                            actplacecell=placecells([currentposition[1],currentposition[2]],parameters[:centres],parameters[:σPC]);  
                                                if !(k==1)
                                                    formeractplacecell=actplacecell; # need storing to compute the self motion estimate
                                                end                      
                                            ### Compute Critic ###
                                            C=dot(criticweights,actplacecell); # current estimation of the future discounted reward         
                                            Xestimate=dot(weightsxcoord,actplacecell); # estimate position 
                                            Yestimate=dot(weightsycoord,actplacecell);
                                            positionestimate=[Xestimate Yestimate];
                                            ####### Take decision and move to new position : ######## 
                                            actactioncell=transpose(actorweights)*actplacecell; # Compute action cells activity actorweights contains place cells in rows and action cells in column 
                                                if maximum(abs.(actactioncell))>=100
                                                    actactioncell=100*actactioncell./maximum(abs.(actactioncell)); 
                                                end
                                            Pactioncell=exp.(parameters[:temperature]*actactioncell)./sum(exp.(parameters[:temperature]*actactioncell)); # Compute probability distribution : 
                                            probaction9=push!(probaction9, Pactioncell[parameters[:NA]]); # we store the probability of action 9, for plot 
                                            indexaction=indice(Pactioncell);# draws a random number according to Pactioncel
                                                if indexaction==parameters[:NA] # if we chose the acoord action
                                                        if platform==0 # if we havent registered the platform position yet 
                                                                secondindexaction=rand(1:(parameters[:NA]-1))
                                                                secondindexaction=secondindexaction[1]
                                                                argdecision=parameters[:angles][Int(secondindexaction)]; # compute the coreesponding angle 
                                                                newdir=[cos(argdecision) sin(argdecision)];
                                                                dir=(newdir./(1.0+parameters[:momentum]).+parameters[:momentum].*prevdir./(1.0+parameters[:momentum]));
                                                                if !(norm(dir)==0)
                                                                    dir=dir./norm(dir);
                                                                end
                                                            elseif platform==1 # if we have registered the platform position
                                                                dir=[Xplatformestimate Yplatformestimate].-positionestimate; # get the vector of displacement 
                                                                if !(norm(dir)==0)
                                                                    dir=dir./norm(dir);
                                                                end                
                                                        end                                        
                                                    else # if indexaction is one of the 8th first indexes
                                                        argdecision=parameters[:angles][Int(indexaction)]; # compute the coreesponding angle 
                                                        newdir=[cos(argdecision) sin(argdecision)];
                                                        dir=(newdir./(1.0+parameters[:momentum]).+parameters[:momentum].*prevdir./(1.0+parameters[:momentum]));
                                                        if !(norm(dir)==0)
                                                            dir=dir./norm(dir);
                                                        end                     
                                                end 
                                            # Store former position 
                                            formerposition=currentposition;
                                            # Compute new position : 
                                            currentposition=currentposition.+parameters[:dt].*parameters[:speed].*dir; 
                                                if currentposition[1]^2+currentposition[2]^2>=parameters[:R]^2
                                                    currentposition = (currentposition./norm(currentposition))*(parameters[:R] - parameters[:R]/50);
                                                    if !(norm(currentposition-formerposition)==0)
                                                        dir=(currentposition-formerposition)./norm(currentposition-formerposition);
                                                    else
                                                        dir=[0 0]
                                                    end
                                                end
                                            prevdir=dir;                            
                                            X=currentposition[1];
                                            Y=currentposition[2];
                                            Xf=formerposition[1];
                                            Yf=formerposition[2];
                                                ####### ####### ####### Updating search preference  ####### ####### #######
                                                if (currentposition[1]-xp)^2+(currentposition[2]-yp)^2<= radiussearchpref^2   
                                                    searchpref=searchpref+1*parameters[:dt];
                                                end
                                                # compute time in the 8 zones :
                                                if ((currentposition[1]-parameters[:Xplatform][1])^2+(currentposition[2]-parameters[:Yplatform][1])^2<= parameters[:radiussearchpref]^2)|((currentposition[1]-parameters[:Xplatform][2])^2+(currentposition[2]-parameters[:Yplatform][2])^2<= parameters[:radiussearchpref]^2)|((currentposition[1]-parameters[:Xplatform][3])^2+(currentposition[2]-parameters[:Yplatform][3])^2<= parameters[:radiussearchpref]^2)|((currentposition[1]-parameters[:Xplatform][4])^2+(currentposition[2]-parameters[:Yplatform][4])^2<= parameters[:radiussearchpref]^2)|((currentposition[1]-parameters[:Xplatform][5])^2+(currentposition[2]-parameters[:Yplatform][5])^2<= parameters[:radiussearchpref]^2)|((currentposition[1]-parameters[:Xplatform][6])^2+(currentposition[2]-parameters[:Yplatform][6])^2<= parameters[:radiussearchpref]^2)|((currentposition[1]-parameters[:Xplatform][7])^2+(currentposition[2]-parameters[:Yplatform][7])^2<= parameters[:radiussearchpref]^2)|((currentposition[1]-parameters[:Xplatform][8])^2+(currentposition[2]-parameters[:Yplatform][8])^2<= parameters[:radiussearchpref]^2)
                                                    searchinzones=searchinzones+1*parameters[:dt];
                                                end
                                            # compute new activity of pace cells :
                                            actplacecell=placecells([currentposition[1],currentposition[2]],parameters[:centres],parameters[:σPC]);
                                            ###  Compute reward ### 
                                            re=reward(currentposition[1],currentposition[2],xp,yp,r);   
                                                if re==1 # if we are on the platform 
                                                   ###  Compute error ###
                                                    Cnext=0;
                                                    platform=1;
                                                    Xplatformestimate=dot(weightsxcoord,placecells([currentposition[1],currentposition[2]],parameters[:centres],parameters[:σPC])); # we register our estimate of the position of the paltform
                                                    Yplatformestimate=dot(weightsycoord,placecells([currentposition[1],currentposition[2]],parameters[:centres],parameters[:σPC]));
                                                    else 
                                                    Cnext=dot(criticweights,actplacecell);# new estimation of the future discounted reward 
                                                end 
                                            #### Compute error  ####
                                            err=re+parameters[:γ]*Cnext-C[1];
                                            push!(TDerrors,err);
                                            k=k+1;
                                            t=parameters[:times][k];
                                              ######### Compute new weights : ########
                                            # we dont learn anything about critic or actor during probe trials 
                                            # update the weight for position estimate 
                                                if k==2
                                                    eligibilitytrace=actplacecell;
                                                end
                                                if !(k==2)
                                                    # self motion estimate : 
                                                    eligibilitytrace=parameters[:λ]*eligibilitytrace+actplacecell;
                                                    deltax=dot(weightsxcoord,actplacecell)-dot(weightsxcoord,formeractplacecell); # how much the former weights estimate my motion         
                                                    deltay=dot(weightsycoord,actplacecell)-dot(weightsycoord,formeractplacecell);          
                                                    weightsxcoord=weightsxcoord.+parameters[:LRxcoord].*(-parameters[:dt].*parameters[:speed].*dir[1]+deltax).*eligibilitytrace;        
                                                    weightsycoord=weightsycoord.+parameters[:LRycoord].*(-parameters[:dt].*parameters[:speed].*dir[2]+deltay).*eligibilitytrace;              
                                                end     
                                        end # end while 
                                    else
                                        searchpref=0;
                                        while t<=T && re==0                                           
                                                if t==T
                                                    X=xp;
                                                    Y=yp;
                                                    currentposition=[X Y];
                                                    timeout=1; # if we have to put the rat on the platform then we dont reinforce the actor but only the critic
                                                    platform=1;
                                                    Xplatformestimate=dot(weightsxcoord,placecells([currentposition[1],currentposition[2]],centres,parameters[:σPC])); # we register our estimate of the position of the paltform
                                                    Yplatformestimate=dot(weightsycoord,placecells([currentposition[1],currentposition[2]],centres,parameters[:σPC]));
                                                end                     
                                            # Store former position to be able to draw trajectory
                                            push!(historyX,currentposition[1]); 
                                            push!(historyY,currentposition[2]);
                                            # compute new activity of pace cells :
                                            actplacecell=placecells([currentposition[1],currentposition[2]],centres,parameters[:σPC]);  
                                            if !(k==1)
                                                formeractplacecell=actplacecell; # need storing to compute the self motion estimate
                                            end
                                            ### Compute Critic ###
                                            C=dot(criticweights,actplacecell); # current estimation of the future discounted reward                                     
                                            # estimate position 
                                            Xestimate=dot(weightsxcoord,actplacecell);
                                            Yestimate=dot(weightsycoord,actplacecell);
                                            positionestimate=[Xestimate Yestimate];
                                            ####### Take decision and move to new position : ########
                                            #  Compute action cell activity    
                                            actactioncell=transpose(actorweights)*actplacecell; # careful actorweights contains place cells in rows and action cells in column 
                                                if maximum(abs.(actactioncell))>=100
                                                    actactioncell=100*actactioncell./maximum(abs.(actactioncell)); 
                                                end
                                            # Compute probability distribution : 
                                            Pactioncell=exp.(parameters[:temperature]*actactioncell)./sum(exp.(parameters[:temperature]*actactioncell)); 
                                            probaction9=push!(probaction9, Pactioncell[parameters[:NA]]); # we store the probability of action 9 
                                            indexaction=indice(Pactioncell); 
                                            if indexaction==parameters[:NA] # if we chose the acoord action
                                                if platform==0 # if we havent registered the platform position yet 
                                                    secondindexaction=rand(1:(parameters[:NA]-1))
                                                    secondindexaction=secondindexaction[1]
                                                    argdecision=parameters[:angles][Int(secondindexaction)]; # compute the coreesponding angle 
                                                    newdir=[cos(argdecision) sin(argdecision)];
                                                    dir=(newdir./(1.0+parameters[:momentum]).+parameters[:momentum].*prevdir./(1.0+parameters[:momentum]));
                                                    if !(norm(dir)==0)
                                                        dir=dir./norm(dir);
                                                    end
                                                elseif platform==1 # if we have registered the platform position
                                                    dir=[Xplatformestimate Yplatformestimate].-positionestimate; # get the vector of displacement 
                                                    if !(norm(dir)==0)
                                                        dir=dir./norm(dir);
                                                    end                
                                                end                                        
                                            else # if indexaction is one of the 8th first indexes
                                                argdecision=parameters[:angles][Int(indexaction)]; # compute the coreesponding angle 
                                                newdir=[cos(argdecision) sin(argdecision)];
                                                dir=(newdir./(1.0+parameters[:momentum]).+parameters[:momentum].*prevdir./(1.0+parameters[:momentum]));
                                                if !(norm(dir)==0)
                                                    dir=dir./norm(dir);
                                                end                     
                                            end   
                                            # Store former position 
                                            formerposition=currentposition;
                                            # Compute new position : 
                                            currentposition=currentposition.+parameters[:dt].*parameters[:speed].*dir; 
                                                if currentposition[1]^2+currentposition[2]^2>=parameters[:R]^2
                                                    currentposition = (currentposition./norm(currentposition))*(parameters[:R] - parameters[:R]/50);
                                                    if !(norm(currentposition-formerposition)==0)
                                                        dir=(currentposition-formerposition)./norm(currentposition-formerposition);
                                                    else
                                                        dir=[0 0]
                                                    end
                                                end
                                            prevdir=dir;                            

                                            X=currentposition[1];
                                            Y=currentposition[2];
                                            Xf=formerposition[1];
                                            Yf=formerposition[2];
                                            # compute new activity of pace cells :
                                            actplacecell=placecells([currentposition[1],currentposition[2]],parameters[:centres],parameters[:σPC]);
                                            ###  Compute reward ### 
                                            re=reward(currentposition[1],currentposition[2],xp,yp,parameters[:r]);   
                                            if re==1 # if we are on the platform 
                                               ###  Compute error ###
                                                Cnext=0;
                                                platform=1;
                                                Xplatformestimate=dot(weightsxcoord,placecells([currentposition[1],currentposition[2]],parameters[:centres],parameters[:σPC])); # we register our estimate of the position of the paltform
                                                Yplatformestimate=dot(weightsycoord,placecells([currentposition[1],currentposition[2]],parameters[:centres],parameters[:σPC]));
                                            else 
                                                Cnext=dot(criticweights,actplacecell);# new estimation of the future discounted reward 
                                            end 
                                            #### Compute error  ####
                                            err=re+parameters[:γ]*Cnext-C[1];
                                    
                                            push!(TDerrors,err);
                                            ######### Compute new weights : ########
                                                if timeout==0
                                                    G=zeros(parameters[:NA],1);
                                                    G[indexaction]=1;
                                                    # weights between action cells and place cells only reinforced when the rats actually found the platform
                                                    # z[:,indexaction]=z[:,indexaction]+Z.*err.*actplacecell; # only the weights between place cells and the action taken are updated
                                                    if indexaction==parameters[:NA] # if we chose the last action we dont update it with the product of place cella ctivity but only with the error 
                                                        # we update this only if it has been specifically using the stored location of platform
                                                        if platform==1
                                                            actorweights=actorweights+parameters[:actor2LR].*err.*ones(size(actplacecell,1),size(actplacecell,2))*transpose(G); 
                                                        end
                                                    else 
                                                    
                                                        actorweights=actorweights+parameters[:actorLR].*err.*actplacecell*transpose(G); 
                                                    end
                                                    
                                                end
                                            
                                            # weights between critic and place cells :
                                            # Save value to draw valuemap
                                            # push!(valuemap,w);
                                            criticweights=criticweights+parameters[:criticLR].*err.*actplacecell;
                                            k=k+1;
                                            t=times[k];
                                            # update the weight for position estimate 
                                            if (k==2)
                                                eligibilitytrace=actplacecell;
                                            end
                                            if !(k==2)
                                                # self motion estimate : 
                                                eligibilitytrace=parameters[:λ]*eligibilitytrace+actplacecell;
                                                deltax=dot(weightsxcoord,actplacecell)-dot(weightsxcoord,formeractplacecell); # how much the former weights estimate my motion         
                                                deltay=dot(weightsycoord,actplacecell)-dot(weightsycoord,formeractplacecell);          
                                                weightsxcoord=weightsxcoord.+LRxcoord.*(-parameters[:dt].*parameters[:speed].*dir[1]+deltax).*eligibilitytrace;     
                                                weightsycoord=weightsycoord.+LRycoord.*(-parameters[:dt].*parameters[:speed].*dir[2]+deltay).*eligibilitytrace;                   
                                            end
                                        end
                                    end
                                 ########## ##########  END TRIAL ########## ##########             
                                    push!(historyX,currentposition[1]); # Store the last position visited 
                                    push!(historyY,currentposition[2]);
                                    ############### SAVING THE THINGS IN THE DIFFERENT CLASS ################
                                    ## in creating a new trial type one should write Trial(Trajectory, latency, searchpreference, actionmap) # action map atm is just z, then it will be improved adding a new attribute being value map 
                                    if searchinzones==0 # for searchpref 
                                        currenttrial=(trajectory=hcat(historyX,historyY),latency=t,actionmap=actorweights,valuemap=criticweights,TDerror=TDerrors,Xcoord=weightsxcoord,Ycoord=weightsycoord,SearchPref=searchpref)#, meanprobaction9=mean(probaction9),historyxmotion=historyXmotion,historydeltax=historydeltaX,historyeligibility=historyeligibilitytrace,historyymotion=historyYmotion,historydeltay=historydeltaY); # Creating the current trial with all its fields
                                    else 
                                       currenttrial=(trajectory=hcat(historyX,historyY),latency=t,actionmap=actorweights,valuemap=criticweights,TDerror=TDerrors,Xcoord=weightsxcoord,Ycoord=weightsycoord,SearchPref=searchpref/searchinzones)#, meanprobaction9=mean(probaction9),historyxmotion=historyXmotion,historydeltax=historydeltaX,historyeligibility=historyeligibilitytrace,historyymotion=historyYmotion,historydeltay=historydeltaY); # Creating the current trial with all its fields
                                    end
                                    push!(currentday,currenttrial) # Storing it in the current day 
                                end # end let all trial  related stuff
                            end 
                            ########## ##########  END DAY ########## ##########
                            dayc=(day=currentday, platformposition=[xp, yp]);
                            push!(currentexperiment,dayc);                    
                        ##################################################     
                    end # end let platform related variables 
                end # end loop number of days  
                ########## ##########  END EXPERIMENT ########## ##########
            push!(experiment,currentexperiment); # Storing the current experiment 
            ##################################################     
        end # end let currentexperiments, weights 
    end 
    save(NameOfFile, "parameters",parameters,"features",featuresexperiment,"data",experiment);
end # end experiment 

