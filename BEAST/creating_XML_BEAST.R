#/*=====Generating XML file for BEAST2 BDSKY=================*/
install.packages("XML")
install.packages("ape")

library(XML)
library(ape)
#setwd("E:/ARATA/Documents/Research/C_run/Phylodynamic/FMD_phylodynamics_2018April/READ/F") #for my laptop
#or
wd <- "C:/BEAST_with_jre.v2.5.2.Windows/BEAST/Nexus"
setwd("C:/BEAST_with_jre.v2.5.2.Windows/BEAST/Nexus")
i = 1
#make Nexus name
Nexus <- list()
for(i in 1:100)
{
  Nexus[[i]] <- paste0("FMD_Nexus",i)
}

# Read Nexus file and extract seq name and dates
lf_r<-list.files(path=wd,full.names=F,pattern="nex")

data_e<-lapply(lf_r,function(i){
  read.nexus.data(i)
})


for(i in seq(data_e))
{
  tem_rows <- c()
  num<-gsub("\\D","",substring(lf_r[i],10,12)) #extract number in file name to make sure data name matches to the actual file name
  #seq <- data_e[[i]][1][[1]]
  assign(paste0("df",num),data_e[[i]])
  #assign(paste0("seq",num),seq)

}





#===========Loop begins====================================#
days_to_add <- 1900
run_length <- 50000000

ite <- 100
# XML STRING 
for(i in 1:ite)
{
doc = newXMLDoc()
beast = newXMLNode("beast", attrs=c(beautitemplate='Standard', beautistatus='',
namespace="beast.core:beast.core.parameter:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood",
required="BEAST v2.5.2:BDSKY v1.4.5",
version="2.5"
)
,doc = doc)
#print(doc)

#========Create nodes========================================================================================================*/
#########Inside XML 3 big parts: I.data, II.function and map, III.run

#========I. data======================================================================================================================================================#
#data
nexus_name = Nexus[[i]]
datanode = newXMLNode("data", attrs=c(id=nexus_name,name="alignment"), parent=beast)
isolated_date_char <- ""
tem<-get(paste0("df",i))
latest_sample <- 0
oldest_sample <- 5000
# sequence node within data
for(k in seq(length(names(tem))))
{  
 
  isolate_name <- names(tem)[k]
  isolate_date<-as.numeric(gsub("^(?:[^_]+_){2}([^_]+).*", "\\1", names(tem)[k]))+days_to_add
  if(isolate_date>latest_sample)
  {
    latest_sample <- isolate_date
  }
  if(isolate_date<oldest_sample)
  {
    oldest_sample <- isolate_date
  }
  isolated_date_char <- paste0(isolated_date_char,isolate_name,"=",isolate_date,",")
  #Add sequence node
  seqnode = newXMLNode("sequence", attrs=c(id=paste0("seq_",isolate_name),
                                           taxon= isolate_name,
                                           totalcount = "4",
                                           value=toupper(paste0(tem[[k]],collapse=""))
  ), parent=datanode)
  
}
#output into value data in trait clause
isolated_date_char <- substr(isolated_date_char,1,nchar(isolated_date_char)-1)
time_to_oldest <- ceiling(latest_sample-oldest_sample) 


  

#=========II.function and map==========================================================================================================================================#
#function
functionnode = newXMLNode("function", 
attrs=c(spec="beast.core.util.Slice", 
        id="samplingProportionSlice", 
        arg=paste0("@samplingProportion_BDSKY_Serial.t:",nexus_name), 
        index="1",
        count="1"),
      parent=beast)
#map
mapnode1 = newXMLNode("map", attrs=c(name="Uniform"), "beast.math.distributions.Uniform",parent=beast)
mapnode2 = newXMLNode("map", attrs=c(name="Exponential"), "beast.math.distributions.Exponential",parent=beast)
mapnode3 = newXMLNode("map", attrs=c(name="LogNormal"), "beast.math.distributions.LogNormalDistributionModel",parent=beast)
mapnode4 = newXMLNode("map", attrs=c(name="Normal"), "beast.math.distributions.Normal",parent=beast)
mapnode5 = newXMLNode("map", attrs=c(name="Beta"), "beast.math.distributions.Beta",parent=beast)
mapnode6 = newXMLNode("map", attrs=c(name="LaplaceDistribution"), "beast.math.distributions.LaplaceDistribution",parent=beast)
mapnode7 = newXMLNode("map", attrs=c(name="prior"), "beast.math.distributions.Prior",parent=beast)
mapnode8 = newXMLNode("map", attrs=c(name="InverseGamma"), "beast.math.distributions.InverseGamma",parent=beast)
mapnode9 = newXMLNode("map", attrs=c(name="OneOnX"), "beast.math.distributions.OneOnX",parent=beast)

#========III. run========================================================================================================================================================#
#A. Creating children of run node===================================================================================#
#run node
runnode = newXMLNode("run", attrs=c(id="mcmc",spec="MCMC",chainLength="50000000", preBurnin="0"), parent=beast)

#####Inside run, A1.state, A2.init, A3.Posterior, A4.Operator, A5.Logger exist================================

#=============A1 State node=================================================================================#
#Creating children of state node========#
###state node - children nodes of run node
statenode = newXMLNode("state", attrs=c(id="state",storeEvery="5000"), parent=runnode)
 
######B. inside state node, B1. tree and B2. 5 parameters
#=========B1.tree node - children node of statenode
        treenode = newXMLNode("tree", attrs=c(id=paste0("Tree.t:",nexus_name),name="stateNode"), parent=statenode)
          ##C1. trait node - children of treenode
          #######Create sampling date list feeding into traitnode
          traitnode = newXMLNode("trait", attrs=c(id=paste0("dateTrait.t:",nexus_name),
                                                  spec="beast.evolution.tree.TraitSet",
                                                  traitname="date",
                                                  units="day",
                                                  value= isolated_date_char
                                                  ), parent=treenode)
              #child of trait
              taxanode = newXMLNode("taxa",attrs=c(id=paste0("TaxonSet.",nexus_name),spec="TaxonSet"), parent=traitnode)
                #child of taxa
              alignment = newXMLNode("alignment",attrs= c(idref=nexus_name), parent=taxanode)
                
          ##C2. taxonset - child of tree
          taxonset = newXMLNode("taxonset", attrs=c(idref=paste0("TaxonSet.",nexus_name)), parent=treenode)
#=========B2. parameter node
          parameter1 = newXMLNode("parameter", attrs=c(id=paste0("clockRate.c:",nexus_name),
                                                       name = "stateNode"),
                                                       "3.0E-5",parent=statenode)
          parameter2 = newXMLNode("parameter", attrs=c(id=paste0("origin_BDSKY_Serial.t:",nexus_name),
                                                       lower="0.0",
                                                       name="stateNode",
                                                       upper="Infinity"
                                                       ),
                                  "1500.0",parent=statenode)
          parameter3 = newXMLNode("parameter", attrs=c(id=paste0("becomeUninfectiousRate_BDSKY_Serial.t:",nexus_name),
                                                       lower="0.0",
                                                       name="stateNode",
                                                       upper="Infinity"
                                                          ),
                                                        "1.0",parent=statenode)
          parameter4 = newXMLNode("parameter", attrs=c(id=paste0("reproductiveNumber_BDSKY_Serial.t:",nexus_name),
                                                       dimension="3",
                                                       lower="0.0",
                                                       name="stateNode",
                                                       upper="Infinity"
          ),
          "2.0",parent=statenode)
          parameter5 = newXMLNode("parameter", attrs=c(id=paste0("samplingProportion_BDSKY_Serial.t:",nexus_name),
                                                       dimension="2",
                                                       lower="0.0",
                                                       name="stateNode",
                                                       upper="1.0"
          ),
          "0 1.0E-5",parent=statenode)
          
#=========A1 node done========#         
          
#========A2 init======================================================================================#
#create init inside run
          initnode = newXMLNode("init", attrs=c(id=paste0("RandomTree.t:",nexus_name),
                                                spec="beast.evolution.tree.RandomTree",
                                                estimate="false",
                                                initial=paste0("@Tree.t:",nexus_name),
                                                taxa=paste0("@",nexus_name)
                                                )
                                                , parent=runnode)
          #populationModel nested in init
          populationModel = newXMLNode("populationModel", attrs=c(id=paste0("ConstantPopulation0.t:",nexus_name),
                                                                  spec="ConstantPopulation"
                                                                  ), parent=initnode)
                     # parameter node nested in populationModel
                      parameter = newXMLNode("parameter", attrs=c(id=paste0("randomPopSize.t:",nexus_name),
                                                                  name = "popSize"
                                                                  ),"1.0",parent=populationModel)
#========A2 done===========#
                      
#============A3 Creating Posterior=========================================================================================#
# Two nested nodes, 1.prior, 2.likelihood
                      
posterior = newXMLNode("distribution", attrs=c(id="posterior", spec="util.CompoundDistribution"), parent = runnode)

#=======1. Setting up prior node within posterior==========================================================
  #prior nested in posterior
    prior = newXMLNode("distribution", attrs=c(id="prior", spec="util.CompoundDistribution"), parent = posterior)
######## 1.1 Nodes specify what types of distributions exist
     # Parameter distribution nested in prior
      para_dis = newXMLNode("distribution", attrs=c(id=paste0("BDSKY_Serial.t:",nexus_name),
                                                    spec = "beast.evolution.speciation.BirthDeathSkylineModel",
                                                    becomeUninfectiousRate=paste0("@becomeUninfectiousRate_BDSKY_Serial.t:",nexus_name),
                                                    origin=paste0("@origin_BDSKY_Serial.t:",nexus_name), 
                                                    reproductiveNumber=paste0("@reproductiveNumber_BDSKY_Serial.t:",nexus_name),
                                                    samplingProportion=paste0("@samplingProportion_BDSKY_Serial.t:",nexus_name),
                                                    tree=paste0("@Tree.t:",nexus_name)
                                                    ), parent = prior)
            # samplingRateChangeTimes
            samplingRateChangeTimes = newXMLNode("samplingRateChangeTimes", attrs=c(spec = "RealParameter",
                                                                                    value = paste0("0 ",time_to_oldest+0.0000001)
                                                                                    ), parent = para_dis)
            # reverseTimeArrays
            reverseTimeArrays = newXMLNode("reverseTimeArrays", attrs=c(spec ="BooleanParameter",
                                                                        value ="false false true false false false"),
                                           parent = para_dis)
########1.2 Nodes to specify priors for each distribution                                                    
# ========1.2.1 becomeUninfectiousRatePrior_BDSKY_Serial nested in prior
    becomeUninfectiousRatePrior_BDSKY_Serial =  newXMLNode("prior", attrs=c(id =paste0("becomeUninfectiousRatePrior_BDSKY_Serial.t:",nexus_name),
                                                                            name="distribution",
                                                                            x=paste0("@becomeUninfectiousRate_BDSKY_Serial.t:",nexus_name)
                                                                            ),
                                                                      parent = prior)
    
          # Distribution of becomeUninfectious
    becomeUninfectiousRatePrior_BDSKY_Serial_dis = newXMLNode("LogNormal", attrs=c(id ="LogNormalDistributionModel.0",
                                                                                   name="distr"), parent = becomeUninfectiousRatePrior_BDSKY_Serial)
          # parameter of this distribution
            parameter_dis1 = newXMLNode("parameter", attrs=c(id ="RealParameter.3",
                                                             estimate="false",
                                                             name="M"),"-5.1",parent = becomeUninfectiousRatePrior_BDSKY_Serial_dis)
            parameter_dis2 = newXMLNode("parameter", attrs=c(id ="RealParameter.4",
                                                             estimate="false",
                                                             name="S"),"0.7",parent = becomeUninfectiousRatePrior_BDSKY_Serial_dis)
                                                                              
  
#========1.2.2 Clockrate
  ClockPrior =  newXMLNode("prior",attrs=c(id =paste0("ClockPrior.c:",nexus_name),
                                   name="distribution",
                                    x=paste0("@clockRate.c:",nexus_name)),parent = prior)
            # prior for clock
            prior_clock = newXMLNode("Uniform",attrs=c(id ="Uniform.0",lower="2.0E-5",name="distr",upper="5.0E-5"), parent = ClockPrior)
                                                       
#=======1.2.3 origin
      originPrior =  newXMLNode("prior",attrs=c(id =paste0("originPrior_BDSKY_Serial.t:",nexus_name),
                                     name="distribution",
                                     x=paste0("@origin_BDSKY_Serial.t:",nexus_name)),parent = prior)
            # prior for origin
            prior_origin = newXMLNode("Uniform",attrs=c(id ="Uniform.3",name="distr",upper="2500.0"), parent = originPrior)     
   
#=======1.2.4 reproductiveNumberPrior_BDSKY_Serial
          
reproductiveNumberPrior_BDSKY_Serial =  newXMLNode("prior", attrs=c(id =paste0("reproductiveNumberPrior_BDSKY_Serial.t:",nexus_name),
                                                                                    name="distribution",
                                                                                    x=paste0("@reproductiveNumber_BDSKY_Serial.t:",nexus_name)
            ),
            parent = prior)
            
            # Distribution of reproductiveNumberPrior_BDSKY_Serial
            reproductiveNumberPrior_BDSKY_Serial_dis = newXMLNode("LogNormal", attrs=c(id ="LogNormalDistributionModel.1",
                                                                                           name="distr"), parent = reproductiveNumberPrior_BDSKY_Serial)
            # parameter of this distribution
            parameter_dis3 = newXMLNode("parameter", attrs=c(id ="RealParameter.5",
                                                             estimate="false",
                                                             name="M"),"0.0",parent = reproductiveNumberPrior_BDSKY_Serial_dis)
            parameter_dis4 = newXMLNode("parameter", attrs=c(id ="RealParameter.6",
                                                             estimate="false",
                                                             name="S"),"1.25",parent = reproductiveNumberPrior_BDSKY_Serial_dis)            
            
#=====1.2.5 samplingProportionPrior_BDSKY_Serial
            
samplingProportionPrior_BDSKY_Serial =  newXMLNode("prior", attrs=c(id =paste0("samplingProportionPrior_BDSKY_Serial.t:",nexus_name),
                                                                                name="distribution",
                                                                                x="@samplingProportionSlice"),            parent = prior)     
    # samplingProportionPrior_BDSKY_Serial_dis
    samplingProportionPrior_BDSKY_Serial =  newXMLNode("Uniform",attrs=c(id ="Uniform.4",lower="1.0E-5",name="distr",upper="0.05"),parent = samplingProportionPrior_BDSKY_Serial)
            
#==========1. Prior Done============================================================================#
    
#======2. Likelihood=================================================================================#
#####2.1 Set up likelihood
likelihood = newXMLNode("distribution", attrs=c(id="likelihood", spec="util.CompoundDistribution",useThreads="true"), parent = posterior)   
   
#####2.2 Setup tree likelihood
treelikelihood = newXMLNode("distribution", attrs=c(id=paste0("treeLikelihood.",nexus_name), 
                                                    spec="ThreadedTreeLikelihood",
                                                    data=paste0("@",nexus_name),
                                                    tree=paste0("@Tree.t:",nexus_name)
                                                    ), parent = likelihood) 
######2.3 Inside tree likelihood, set up siteModel and branchRateModel
###########2.3.1 site model
siteModel = newXMLNode("siteModel", attrs=c(id=paste0("SiteModel.s:",nexus_name), spec="SiteModel"), parent = treelikelihood)   

###########2.3.1.1 For site model, specify 4 parameters
site_parameter1 = newXMLNode("parameter", attrs=c(id=paste0("mutationRate.s:",nexus_name), estimate="false", name="mutationRate"), "1.0",parent = siteModel)   
site_parameter2 = newXMLNode("parameter", attrs=c(id=paste0("gammaShape.s:",nexus_name), estimate="false", name="shape"), "1.0",parent = siteModel)   
site_parameter3 = newXMLNode("parameter", attrs=c(id=paste0("proportionInvariant.s:",nexus_name), estimate="false",lower="0.0", name="proportionInvariant",upper="1.0"), "0.0",parent = siteModel)   
site_parameter4 = newXMLNode("substModel", attrs=c(id=paste0("JC69.s:",nexus_name), spec="JukesCantor"),parent = siteModel)   

###########2.3.2 branchRateModel
branchRateModel = newXMLNode("branchRateModel", attrs=c(id=paste0("StrictClock.c:",nexus_name), spec="beast.evolution.branchratemodel.StrictClockModel",
                                                        clock.rate=paste0("@clockRate.c:",nexus_name)
                                                      ),parent = treelikelihood)
#=========A3 Posterior done========#


#======A4 14 operators===============================================================================#
#A4.1 StrictClockRateScaler.c
StrictClockRateScaler.c = newXMLNode("operator", attrs=c(id=paste0("StrictClockRateScaler.c:",nexus_name),
 spec="ScaleOperator",
 parameter=paste0("@clockRate.c:",nexus_name),
 scaleFactor="0.75",
 weight="3.0"), parent = runnode)

#A4.2 StrictClockRateScaler.c
strictClockUpDownOperator.c = newXMLNode("operator", attrs=c(id=paste0("strictClockUpDownOperator.c:",nexus_name),
                                                         spec="UpDownOperator",
                                                        scaleFactor="0.75",
                                                         weight="3.0"), parent = runnode)
    # up and down
    up = newXMLNode("up", attrs=c(idref=paste0("clockRate.c:",nexus_name)),parent = strictClockUpDownOperator.c)
    down = newXMLNode("down", attrs=c(idref=paste0("Tree.t:",nexus_name)),parent = strictClockUpDownOperator.c)
                                     

#A4.3 BDSKY_SerialTreeScaler.t
BDSKY_SerialTreeScaler.t = newXMLNode("operator", attrs=c(id=paste0("BDSKY_SerialTreeScaler.t:",nexus_name),
                                                             spec="ScaleOperator",
                                                             scaleFactor="0.5",
                                                          tree=paste0("@Tree.t:",nexus_name),
                                                             weight="3.0"), parent = runnode)

#A4.4 BDSKY_SerialTreeRootScaler.t
BDSKY_SerialTreeRootScaler.t= newXMLNode("operator", attrs=c(id=paste0("BDSKY_SerialTreeRootScaler.t:",nexus_name),
                                                             spec="ScaleOperator",
                                                             rootOnly="true",
                                                             scaleFactor="0.5",
                                                             tree=paste0("@Tree.t:",nexus_name),
                                                             weight="3.0"), parent = runnode)

#A4.5 BDSKY_SerialUniformOperator.t
BDSKY_SerialUniformOperator.t = newXMLNode("operator", attrs=c(id=paste0("BDSKY_SerialUniformOperator.t:",nexus_name),
                                                               spec="Uniform",
                                                               tree=paste0("@Tree.t:",nexus_name),
                                                               weight="30.0"), parent = runnode)


#A4.6 BDSKY_SerialSubtreeSlide.t
BDSKY_SerialSubtreeSlide.t = newXMLNode("operator", attrs=c(id=paste0("BDSKY_SerialSubtreeSlide.t:",nexus_name),
                                                            spec="SubtreeSlide",
                                                            tree=paste0("@Tree.t:",nexus_name),
                                                            weight="15.0"), parent = runnode)
  
#A4.7 BDSKY_SerialNarrow.t
BDSKY_SerialNarrow.t = newXMLNode("operator", attrs=c(id=paste0("BDSKY_SerialNarrow.t:",nexus_name),
                                                     spec="Exchange",
                                                     tree=paste0("@Tree.t:",nexus_name),
                                                     weight="15.0"), parent = runnode)
  
#A4.8 BDSKY_SerialWide.t
BDSKY_SerialWide.t = newXMLNode("operator", attrs=c(id=paste0("BDSKY_SerialWide.t:",nexus_name),
                                                    spec="Exchange",
                                                    isNarrow="false",
                                                    tree=paste0("@Tree.t:",nexus_name),
                                                    weight="3.0"), parent = runnode)


#A4.9 BDSKY_SerialWilsonBalding.t
BDSKY_SerialWilsonBalding.t = newXMLNode("operator", attrs=c(id=paste0("BDSKY_SerialWilsonBalding.t:",nexus_name),
                                                            spec="WilsonBalding",
                                                            tree=paste0("@Tree.t:",nexus_name),
                                                            weight="3.0"), parent = runnode)
  
#A4.10 becomeUninfectiousRateScaler_BDSKY_Serial.t
becomeUninfectiousRateScaler_BDSKY_Serial.t =newXMLNode("operator", attrs=c(id=paste0("becomeUninfectiousRateScaler_BDSKY_Serial.t:",nexus_name),
                                                                            spec="ScaleOperator",
                                                                            parameter=paste0("@becomeUninfectiousRate_BDSKY_Serial.t:",nexus_name),
                                                                            scaleFactor="0.75",
                                                                            weight="2.0"), parent = runnode)
  
#A4.11 reproductiveNumberScaler_BDSKY_Serial.t
reproductiveNumberScaler_BDSKY_Serial.t =newXMLNode("operator", attrs=c(id=paste0("reproductiveNumberScaler_BDSKY_Serial.t:",nexus_name),
                                                                        spec="ScaleOperator",
                                                                        parameter=paste0("@reproductiveNumber_BDSKY_Serial.t:",nexus_name),
                                                                        scaleFactor="0.75",
                                                                        weight="10.0"), parent = runnode)

#A4.12 samplingProportionScaler_BDSKY_Serial.t
samplingProportionScaler_BDSKY_Serial.t = newXMLNode("operator", attrs=c(id=paste0("samplingProportionScaler_BDSKY_Serial.t:",nexus_name),
                                                                         spec="ScaleOperator",
                                                                         parameter=paste0("@samplingProportion_BDSKY_Serial.t:",nexus_name),
                                                                         scaleFactor="0.75",
                                                                         weight="10.0"), parent = runnode)

#A4.13 updownBD_BDSKY_Serial.t
updownBD_BDSKY_Serial.t = newXMLNode("operator", attrs=c(id=paste0("updownBD_BDSKY_Serial.t:",nexus_name),
                                                         spec="UpDownOperator",
                                                         scaleFactor="0.75",
                                                         weight="2.0"), parent = runnode)

     # up and down
        up = newXMLNode("up", attrs=c(idref=paste0("reproductiveNumber_BDSKY_Serial.t:",nexus_name)),parent = updownBD_BDSKY_Serial.t)
        down = newXMLNode("down", attrs=c(idref=paste0("becomeUninfectiousRate_BDSKY_Serial.t:",nexus_name)),parent = updownBD_BDSKY_Serial.t)
      
#A4.14 origScaler_BDSKY_Serial.t
origScaler_BDSKY_Serial.t   =    newXMLNode("operator", attrs=c(id=paste0("origScaler_BDSKY_Serial.t:",nexus_name),
                                                                spec="ScaleOperator",
                                                                parameter=paste0("@origin_BDSKY_Serial.t:",nexus_name),
                                                                scaleFactor="0.75",
                                                                weight="1.0"), parent = runnode)  
        
#======A4 done======#

#======A5 Logger=======================================================================================#
# Three loggers
#A5.1 tracelog
tracelog   =    newXMLNode("logger", attrs=c(id="tracelog",
                                             fileName=paste0(nexus_name,".log"),
                                             logEvery="1000",
                                             model="@posterior",
                                             sanitiseHeaders="true",
                                             sort="smart"), parent = runnode) 
  #Inside tracelog 11 logs
  log1 = newXMLNode("log", attrs=c(idref="posterior"),parent = tracelog)
  log2 = newXMLNode("log", attrs=c(idref="likelihood"),parent = tracelog)
  log3 = newXMLNode("log", attrs=c(idref="prior"),parent = tracelog)
  log4 = newXMLNode("log", attrs=c(idref=paste0("treeLikelihood.",nexus_name)),parent = tracelog)
  log5 = newXMLNode("log", attrs=c(id=paste0("TreeHeight.t:",nexus_name),
                                   spec="beast.evolution.tree.TreeHeightLogger",
                                   tree=paste0("@Tree.t:",nexus_name)
                                   ),parent = tracelog)
 
  log6 = newXMLNode("log", attrs=c(idref=paste0("clockRate.c:",nexus_name)),parent = tracelog)
  log7 = newXMLNode("log", attrs=c(idref=paste0("BDSKY_Serial.t:",nexus_name)),parent = tracelog)
  log8 = newXMLNode("log", attrs=c(idref=paste0("origin_BDSKY_Serial.t:",nexus_name)),parent = tracelog)
  log9 = newXMLNode("log", attrs=c(idref=paste0("becomeUninfectiousRate_BDSKY_Serial.t:",nexus_name)),parent = tracelog)
  log10 = newXMLNode("log", attrs=c(idref=paste0("reproductiveNumber_BDSKY_Serial.t:",nexus_name)),parent = tracelog)
  log11 = newXMLNode("log", attrs=c(idref=paste0("samplingProportion_BDSKY_Serial.t:",nexus_name)),parent = tracelog)
  
#A5.2 screenlog
  screenlog   =    newXMLNode("logger", attrs=c(id="screenlog",logEvery="1000"),parent = runnode) 
  # inside 3 logs
  log1_1 = newXMLNode("log", attrs=c(idref="posterior"),parent = screenlog)
  log1_2 = newXMLNode("log", attrs=c(idref="likelihood"),parent = screenlog)
  log1_3 = newXMLNode("log", attrs=c(idref="prior"),parent = screenlog)    
  
#A5.3 treelog.t
  treelog.t   =    newXMLNode("logger", attrs=c(id=paste0("treelog.t:",nexus_name),
                                                fileName="$(tree).trees",
                                                logEvery="1000",
                                                mode ="tree"),parent = runnode) 
  # inside 1 log
  log1_1_1 = newXMLNode("log", attrs=c(id=paste0("TreeWithMetaDataLogger.t:",nexus_name),
                        spec="beast.evolution.tree.TreeWithMetaDataLogger",
                        tree=paste0("@Tree.t:",nexus_name)
                          ),
                        parent = treelog.t)
 

#======A5 done======#
#========Create node done================*/
  saveXML(doc, paste0("C:/BEAST_with_jre.v2.5.2.Windows/BEAST/xml2/nexus",i,".xml"), prefix='<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
}
#print(doc)
#Export XML file


#======Preprare batch file to process xml into BEAST
k = 0
for(i in k:19)
{
  
  fileConn<-file(paste0("C:\\BEAST_with_jre.v2.5.2.Windows\\BEAST\\run1_",i*5+1,"_",i*5+5,".bat"))
  writeLines(c("c:",
                "cd C:\\BEAST_with_jre.v2.5.2.Windows\\BEAST",
               paste0("java -Dbeast.load.jars=true -jar lib\\beast.jar -overwrite xml2\\nexus",i*5+1,".xml"),
               paste0("java -Dbeast.load.jars=true -jar lib\\beast.jar -overwrite xml2\\nexus",i*5+2,".xml"),
               paste0("java -Dbeast.load.jars=true -jar lib\\beast.jar -overwrite xml2\\nexus",i*5+3,".xml"),
               paste0("java -Dbeast.load.jars=true -jar lib\\beast.jar -overwrite xml2\\nexus",i*5+4,".xml"),
               paste0("java -Dbeast.load.jars=true -jar lib\\beast.jar -overwrite xml2\\nexus",i*5+5,".xml")
               ), fileConn)
  close(fileConn)
  
}

k = 0
for(i in k:9)
{
  
  fileConn<-file(paste0("C:\\BEAST_with_jre.v2.5.2.Windows\\BEAST\\run1_",i*10+1,"_",i*10+10,".bat"))
  writeLines(c("c:",
               "cd C:\\BEAST_with_jre.v2.5.2.Windows\\BEAST",
               paste0("java -Dbeast.load.jars=true -jar lib\\beast.jar -overwrite xml2\\nexus",i*10+1,".xml"),
               paste0("java -Dbeast.load.jars=true -jar lib\\beast.jar -overwrite xml2\\nexus",i*10+2,".xml"),
               paste0("java -Dbeast.load.jars=true -jar lib\\beast.jar -overwrite xml2\\nexus",i*10+3,".xml"),
               paste0("java -Dbeast.load.jars=true -jar lib\\beast.jar -overwrite xml2\\nexus",i*10+4,".xml"),
               paste0("java -Dbeast.load.jars=true -jar lib\\beast.jar -overwrite xml2\\nexus",i*10+5,".xml"),
               paste0("java -Dbeast.load.jars=true -jar lib\\beast.jar -overwrite xml2\\nexus",i*10+6,".xml"),
               paste0("java -Dbeast.load.jars=true -jar lib\\beast.jar -overwrite xml2\\nexus",i*10+7,".xml"),
               paste0("java -Dbeast.load.jars=true -jar lib\\beast.jar -overwrite xml2\\nexus",i*10+8,".xml"),
               paste0("java -Dbeast.load.jars=true -jar lib\\beast.jar -overwrite xml2\\nexus",i*10+9,".xml"),
               paste0("java -Dbeast.load.jars=true -jar lib\\beast.jar -overwrite xml2\\nexus",i*10+10,".xml")
                 ), fileConn)
  close(fileConn)
  
}
