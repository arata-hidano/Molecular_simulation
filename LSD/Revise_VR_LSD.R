#=========================================================================#
# Converting NEXUS to Phylip format
#=========================================================================#

#1. Install packages
install.packages("phangorn")
install.packages("ape")

library(phangorn); library(ape)
sessionInfo()
setwd("E:/ARATA/Documents/Research/C_run/Phylodynamic/FMD_phylodynamics_2018April/READ/F")

#2. Read data
lf_r<-list.files(path="E:/ARATA/Documents/Research/C_run/Phylodynamic/FMD_phylodynamics_2018April/READ/F",full.names=F,pattern="nex")


data_e<-lapply(lf_r,function(i){
  read.nexus.data(i)
})

k <- 0
for(i in seq(data_e))
{
  num<-gsub("\\D","",substring(lf_r[i],10,12))
  data_e[[i]] <- phyDat(data_e[[i]], type = "DNA")
  assign(paste0("df",num),data_e[[i]])
  
}
#3. Run modeltest
#mt<-modelTest(dat, tree = NULL, model = "all")
#bestmodel <- mt$Model[which.min(mt$AICc)]

trees <- list()
#4. Make a NJ tree using ape's NJ function then ML tree using pml function
# check http://www.molecularecologist.com/2016/02/quick-and-dirty-tree-building-in-r/

for(i in 1:100)
{
  tem<-get(paste0("df",i))
  tem1 <- dist.ml(tem)
  treeNJ <- NJ(tem1)
  fit <- pml(treeNJ, tem)
  fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")
  trees <- c(trees,fitJC)
}
#plot(fitJC)

#5. Export the tree in Newick format
for(i in 1:100)
{
    exp_tree <- trees[[17+22*(i-1)]]
    num<-gsub("\\D","",substring(lf_r[i],10,12))
  write.tree(exp_tree,file=paste0("E:/ARATA/Documents/Research/C_run/Phylodynamic/FMD_phylodynamics_2018April/READ/F/tree_file/tree",num,".tre"))
  #needs to create a text file that contains number of leaves in the first line then name of the tip with sampling date
  tem<-get(paste0("df",i))
  tem_rows <- c(length(names(tem)),"")
  for(k in seq(length(names(tem))))
  {
    new_rows <- c(names(tem)[k],gsub("^(?:[^_]+_){2}([^_]+).*", "\\1", names(tem)[k]))
    tem_rows <- rbind(tem_rows,new_rows)
    
  }
  tem_rows <- as.data.frame(tem_rows)
  write.table(tem_rows, paste0("E:/ARATA/Documents/Research/C_run/Phylodynamic/FMD_phylodynamics_2018April/READ/F/date_file/dates",num,".txt"), sep="\t",col.names =F,row.names = F,quote=FALSE)
}

#6. Run command line from R environment using system2

##try lsd3 beta
setwd("E:/lsd-0.3beta-master/lsd-0.3beta-master/src")
#system2('lsd', args = c('-i',paste0('E:/ARATA/Documents/Research/C_run/Phylodynamic/FMD_phylodynamics_2018April/READ/F/tree_file/tree',"7",'.tre'),'-d',paste0('E:/ARATA/Documents/Research/C_run/Phylodynamic/FMD_phylodynamics_2018April/READ/F/date_file/dates',"7",'.txt'), '-c', '-r','a'))

for(i in 1:100)
{
  system2('lsd', args = c('-i',paste0('E:/ARATA/Documents/Research/C_run/Phylodynamic/FMD_phylodynamics_2018April/READ/F/tree_file/tree',i,'.tre'),'-d',paste0('E:/ARATA/Documents/Research/C_run/Phylodynamic/FMD_phylodynamics_2018April/READ/F/date_file/dates',i,'.txt'), '-c', '-r','a'))
  
}


#7. Now read result file and extract tMRCA value
setwd("E:/ARATA/Documents/Research/C_run/Phylodynamic/FMD_phylodynamics_2018April/READ/F/tree_file")
time_list <- c()
for(i in 1:100)
{
  result_file <- read.delim(paste0("tree",i,".tre.result"),header=F,sep = "\t")
  ##row 17 is always containing tMRCA
  whole_string <-as.character(result_file[20,1])
  time <- as.numeric(sub(" ,.+","",sub("^.+tMRCA ", "",whole_string)))
  time_list <- cbind(time_list,time)
}

hist(time_list,breaks=50)
dat_time <- as.data.frame(t(time_list))
write.table(dat_time,"E:/ARATA/Documents/Research/C_run/Phylodynamic/FMD_phylodynamics_2018April/time_LSD.csv",sep=",",col.names = F,row.names = F)
#8. Make it to ggplot 
#Do it in Remaking_figures_18Feb2019.R 
