#Working with my CNV data
setwd('/home/joyanta/Documents/Projects/CNV/test/Group1/simul1')
library(DNAcopy)
#Load DNA Pos Coverage data
path_ref<- file.path("/home/joyanta/Documents/Projects/CNV/test/Group1/simul1","readsInPerPositionGenome_ref.txt")
path_ref

path_test<- file.path("/home/joyanta/Documents/Projects/CNV/test/Group1/simul1","readsInPerPositionGenome_test.txt")
path_test

ref<- read.table(path_ref, 
                     header = FALSE, sep = '\t',
                     stringsAsFactors = FALSE)
ref<- ref[,1]

test<- read.table(path_test, 
                 header = FALSE, sep = '\t',
                 stringsAsFactors = FALSE)
test<-test[,1]


#Returns the Mean converages. Each data represents one mean value for window length positions
slideFunct <- function(data, window, step){
  total <- length(data)
  spots <- seq(from=1, to=total, by=step)
  #print(spots)
  result <- vector(length = length(spots))
  sz = length(spots)
  for(i in 1:sz){
    if(spots[i]+window-1<=total){
      result[i] <- mean(data[spots[i]:(spots[i]+window-1)])
    }else{
      result[i] <- mean(data[spots[i]:total])  
    }
  }
  return(result)
}


# log(1) = 0, (ref, ref=0,1) replaces ref val 0 with 1
my_data_vec<- log2(1+(test/replace(ref,ref==0,1))) #Actual formula
#my_data_vec<- log2(1+test)
class(my_data_vec)
length(my_data_vec)
#############PARAMETERS###############
windowsize = 50

my_data_vec<- slideFunct(my_data_vec,windowsize,windowsize)
############################

maplocation<-(1:length(my_data_vec))
length(maplocation)

#As API can work on multiple Chromosome, chomosome number needs to be specified
chrom<-rep(22,length(my_data_vec))
length(chrom)

#DNAcopy API
my_seg_input<- CNA(my_data_vec, chrom, maplocation, data.type="logratio")
#print("Checking Input Class")
#class(my_seg_input)

smooth.my_seg_input<- smooth.CNA(my_seg_input)
length((smooth.my_seg_input[,1]))

#?segment
start.time<-Sys.time()
my_seg_output<-segment(smooth.my_seg_input,alpha= 0.0005, undo.splits="sdundo", 
                       undo.SD=3, verbose=1)

end.time<-Sys.time()
duration<-end.time-start.time
duration



my_seg_output2<- my_seg_output

#print("Checking Output Class")
#class(my_seg_output)
#class(my_seg_output$output)

length(my_seg_output2$output[,3])
my_seg_output2$output[,3:5]<-my_seg_output2$output[,3:5]*windowsize
my_seg_output2$output[,3][1]<-1 #As previously multiplied by windowsize
my_seg_output2$output[,4][length(my_seg_output2$output[,4])]<-length(test) #Correct last Length
temp<-my_seg_output2$output[,4][length(my_seg_output2$output[,4])]-my_seg_output2$output[,3][length(my_seg_output2$output[,3])]
temp
my_seg_output2$output[,5][length(my_seg_output2$output[,5])]<-temp #Corrected Last region length
#?write.table

####### Merge small regions to neighbour ############
print("Merging small regions")
cutoff = 500
last_i = length(my_seg_output2$output[,5])
print(last_i)
rows_to_remove <- vector("list",last_i)
rr_counter= 0
for(i in 1:(last_i-1)){
  if(my_seg_output2$output[,5][i] < cutoff) {
    my_seg_output2$output[,5][i+1] = my_seg_output2$output[,5][i] + my_seg_output2$output[,5][i+1] #increase region size
    my_seg_output2$output[,3][i+1] = my_seg_output2$output[,3][i] #update start
    my_seg_output2$output[,6][i+1] = ((my_seg_output2$output[,6][i]*my_seg_output2$output[,5][i]) + (my_seg_output2$output[,6][i+1]*(my_seg_output2$output[,5][i+1] - my_seg_output2$output[,5][i])))/my_seg_output2$output[,5][i+1] #update mean
    rr_counter = rr_counter + 1
    rows_to_remove[[rr_counter]] <- i
  }
}

#for last row add to prev
if (my_seg_output2$output[,5][last_i] < cutoff) {
  my_seg_output2$output[,5][last_i-1] = my_seg_output2$output[,5][last_i] + my_seg_output2$output[,5][last_i-1] #increase region size
  my_seg_output2$output[,4][last_i-1] = my_seg_output2$output[,4][last_i] #update end
  my_seg_output2$output[,6][last_i-1] = ((my_seg_output2$output[,6][last_i]*my_seg_output2$output[,5][last_i]) + (my_seg_output2$output[,6][last_i-1]*(my_seg_output2$output[,5][last_i-1] - my_seg_output2$output[,5][last_i])))/my_seg_output2$output[,5][last_i-1] #update mean
  rr_counter = rr_counter + 1
  rows_to_remove[[rr_counter]] <- last_i
} 
remove_list <- rows_to_remove[!sapply(rows_to_remove[],is.null)]
remove_list <-as.numeric(remove_list[])
my_seg_output2$output <- my_seg_output2$output[-remove_list,]
print(nrow(my_seg_output2$output))


#######Save file ###########
print("Writing")
write.csv(my_seg_output2$output, file="/home/joyanta/Documents/Projects/CNV/test/Group1/simul1/segmented_result_Win_50_RCutoff_500_LR_Alpha=0.0005.csv", row.names = FALSE)

#plot(my_seg_output, plot.type="w")
print("DONE")
