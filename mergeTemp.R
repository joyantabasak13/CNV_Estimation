setwd('/home/joyanta/Desktop/copy-number-variation/Code')
load('data.Rdata')
####### Merge small regions to neighbour ############
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
write.csv(my_seg_output2$output, file="/home/joyanta/Desktop/copy-number-variation/Code/segmented_result_WithSmallReg_LR_Alpha=0.0005.csv")
a <-c(1,2,3)
class(a)
class(remove_list)
remove_list <-as.numeric(remove_list[])
class(remove_list)
print(nrow(my_seg_output2$output))
print(my_seg_output2$output[1,])
my_seg_output3 <- my_seg_output2
my_seg_output3$output <- my_seg_output3$output[-remove_list,]
print(nrow(my_seg_output3$output))

#my_seg_output2$output <- my_seg_output2$output[]
#print(nrow(my_seg_output2$output))

#my_seg_output2$output
write.csv(my_seg_output3$output, file="/home/joyanta/Desktop/copy-number-variation/Code/segmented_result_WithOutSmallReg_LR_Alpha=0.0005.csv", row.names = FALSE)

plot(my_seg_output, plot.type="w")