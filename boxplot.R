set.seed(1)
df <- data.frame(
  gp = factor(rep(letters[1:3], each = 10)),
  y = rnorm(30)
)
ds <- do.call(rbind, lapply(split(df, df$gp), function(d) {
  data.frame(mean = mean(d$y), sd = sd(d$y), gp = d$gp)
}))
ggplot(df, aes(gp, y)) +
  geom_point() +
  geom_point(data = ds, aes(y = mean), colour = 'red', size = 3)

ggplot(df) +
  geom_point(aes(gp, y)) +
  geom_point(data = ds, aes(gp, mean), colour = 'red', size = 3)

ggplot() +
  geom_point(data = df, aes(gp, y)) +
  geom_point(data = ds, aes(gp, mean), colour = 'red', size = 3) +
  geom_errorbar(
    data = ds,
    aes(gp, mean, ymin = mean - sd, ymax = mean + sd),
    colour = 'red',
    width = 0.4
  )

###########

library(ggplot2)
library(viridis)
library(hrbrthemes)
# ggplot() +
#   geom_boxplot(data = seq.r1.after, aes(x = sampleID, y = sequences)) +
#   geom_boxplot(data = seq.r1.after) +
#   geom_boxplot(data = seq.r2.after) +
#   scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
#   geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
#   theme_ipsum() +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(size=11)
#   ) +
#   ggtitle("Number of segments in Primary Vs. Brain Mets") +
#   xlab("")
# code above doesn't work, commenting out and starting new
final.table$data[which(final.table$data == "OCT.")] <- "Control"
library(tidyverse)
ggplot(final.table, aes(x=data, y=reads, fill = data)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9, position = position_jitter(0.1)) +
  labs(title = "Reads per cohort", x = "Cohort", y = "Reads(n)")

############changed mind, making a bar chart
seq.r1.before <- as_tibble(seq.r1.before)
seq.r2.before <- as_tibble(seq.r2.before)

colnames(seq.r1.before)[2] <- "sequences"
seq.r1.before <- seq.r1.before %>%
  add_column(status = "paired", .before = "sequences")
seq.r1.before <- seq.r1.before %>%
  add_column(time = "before", .after = "status")
seq.r1.before <- seq.r1.before %>%
  add_column(read = "r1", .after = "time")

colnames(seq.r2.before)[2] <- "sequences"
seq.r2.before <- seq.r2.before %>%
  add_column(status = "paired", .before = "sequences")
seq.r2.before <- seq.r2.before %>%
  add_column(time = "before", .after = "status")
seq.r2.before <- seq.r2.before %>%
  add_column(read = "r2", .after = "time")

seq.r1.after <- as_tibble(seq.r1.after)
seq.r2.after <- as_tibble(seq.r2.after)

seq.r1.after <- seq.r1.after %>%
  add_column(time = "after", .after = "status")+
  add_column(read = "r1", .after = "time")

seq.r2.after <- seq.r2.after %>%
  add_column(time = "after", .after = "status")
seq.r2.after <- seq.r2.after %>%
  add_column(read = "r2", .after = "time")

seq <- bind_rows(
  seq.r1.before,
  seq.r2.before,
  seq.r1.after,
  seq.r2.after)

seq

seq %>%
  mutate(time = factor(time, levels = c('before', 'after'))) %>%
  filter(status == "paired") %>%
  group_by(time, read) %>%
  mutate(sequences = mean(sequences)) %>%
  ggplot(aes(x = time, fill = read, y = sequences)) + 
  geom_col()

#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summarized
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

seq.sd <- data_summary(seq, varname = "")

colnames(seq.r1.before)[1] <- "sampleID"

seq.r1.before <- seq.r1.before %>%
  add_column()

list(
  r2.before = seq.r2.before,
  r2.after = seq.r2.after,
  
)

##############Practice barplot with wrror bars################
library(ggplot2)

data <- data.frame(
  name = letters[1:5],
  value = sample(seq(4,15),5),
  sd=c(1,0.2,3,2,4)
)

data <- as_tibble(data)
data
ggplot(data) +
  geom_bar( aes(x=name, y=value), stat = "identity", fill = "skyblue", alpha = 0.7) +
  geom_errorbar( aes(x = name, ymin = value-sd, ymax = value+sd), width = 0.4, colour = "orange", alpha = 0.9, linewidth = 1.3)
##################################
# new bar graph
# need vectors for mean sequences of before, after_rw, after_r2, and after_all
before
before <- seq[grepl("before", seq$time),]
mean.before <- mean(before$sequences)
sd.before <- sd(before$sequences)
before.vec <- data.frame(mean = mean.before, sd = sd.before)
rownames(before.vec)[1] <- "before"

after.r1 <- subset(seq, time == "after" & read == "r1")
mean.after.r1 <- mean(after.r1$sequences)
sd.after.r1 <- sd(after.r1$sequences)
after.r1.vec <- data.frame(mean = mean.after.r1, sd = sd.after.r1, row.names = "after.r1")
before.vec <- before.vec %>%
  add_row(after.r1.vec)
rownames(before.vec)[2] <- "after.r1"
total.mean.seq <- before.vec

after.r2 <- subset(seq, time =="after" & read == "r2")
mean.after.r2 <- mean(after.r2$sequences)
sd.after.r2 <- sd(after.r2$sequences)
after.r2.vec <- data.frame(mean = mean.after.r2, sd = sd.after.r2, row.names = "after.r2")
total.mean.seq <- total.mean.seq %>%
  add_row(after.r2.vec)
rownames(total.mean.seq)[3] <- "after.r2"

after.all <- subset(seq, time == "after" & status == "paired")
mean.after.all <- mean(after.all$sequences)
sd.after.all <- sd(after.all$sequences)
after.all.vec <- data.frame(mean = mean.after.all, sd = sd.after.all, row.names = "after.all")
total.mean.seq <- total.mean.seq %>%
  add_row(after.all.vec)
rownames(total.mean.seq)[4] <- "after.all"
total.mean.seq <- as_tibble(total.mean.seq)
row.seq <- c("before", "after.r1", "after.r2", "after.all")
total.mean.seq <- total.mean.seq %>%
  add_column(row.seq, .before = "mean")

#plot
sequences <- total.mean.seq$mean
sd <- total.mean.seq$sd

total.mean.seq.2 <- total.mean.seq
total.mean.seq.2$mean <- (total.mean.seq.2$mean) / 1000000
total.mean.seq.2$sd <- (total.mean.seq.2$sd) / 1000000
sequences <- sequences/1000000
sd <- sd/ 1000000
ggplot(total.mean.seq.2) +
  geom_bar( aes(x=row.seq, y=mean), stat = "identity", fill = "skyblue", alpha = 0.7) +
  geom_errorbar( aes(x = row.seq, ymin = sequences-sd, ymax = sequences+sd), width = 0.4, colour = "orange", alpha = 0.9, linewidth = 1.3) +
  scale_x_discrete(limits = c("before", "after.r1", "after.r2", "after.all")) +
  labs(title = "Number of sequences before and after trimming", x = "Timepoint", y = "Sequences (per million)")

pearson <- cor.test(final.table.run1.run2$reads, final.table.run1.run2$n_segments)
p.adjust(pearson, method = "fdr")
