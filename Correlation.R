# Correlation Plot
final.table
library(ggplot2)
tibble::glimpse(final.table)
final.table.run1.run2 <- final.table[-which(final.table$data == "OCT."),]
ggplot(final.table.run1.run2, aes(x = reads, y = n_segments, color = data)) +
  geom_point() +
  geom_abline(intercept = 8.66, slope = 3.517e-05)
cor(final.table.run1.run2$reads, final.table.run1.run2$n_segments)
# 0.614463
# spearmanss
cor.test(final.table.run1.run2$reads, final.table.run1.run2$n_segments, method = 'spearman')
cor.test(final.table.run1.run2$reads, final.table.run1.run2$n_segments)
# Pearson's product-moment correlation
# 
# data:  final.table.run1.run2$reads and final.table.run1.run2$n_segments
# t = 5.5072, df = 50, p-value = 1.263e-06
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.4103777 0.7599324
# sample estimates:
#      cor 
# 0.614463 

# regression line
require(stats)
reg <- lm(n_segments ~ reads, data = final.table.run1.run2)
reg
# Call:
#   lm(formula = n_segments ~ reads, data = final.table.run1.run2)
# 
# Coefficients:
#   (Intercept)        reads  
# 8.666e+00    3.517e-05  
