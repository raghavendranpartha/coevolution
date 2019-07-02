require(precrec)
require(ggplot2)
s1 <- c(1, 3, 2, 4)
s2 <- c(5, 6, 7, 8)
s3 <- matrix(1:8, 4, 2)

scores1 <- join_scores(s1, s2)
scores2 <- join_scores(s1, s2, s3)

l1 <- c(1, 0, 1, 1)
l2 <- c(1, 0, 1, 1)
l3 <- c(1, 0, 1, 0)

labels1 <- join_labels(l1, l2)
labels2 <- join_labels(l1, l3)

mmmdat <- mmdata(scores1, labels1, modnames = c("mod1", "mod2"), dsids = c(1,2))
#mmmdat <- mmdata(scores1, labels1, dsids = c(1,2))
mmcurves <- evalmod(mmmdat,cb_alpha = 0.05)
autoplot(mmcurves, "PRC", show_cb = TRUE )