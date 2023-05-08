## test the TCA object constructor
d1<- data.frame(timepoint = rep(c("0h", "24h", "48h", "72h"), 2), group = rep(c(1, 2, 3, 4), 2))

d3 <- data.frame(sampleid = 1:8, timepoint = rep(c("0h", "24h", "48h", "72h"), 2),
                 group = rep(c(1, 2, 3, 4), 2))

gf1 <- data.frame(chr = c(rep("chr1", 3), rep("chr2", 2), rep("chr4", 2)),
                  start = seq(100, 2000, by = 300), end = seq(100, 2000, by = 300) + 150)

gf2 <- data.frame(CHR = c(rep("chr1", 3), rep("chr2", 2), rep("chr4", 2)),
                  start = seq(100, 2000, by = 300), end = seq(100, 2000, by = 300) + 150,
                  id = paste0("peak", 1:7))

gf3 <- data.frame(chr = c(rep("chr1", 3), rep("chr2", 2), rep("chr4", 2)),
                  start = seq(100, 2000, by = 300), end = seq(100, 2000, by = 300) + 150,
                  id = paste0("peak", 1:7))

tca <- TCA(design = d3, genomicFeature = gf3)
expect_error({
  tca <- TCA(design = d1, genomicFeature = gf3)
})
expect_error({
  tca <- TCA(design = d3, genomicFeature = gf1)
})
expect_warning({
  tca <- TCA(design = d3, genomicFeature = gf2)
})

c1 <- matrix(sample(500, 56), nrow = 7, dimnames = list(paste0("peak",
                                                               1:7), 1:8))
c2 <- matrix(sample(500, 48), nrow = 6, dimnames = list(paste0("peak",
                                                               1:6), 1:8))
c3 <- matrix(sample(500, 49), nrow = 7, dimnames = list(paste0("peak",
                                                               1:7), 1:7))
tca <- TCA(design = d3, counts = c1, genomicFeature = gf3)
expect_error({
  TCA(design = d3, counts = c2, genomicFeature = gf3)
})
expect_error({
  TCA(design = d3, counts = c3, genomicFeature = gf3)
})

## test the correctness of the merge result results
peaks <- data.frame(chr = c(rep("chr1",4),rep("chr2", 3), rep("chr3",2)),
                    start = c(100,148,230,300,330,480,1000,700,801),
                    end = c(150,220,500,450,600,900,1050,760,900))

merged_peaks <- peakreference(data = peaks, merge = T, overlap = 1)

peaks_expect <- data.frame(chr = c(rep("chr1",2),rep("chr2", 2), rep("chr3",2)),
                           start = c(100, 230, 330, 1000, 700, 801),
                           end = c(220, 500, 900, 1050, 760, 900),
                           id = paste0("peak", 1:6))

expect_equal(merged_peaks, peaks_expect)

merged_peaks2 <- peakreference(data = peaks, merge = T, ratio = 0.2)
peaks_expect2 <- data.frame(chr = c(rep("chr1",3),rep("chr2", 2), rep("chr3",2)),
                           start = c(100,148, 230, 330, 1000, 700, 801),
                           end = c(150, 220, 500, 900, 1050, 760, 900),
                           id = paste0("peak", 1:7))

expect_equal(merged_peaks2, peaks_expect2)

