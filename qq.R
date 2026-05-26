# Produceert de data achter qq.Plot()
# Neemt dezelfde parameter maar returns a dataframe with colums:
# Chrom   model   y   x
# The default plotting command is the following:
# ggplot(data = tmp, aes(x = .data$x, y = .data$y, colour = .data$model)) + 
#   facet_wrap(~Chrom) + 
#   geom_point() + 
#   theme_bw() + 
#   xlab(expression(paste("Expected -log"[10], "(p)", sep = ""))) + 
#   ylab(expression(paste("Observed -log"[10], "(p)", sep = ""))) + 
#   scale_colour_brewer(palette = "Set1") + 
#   geom_abline(slope = 1, intercept = 0, linetype = 2) + 
#   theme(text = element_text(size = 15))

# An example snippet
#qq(data, trait="TRAIT") %>%
#  filter(Chrom == c("chr03", "chr04", "chr05")) %>%
#    ggplot(data = tmp, aes(x = .data$x, y = .data$y, colour = .data$model)) + 
#    facet_wrap(~Chrom) + 
#    geom_point() + 
#    theme_bw() + 
#    xlab(expression(paste("Expected -log"[10], "(p)", sep = ""))) + 
#    ylab(expression(paste("Observed -log"[10], "(p)", sep = ""))) + 
#    scale_colour_brewer(palette = "Set1") + 
#    geom_abline(slope = 1, intercept = 0, linetype = 2) + 
#    theme(text = element_text(size = 15))

# Depents on
library(dplyr)

qq <- function(data, trait, model = NULL){
  stopifnot(inherits(data, "GWASpoly.fitted"))
  stopifnot(is.element(trait, names(data@scores)))
  all.models <- colnames(data@scores[[trait]])
  if (is.null(models)) {
    models <- all.models
  } else {
    dom.models <- models[grep("dom", models, fixed = T)]
    models <- setdiff(models, dom.models)
    if (length(dom.models) > 0) {
      dom.models <- unlist(lapply(as.list(dom.models), 
                                  function(x) {
                                    all.models[grep(x, all.models, fixed = T)]
                                  }))
    }
    models <- union(models, dom.models)
    stopifnot(all(is.element(models, all.models)))
  }
  scores <- as.data.frame(data@scores[[trait]][, models])
  colnames(scores) <- models
  scores$Chrom <- data@map$Chrom
  tmp <- pivot_longer(data = scores, cols = 1:length(models), 
                      names_to = "model", values_to = "y", values_drop_na = TRUE)
  tmp <- as.data.frame(tmp)
  tmp$model <- factor(tmp$model, levels = models, ordered = T)
  tmp <- tmp[order(tmp$model, tmp$Chrom, tmp$y, decreasing = c(FALSE, 
                                                               FALSE, TRUE)), ]
  tmp2 <- tapply(tmp$y, list(tmp$Chrom, tmp$model), function(x) {
    n <- length(x)
    unif.p <- -log10(ppoints(n))
  })
  tmp$x <- unlist(tmp2)
  return(tmp) 
}
