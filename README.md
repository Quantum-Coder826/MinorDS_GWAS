# MinorDS_GWAS

## potatos.R
Dit script is een one stop voor onze data analyse. Produceer alleen de basis e.g.

- Losse QQ.plots voor elke trait
- Losse manhattan.plots voor elkte trait
- QTL table

QTL table staat in outputs naast de data.loco.scan object. 
Deze kan je met de onderstaande code inlezen zodat je niet de volledige analyse hoeft te doen.
```r
  data.loco.scan <- readRDS("./outputs/data_loco_scan.rds")
```

Je kan het hele scrip runnen door hem te openen en op "Source" te klicken 
of het onderstaande in de rstudio console te runnen.
```r
  source("~/Documents/opleiding/J3/MinorDS/MinorDS_GWAS/potatoes.R", echo = TRUE)
```

## qq.R
Bevat de functie `qq(data, trait, model = NULL)` is een gehackte versie van `qq.PLot()`.
Heeft de tidyverse nodig je kan de functie laden door `source("./qq.R")` dit load ook de libs.
deze returnd de data frame die gebuikt wordt om de plot de maken. Gestructureerd als volgt:

| Chrom | model | x | y |
|-------|-------|---|---|
| chr0  | additive| 0 | 0

Je moet zelf met ggplot het figuur maken de default plot is als volgt:

```r
qq(data, trait, model = NULL) %>%
  ggplot(aes(x, y, colour = model)) + 
    facet_wrap(~Chrom) + 
    geom_point() + 
    theme_bw() + 
    xlab(expression(paste("Expected -log"[10], "(p)", sep = ""))) + 
    ylab(expression(paste("Observed -log"[10], "(p)", sep = ""))) + 
    scale_colour_brewer(palette = "Set1") + 
    geom_abline(slope = 1, intercept = 0, linetype = 2) + 
    theme(text = element_text(size = 15))
```

Omdat we toegang hebben tot het onderstaande data structuur kan je meerdere chromosomen selecteren en complexere figuren maken.
Zie het onderstaande voorbeeld met drie chromosomen en een `facet_wrap()`.
```r
qq(data.loco.scan, trait="TRAIT") %>%   # Produceer de qq data, pipe in dplyr
  filter(Chrom == c("chr03", "chr04", "chr05")) %>%   # Filter voor chromesomen 3, 4 & 5
    ggplot(aes(x, y, colour = model)) +   # boots qq.PLot() na
    facet_wrap(~Chrom) +    # Maak 3 panelen in de plot voor elke chromesoom, qq.Plot() kan dit *niet*.
    geom_point() + 
    theme_bw() + 
    xlab(expression(paste("Expected -log"[10], "(p)", sep = ""))) +   # Dit doen we zo om de text fancy te printen.
    ylab(expression(paste("Observed -log"[10], "(p)", sep = ""))) + 
    scale_colour_brewer(palette = "Set1") + 
    geom_abline(slope = 1, intercept = 0, linetype = 2) + 
    theme(text = element_text(size = 15))
```

Als je verschillende losse plots wilt combineren in een dan raad ik gridExtra aan.
Zie de officiële [wiki](https://github.com/baptiste/gridextra/wiki) voor alle voorbeelden.