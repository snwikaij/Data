assess_curve <- function(b1=-0.3, b0=3, nsim=3000, cut=150){

  df <- array(NA, dim=c(nsim,18))

  for(i in 1:nrow(df)){

  minx <- runif(1, 20, 80)
  maxx <- runif(1, 500, 800)
  x    <- runif(100, minx, maxx)
  y    <- exp(b1*log(x)+b0)
  sdg  <- 1

  y <- rgamma(length(x), y^2/sdg^2, y/sdg^2)

  mod1 <- glm(y ~ log(x), family = Gamma(link = "log"))

  dat    <- data.frame(x=x, y=y)
  dat    <-dat[order(dat$x, decreasing = F),]
  low_g  <- dat[dat$x < cut,]
  hig_g  <- dat[dat$x >= cut,]
  mod2   <- glm(y ~ log(x), data = low_g, family = Gamma(link = "log"))
  mod3   <- glm(y ~ log(x), data = hig_g, family = Gamma(link = "log"))
  mod4   <- lm(y ~ x)
  mod5   <- lm(y ~ x, data = low_g)
  mod6   <-lm(y ~ x, data = hig_g)

  extr <- function(mod, true){

    ci        <- confint(mod)[2,]
    covered   <- ci[1] <= true & ci[2] >= true
    precision <- coef(mod)[2]-true

    c(covered, precision,summary(mod)$coef[2,4])}

  df[i,] <- c(extr(mod1, b1), extr(mod2, b1), extr(mod3, b1), extr(mod4, b1), extr(mod5, b1), extr(mod6, b1))

}

  df        <- as.data.frame(df)
  colnames(df) <- c("coverage_glm_full", "precision_glm_full", "sig_glm_full",
                    "coverage_glm_lower", "precision_glm_lower", "sig_glm_lower",
                    "coverage_glm_higher", "precision_glm_higher", "sig_glm_higher",
                    "coverage_lm_full", "precision_lm_full", "sig_lm_full",
                    "coverage_lm_lower", "precision_lm_lower", "sig_lm_lower",
                    "coverage_lm_higher", "precision_lm_higher", "sig_lm_higher")
  precision <- tidyr::gather(df[c(2,5,8,11,14,17)])

  pl1 <- ggplot(precision, aes(value))+
    geom_histogram(bins = 10, fill="grey", col="black")+
    facet_wrap(.~key, scales="free")+ylab("count")+xlab("estimate")+
    geom_vline(xintercept = 0, col="red", lty=2)+
    theme_classic()

  precision <- aggregate(precision, value~key, function(x) c(mean(x), sd(x)))
  precision <- setNames(cbind.data.frame(precision$key, precision$value), c("key", "mu", "sd"))

  pardata <- ggplot_build(pl1)
  textbox <- aggregate(data=pardata$data[[1]][c("PANEL", "y")],
                       y~PANEL, max)[,-1]
  textbox <- textbox-(0.3*as.numeric(textbox))

  coverage <- tidyr::gather(df[c(1,4,7,10,13,16)])
  coverage <- prop.table(table(coverage), margin = 1)

  sigtab       <- tidyr::gather(df[c(3,6,9,12,15,18)])
  sigtab$value <- ifelse(sigtab$value<0.05,1,0)
  sigtab       <- prop.table(table(sigtab), margin = 1)

  textbox  <- cbind.data.frame(precision, coverage[,2], sigtab[,2], textbox)
  colnames(textbox)[4:6] <- c("coverage", "power", "value")
  textbox$text <- paste("mean deviation \nfrom true estimate=", round(textbox$mu,2),
                        "\nstandard deviation \nof precision=", round(textbox$sd, 2),
                        "\ncoverage probability=", round(textbox$coverage,2),
                        "\npower=", round(textbox$power,2))

  pl1+geom_text(data=textbox, aes(x=0.1, y=value, label=text))}

figure <- assess_curve()

print(figure)

  set.seed(123)
  minx <- runif(1, 20, 80)
  maxx <- runif(1, 500, 800)
  x    <- runif(100, minx, maxx)
  y    <- exp(b1*log(x)+b0)
  sdg  <- 1

  y <- rgamma(length(x), y^2/sdg^2, y/sdg^2)

  xl          <- seq(20, 800, 1)
  points_data <- data.frame(x, y)
  line_full   <- data.frame(xl, yl=exp(log(xl)*-0.3+3))
  line_lower  <- data.frame(xl, yl=exp(log(xl)*-0.29+3))
  line_lower  <- line_lower[line_lower$xl<150,]
  line_higher <- data.frame(xl, yl=exp(log(xl)*-0.31+3))
  line_higher <- line_higher[line_higher$xl>=150,]

  pl2 <- ggplot(points_data, aes(x, y))+
    geom_point()+
    geom_line(data=line_full, aes(xl, yl), lwd=1, inherit.aes = F)+
    geom_line(data=line_lower, aes(xl, yl), lwd=1, col="darkgreen", lty=2, inherit.aes = F)+
    geom_line(data=line_higher, aes(xl, yl), lwd=1, col="tomato3",lty=2, inherit.aes = F)+
    theme_classic()

  line_full_lm    <- data.frame(xl, yl=xl*-0.0043+5.1)
  line_lower_lm   <- data.frame(xl, yl=xl*-0.012+6.45)
  line_lower_lm   <- line_lower_lm[line_lower_lm$xl<150,]
  line_higher_lm  <- data.frame(xl, yl=xl*-0.0031+4.6)
  line_higher_lm  <- line_higher_lm[line_higher_lm$xl>=150,]

  pl3 <- ggplot(points_data, aes(x, y))+
    geom_point()+xlim(20, 580)+
    geom_line(data=line_full_lm, aes(xl, yl), lwd=1, inherit.aes = F)+
    geom_line(data=line_lower_lm, aes(xl, yl), lwd=1, col="darkgreen", lty=2, inherit.aes = F)+
    geom_line(data=line_higher_lm, aes(xl, yl), lwd=1, col="tomato3",lty=2, inherit.aes = F)+
    theme_classic()

    cowplot::plot_grid(pl2, pl3, ncol=1, labels = c("GLM", "LM"))
