library(ggplot2)

asses_bias <- function(b1=-0.3, b0=3, nsim=3000){

df <- array(NA, dim=c(nsim,6))

for(i in 1:nrow(df)){

minx <- runif(1, 30, 300)
maxx <- runif(1, 500, 700)
x    <- runif(100, minx, maxx)
y    <- exp(b1*log(x)+b0)
sdg  <- 1

y <- rgamma(length(x), y^2/sdg^2, y/sdg^2)
minmax <- function(x) (x-min(x))/(max(x)-min(x))
zscale <- function(x) (x-mean(x))/sd(x)

mod1 <- glm(y ~ log(x), family = Gamma(link = "log"))
mod2 <- lm(minmax(y)~minmax(x))
mod3 <- lm(zscale(y)~zscale(x))

extr <- function(mod, true){

  ci        <- confint(mod)[2,]
  covered   <- ci[1] <= true & ci[2] >= true
  precision <- coef(mod)[2]-true

  c(covered, precision)}

df[i,] <- c(extr(mod1, b1), extr(mod2, b1), extr(mod3, b1))}

df        <- as.data.frame(df)
colnames(df) <- c("coverage_glm", "precision_glm",
                  "coverage_minmax", "precision_minmax",
                  "coverage_zscale", "precision_zscale")
precision <- tidyr::gather(df[c(2,4,6)])

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
textbox <- textbox-(0.2*as.numeric(textbox))

coverage <- tidyr::gather(df[c(1,3,5)])
coverage <- prop.table(table(coverage), margin = 1)

textbox  <- cbind.data.frame(precision, coverage[,2], textbox)
colnames(textbox)[4:5] <- c("coverage", "value")
textbox$text <- paste("mean deviation \nfrom true estimate=", round(textbox$mu,2),
                  "\nStandard deviation \nof precision=", round(textbox$sd, 2),
                  "\ncoverage probability=", round(textbox$coverage,2))

pl1+geom_text(data=textbox, aes(x=-0.05, y=value, label=text))}
                       
figure <- asses_bias()

print(figure)
