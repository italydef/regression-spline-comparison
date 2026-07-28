## functions.R

library(ggplot2)   # test-function panels (geom_point / geom_line)
library(cowplot)   # ggdraw() / draw_plot() to compose the 6-panel figure

WORKING_DIR  <- getwd()
OUTPUT_DIR <- file.path(WORKING_DIR, "OUTPUT", "Simulation")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Test functions
# The six piecewise-smooth benchmark functions (f1-f6) compared across
# fitting methods in the paper: each stresses a different kind of local
# irregularity (jump discontinuity, spike/peak, smooth bump, multiple
# jumps, sign-based discontinuities).

# f1: single jump discontinuity at x = 0.5 (flat left of the jump, smooth
# oscillation to the right).
f1 <- function(x) {
  y <- ifelse(x < 0.5, -1.5, 0.25 * sin(x ^ 2 * pi ^ 1.5))
  return(y)
}

# f2: two narrow spikes/peaks (near x = 0.3 and x = 0.65) from the
# 1/(eps + (x-c)^2) terms, joined at x = 0.6.
f2 <- function(x) {
  y <-
    ifelse(x < 0.6, 1 / (0.01 + (x - 0.3) ** 2), 1 / (0.015 + (x - 0.65) **
                                                        2))
  return(y)
}

# f3: smooth, symmetric bump (no discontinuities) - a baseline "easy" case.
f3 <- function(x) {
  y <- 100 / (exp(abs(10 * x - 5))) + ((10 * x - 5) ^ 5) / 500
  return(y)
}

# f4: smooth oscillation with three vertical jumps at x = 0.4, 0.6, 0.8.
f4 <- function(x) {
  y <- sin(15 * x) + 0.3 * x ^ 2

  jump1 <- 0.4
  jump2 <- 0.6
  jump3 <- 0.8

  y[x < jump1] <- y[x < jump1] + 2
  y[x >= jump1 & x < jump2] <- y[x >= jump1 & x < jump2] - 2
  y[x >= jump2 & x < jump3] <- y[x >= jump2 & x < jump3] + 1

  return(y)
}

# f5: oscillation with jumps induced via sign() at specific grid indices
# (x[59], x[78], x[180]) - depends on the length/spacing of x (n = 200
# below), not on fixed x-values.
f5 <- function(x){
      result <- sin(5.5*pi*x) - 4*sign(x[59]-x) -
           2*sign(x[78]-x) +3*sign(x[180]-x)-1.75
      return(result)
}

# f6: oscillation with a jump at x = 0.7 (sign term) plus a non-smooth
# kink at x = 0.4 (fractional-power |x-0.4|^(3/10) term).
f6 <- function(x) {
  result <- 2*sin(4*pi*x)-6*abs(x-0.4)^(3/10)-sign(0.7-x)
  return(result)
}

# Evaluate `func` on grid `x`, then add Gaussian noise scaled so that the
# ratio of signal variance to noise variance equals SNR.
simulate_data <- function(func, x, SNR) {
  f <- func(x)
  sigma<- sqrt(var(f))/SNR
  noise<- rnorm(length(x), mean = 0, sd = sigma)
  y<-f+noise
  data <- data.frame(x=x, y=y, f=f)
  return(data)
}

n<-200
x<-seq(0,1, length=n)
SNR=3
set.seed(123)  # fixed seed -> reproducible noisy samples for the figure

data1<-simulate_data(f1,x,SNR)

p1 <- ggplot() + geom_point(data = data1, aes(x = x, y = y), shape = 16,
      size = 1, col = rgb(0, 0, 1, alpha = 0.8)) +
      geom_line(data = data1, aes(x = x, y = f), color = "black",
      linewidth = 1) + 
      labs(x = "x", y = "(a) f1") + theme_minimal()
  
data2<-simulate_data(f2,x,SNR)
 
p2 <- ggplot() + geom_point(data = data2, aes(x = x, y = y), shape = 16,
      size = 1, col = rgb(0, 0, 1, alpha = 0.8)) +
      geom_line(data = data2, aes(x = x, y = f), color = "black",
      linewidth = 1) + 
      labs(x = "x", y = "(b) f2") +  theme_minimal()


data3<-simulate_data(f3,x,SNR)
 
p3 <- ggplot() + geom_point(data = data3, aes(x = x, y = y), shape = 16,
      size = 1, col = rgb(0, 0, 1, alpha = 0.8)) +
      geom_line(data = data3, aes(x = x, y = f), color = "black",
      linewidth = 1) + 
      labs(x = "x", y = "(c) f3") + theme_minimal()


data4<-simulate_data(f4,x,SNR)
 
p4 <- ggplot() + geom_point(data = data4, aes(x = x, y = y), shape = 16,
      size = 1, col = rgb(0, 0, 1, alpha = 0.8)) +
      geom_line(data = data4, aes(x = x, y = f), color = "black",
      linewidth = 1) + 
      labs(x = "x", y = "(d) f4") + theme_minimal()


data5<-simulate_data(f5,x,SNR)
 
p5 <- ggplot() + geom_point(data = data5, aes(x = x, y = y), shape = 16,
      size = 1, col = rgb(0, 0, 1, alpha = 0.8)) +
      geom_line(data = data5, aes(x = x, y = f), color = "black",
      linewidth = 1) + 
      labs(x = "x", y = "(e) f5") + theme_minimal()


data6<-simulate_data(f6,x,SNR)
 
p6 <- ggplot() + geom_point(data = data6, aes(x = x, y = y), shape = 16,
      size = 1, col = rgb(0, 0, 1, alpha = 0.8)) +
      geom_line(data = data6, aes(x = x, y = f), color = "black",
      linewidth = 1) + 
      labs(x = "x", y = "(f) f6") + theme_minimal()

# -----------------------------------------------------------------------------
# Arrange the six panels (f1-f6) into a single figure and save it to PDF.
# Produces plots of the true and noisy functions used for the simulation.
#
# NOTE: the composed plot must be assigned and then print()-ed to the device.
# When the script is source()-d, auto-printing is off, so writing
# pdf(); ggdraw() + draw_plot(...); dev.off()
# opens the file, draws nothing, and closes it -> an empty PDF that cannot be
# opened. print(combined_plot) forces the figure onto the open pdf() device.
# -----------------------------------------------------------------------------
combined_plot <- ggdraw() +
   draw_plot(p1, x = 0,   y = .5, width = .3, height = .5) +
   draw_plot(p2, x = .3,  y = .5, width = .3, height = .5) +
   draw_plot(p3, x = 0.6, y = .5, width = .3, height = .5) +
   draw_plot(p4, x = 0,   y = 0,  width = .3, height = .5) +
   draw_plot(p5, x = .3,  y = 0,  width = .3, height = .5) +
   draw_plot(p6, x = 0.6, y = 0,  width = .3, height = .5)

out_pdf <- file.path(OUTPUT_DIR, "true_noisy_functions.pdf")
pdf(file = out_pdf, width = 7, height = 5, onefile = TRUE, title = "Simulated examples")
print(combined_plot)
dev.off()
