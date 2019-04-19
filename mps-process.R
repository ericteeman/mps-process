# Options -----------------------------------------------------------------

display.psf = "no"
display.har = "no"
export.data = "yes"
export.plots = "yes"
export.grid = "no"

# Import packages ---------------------------------------------------------

my.packages <- c("rChoiceDialogs", "purrr", "dplyr", "ggplot2", "ggpubr", "scales", "sfsmisc", "signal")
lapply(my.packages, require, character.only = TRUE)


# Functions ---------------------------------------------------------------

scientific_10_full <- function(x) {
  ifelse (x %% 1 == 0,
         parse(text = gsub("e+00", "", scientific_format()(x))),
         ifelse (x > 1 & x < 0.11, as.numeric(as.character(x)), 
                parse(text = gsub("e", " %*% 10^", scientific_format()(x)))))
}

scientific_10 <- function(x) {
  ifelse (x %% 1 == 0, parse(text = gsub("e+00", "", scientific_format()(x))),
         ifelse (x > 1 & x < 0.11, as.numeric(as.character(x)), parse(text = gsub(".*e", "10^", scientific_format()(x)))))
}

read.lvm <- function(flnm) {
  read.csv(flnm, skip = 35, sep = "\t", header = F,
    col.names = c("background.time","background.data", "sample.time", "sample.data", "psf.time","psf.data", 
                  "field.time","field.data","na","peakToPeak", "sample.time.full","sample.data.full", "comment"
    )
  ) %>%
    mutate(filename = flnm) %>%
    select(background.time, field.data, background.data, sample.data, psf.data) %>%
    group_by(background.time) %>%
    summarise_all(funs(mean)) %>%
    rename(time = background.time) %>%
    na.omit()
}

read.lvm.full <- function(flnm) {
  read.csv(flnm, skip = 35, sep = "\t", header = F,
    col.names = c("background.time","background.data", "sample.time", "sample.data", "psf.time","psf.data", 
                  "field.time","field.data","na","peakToPeak", "sample.time.full","sample.data.full", "comment"
    )
  ) %>%
    mutate(filename = flnm) %>%
    group_by(filename) %>%
    mutate(position = 1:n()) %>%
    ungroup() %>%
    select(sample.time.full, sample.data.full, position) 
    # rename(time = sample.time.full) %>%
}

read.lvm.field <- function(flnm) { read.csv( flnm, skip = 35, nrow = 1, sep = "\t", header = F) }

read.lvm.param <- function(flnm) { read.csv(flnm, skip = 24, nrow = 1, sep = "\t", header = F) }

read.conc <- function(flnm) { read.csv(flnm, header = FALSE, skip = 0) }

# Import data -------------------------------------------------------------

setwd("/Users/ericteeman/Google Drive/Research/Data/MPS/")
# setwd("/Users/ericteeman/Google Drive/Research/Data/MPS/Cell Internalization/20170801/1X")

setwd(rchoose.dir(caption = "Select Directory")) # Asks user to choose directory containing data files

nfiles <- length(list.files(pattern = "*\\d.lvm", full.names = T, recursive = F))

dat <- list.files(pattern = "*\\d.lvm", full.names = T, recursive = F) %>% map_df(~ read.lvm(.))
full <- list.files(pattern = "*\\d.lvm", full.names = T, recursive = F) %>% map_df(~ read.lvm.full(.))
field <- list.files(pattern = "*\\d.lvm", full.names = T, recursive = F) %>% map_df(~ read.lvm.field(.))
param <- list.files(pattern = "*\\d.lvm", full.names = T, recursive = F) %>% map_df(~ read.lvm.param(.))


# Spectrometer information ------------------------------------------------
# Sensitivity is receive coil sensitivity in [1/m] from [A/m/A].
# Spectrometer v2, built by Jack-Howell Clark in Fall 2013, S = 1989.44 1/m (0.0025 T/A).
# Frequency, sample rate, and field amplitude are all obtained from .lvm data files.

sensitivity = 1989.44 # 1/m
mu0 = 4 * pi * 1e-7
frequency = mean(param$V2) # Hz
omega = 2 * pi * frequency # Hz
sample.rate = mean(param$V4)
periods = mean(param$V14)
field.amplitude = mean(field$V10) / 2 # mT
vol <- 150 #uL

if (file.exists("conc.txt")) {
  conc <- list.files(pattern = "conc.txt", full.names = T, recursive = F) %>% 
    map_df(~ read.conc(.))
  conc <- conc[, 1] #mgFe/mL
  conc.new = conc / 1000 #gFe/mL
  vol.new = vol / 1000 #mL
  mass = conc.new * vol.new #gFe
  psf.label = expression(paste(chi, " [m" ^ "3", "gFe" ^ "-1", "]"))
  har.label = expression(paste("Amp. [Am" ^ "2", "gFe" ^ "-1", "]"))
} else {
  conc <- 0 #mgFe/mL
  mass <- 1 #filler value to prevent conversion without known conc
  psf.label = expression(paste(chi, " [x10"^"-9","m" ^ "3", "]"))
  har.label = expression(paste("Amp. [Am" ^ "2", "]"))
}


# Interpolate data points in primary data set
int = 7000

dat = dat %>%
  dplyr::filter(row_number() <= (n() / nfiles)) %>%
  add_row(., time = approx(.$time,n = int - nrow(.))$y) %>%
  arrange(time) %>%
  mutate(., field.data = approx(.$field.data, n = nrow(.))$y) %>%
  mutate(., background.data = approx(.$background.data, n = nrow(.))$y) %>%
  mutate(., sample.data = approx(.$sample.data, n = nrow(.))$y) %>%
  mutate(., psf.data = approx(.$psf.data, n = nrow(.))$y)

# Find fit coefficient of field to determine appropriate shift (phi)
model = nls(
  dat$field.data ~ -field.amplitude * cos(omega * (dat$time + phi)),
  data = dat,
  control = list(
    maxiter = 100000,
    minFactor = 1e-3,
    printEval = TRUE
  ),
  start = list(phi = 1),
  algorithm = "port"
)

dat = dat %>%
  mutate(field.fitted = predict(model)) %>%
  mutate(psf.m3 = psf.data / -(sensitivity * field.amplitude * omega * sin(omega * (
    time + coef(model)
  )))) %>%
  mutate(psf.m3 = psf.m3 / mass) %>%
  dplyr::filter(row_number() <= (n() / 2)) %>%
  dplyr::filter(field.data >= 0.95 * min(field.data) &
                  field.data <= 0.95 * max(field.data)) %>%
  mutate(field.data = round(field.data * 1000, 4)) %>%
  mutate(field.fitted = round(field.fitted * 1000, 4)) %>%
  mutate(psf.norm = psf.m3/max(psf.m3)) %>%
  mutate(direction = case_when(row_number() <= (n() / 2) ~ "forward", 
                               row_number() > (n() / 2) ~ "reverse")) %>%
  group_by(direction) %>%
  mutate(fwhm = case_when( length(field.fitted[field.fitted > field.fitted[psf.m3 == max(psf.m3)]][which.min(abs(psf.m3[field.fitted > field.fitted[psf.m3 == max(psf.m3)]] - max(psf.m3) / 2))] - field.fitted[field.fitted < field.fitted[psf.m3 == max(psf.m3)]][which.min(abs(psf.m3[field.fitted < field.fitted[psf.m3 == max(psf.m3)]] - max(psf.m3) / 2))]) == 0 ~ 1,
                           length(field.fitted[field.fitted > field.fitted[psf.m3 == max(psf.m3)]][which.min(abs(psf.m3[field.fitted > field.fitted[psf.m3 == max(psf.m3)]] - max(psf.m3) / 2))] - field.fitted[field.fitted < field.fitted[psf.m3 == max(psf.m3)]][which.min(abs(psf.m3[field.fitted < field.fitted[psf.m3 == max(psf.m3)]] - max(psf.m3) / 2))]) >= 1 ~ field.fitted[field.fitted > field.fitted[psf.m3 == max(psf.m3)]][which.min(abs(psf.m3[field.fitted > field.fitted[psf.m3 == max(psf.m3)]] - max(psf.m3) / 2))] - field.fitted[field.fitted < field.fitted[psf.m3 == max(psf.m3)]][which.min(abs(psf.m3[field.fitted < field.fitted[psf.m3 == max(psf.m3)]] - max(psf.m3) / 2))])) %>%
  # mutate(fwhm = 1) %>%
  mutate(mh.norm = cumsum(psf.m3)) %>%
  mutate(mh.norm = case_when(direction == "reverse" ~ max(mh.norm) - mh.norm,
                            direction == "forward" ~ mh.norm)) %>%
  mutate(mh.norm = 2 * ((mh.norm-min(mh.norm))/(max(mh.norm)-min(mh.norm))) - 1) %>%
  ungroup() %>%
  mutate(har = seq.int(n()) - 1) %>%        #shift to correct index
  mutate(amp = abs(fft(mh.norm)) ) %>%      #fft
  # mutate(amp = amp / n()) %>%
  # mutate(amp = amp * har) %>%               #adjusts harmonic intensities to be consistent with exp data
  mutate(fifthThird = amp[har == 5] / amp[har == 3])
  

full = full %>%
  group_by(position) %>%
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))
  # summarise_all(mean)

delf = sample.rate / nrow(full)

full = full %>%
  mutate(fft_mean = fft(sample.data.full_mean)) %>%
  mutate(fft_sd = fft(sample.data.full_sd)) %>%
  mutate(amp_mean = sqrt(((Re(fft_mean)^2) + (Im(fft_mean)^2)) / n())) %>%
  mutate(amp_sd = sqrt(((Re(fft_sd)^2) + (Im(fft_sd)^2)) / n())) %>%
  slice(1:(n()/2)) %>%
  mutate(frequencies = row_number() * delf) %>%
  mutate(amp_mean = (amp_mean * 2) / (sensitivity * mu0 * row_number() * delf * 2 * pi)) %>%
  mutate(amp_sd = (amp_sd * 2) / (sensitivity * mu0 * row_number() * delf * 2 * pi)) %>%
  slice(periods*row_number()+1) %>%
  mutate(har = seq.int(n())) %>%
  mutate(oddeven = case_when(har %% 2 != 0 ~ "odd",
                             har %% 2 == 0 ~ "even"))

dat <- dat %>%
  mutate(fifthThird = full$amp_mean[full$har == 5] / full$amp_mean[full$har == 3])

# mutate(mh.norm = cumsum(sample.data.full)) %>%
#   mutate(amp = abs(fft(mh.norm))) %>%
#   mutate(har = seq.int(n())) %>%
#   dplyr::filter(har %in% har.subset) %>%
#   slice(row_number()*65)


# Plots -------------------------------------------------------------------

theme_new <- function (base_size=24, base_line_size = 1) {
  theme_bw(base_size=base_size,
           base_family="") %+replace%
    theme(
      axis.text.x = element_text(size = base_size, margin = margin(t = 0.5*base_size)),
      axis.text.y = element_text(size = base_size, margin = margin(r = 0.5*base_size)),
      axis.title = element_text(size = base_size+2),
      axis.ticks.length = unit(-8, "pt"),
      panel.border = element_rect(size = base_line_size, fill=NA),
      panel.grid = element_blank(),
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.position = c(0.01, 0.99),
      legend.justification = c("left", "top"),
      legend.direction = "vertical",
      legend.title = element_blank(),
      legend.title.align = 0.5,
      legend.text = element_text(size = base_size-4),
      legend.text.align = 0.5
    )
}

data.set <- dat 

xlab = expression(paste("t [",mu,"s]"))
ylab = expression(paste(mu[0],"H [mT]"))
xmin = 1e6*(min(data.set$time) - 0.1*max(data.set$time))
xmax = 1e6*(max(data.set$time) + 0.1*max(data.set$time))
ymin = min(data.set$field.data) - 0.1*max(data.set$field.data)
ymax = max(data.set$field.data) + 0.1*max(data.set$field.data)

p1 = ggplot(data.set) +
  geom_line(aes(x = 1e6*time, y = field.data), size = 1, color = "black") +
  # geom_line(aes(x = 1e6*time, y = field.fitted), size = 1, color = "red") +
  scale_x_continuous(breaks = pretty_breaks(n=3), limits = c(xmin,xmax)) +
  scale_y_continuous(breaks = pretty_breaks(n=3)) +
  theme_new() +
  guides(col = guide_legend(ncol = 1)) +
  labs(x = xlab, y = ylab)

ylab = "Amp. [mV]"
ymin = 1e3*(min(data.set$background.data,data.set$sample.data, data.set$psf.data) - 0.1*max(data.set$background.data,data.set$sample.data, data.set$psf.data))
ymax = 1e3*(max(data.set$background.data,data.set$sample.data, data.set$psf.data) + 0.1*max(data.set$background.data,data.set$sample.data, data.set$psf.data))

p2 = ggplot(data.set) +
  theme_new() +
  # geom_point(aes(x = 10e5*time, y = background.data), size = 3, stroke = 1, color = "blue") +
  geom_line(aes(x = 1e6*time, y = 1e3*background.data, color = "black"), size = 1) +
  geom_line(aes(x = 1e6*time, y = 1e3*sample.data, color = "blue"), size = 1) +
  geom_line(aes(x = 1e6*time, y = 1e3*psf.data, color = "red"), size = 1) +
  scale_x_continuous(breaks = pretty_breaks(n=3), limits = c(xmin,xmax)) +
  scale_y_continuous(limits = c(ymin,ymax)) +
  scale_color_manual(values = c("black", "blue", "red"), labels = c("bkg", "spl", "dif")) +
  guides(col = guide_legend(ncol = 1)) +
  labs(x = xlab, y = ylab)

data.set <- dat #%>% dplyr::filter(field.fitted < 0.95*max(field.fitted) & field.fitted > 0.95*min(field.fitted))
xlab = expression(paste(mu[0],"H [mT]"))
ylab = psf.label
xmin = min(data.set$field.fitted) - 0.1*max(data.set$field.fitted)
xmax = max(data.set$field.fitted) + 0.1*max(data.set$field.fitted)
ymin = 1e9*(min(data.set$psf.m3) - 0.1*max(data.set$psf.m3))
ymax = 1e9*(max(data.set$psf.m3) + 0.1*max(data.set$psf.m3))

p3 = ggplot(data.set) +
  geom_line(aes(x = field.fitted, y = 1e9*psf.m3, group = direction), size = 1) +
  # geom_point(aes(x = field.fitted, y = psf.m3, group = direction, color = direction), size = 3) +
  scale_x_continuous(breaks = pretty_breaks(n=3), limits = c(xmin,xmax)) +
  scale_y_continuous(breaks = pretty_breaks(n=3), limits = c(ymin,ymax)) +
  theme_new() +
  guides(col = guide_legend(ncol = 1)) +
  labs(x = xlab, y = ylab)


ylab = "Amp. [a.u.]"
ymin = min(data.set$psf.norm) - 0.1*max(data.set$psf.norm)
ymax = max(data.set$psf.norm) + 0.1*max(data.set$psf.norm)

p4 = ggplot(data.set) +
  geom_line(aes(x = field.fitted, y = psf.norm, group = direction), size = 1) +
  # geom_point(aes(x = field.fitted, y = psf.norm, group = direction, color = direction), size = 2) +
  scale_x_continuous(breaks = pretty_breaks(n=3), limits = c(xmin,xmax)) +
  scale_y_continuous(breaks = pretty_breaks(n=3), limits = c(ymin,ymax)) +
  theme_new() +
  guides(col = guide_legend(ncol = 1)) +
  labs(x = xlab, y = ylab)

ymin = min(data.set$mh.norm) - 0.1*max(data.set$mh.norm)
ymax = max(data.set$mh.norm) + 0.1*max(data.set$mh.norm)

p5 = ggplot(data.set) +
  geom_line(aes(x = field.fitted, y = mh.norm, group = direction), size = 1) +
  scale_x_continuous(breaks = pretty_breaks(n=3), limits = c(xmin,xmax)) +
  scale_y_continuous(breaks = pretty_breaks(n=3), limits = c(ymin,ymax)) +
  theme_new() +
  guides(col = guide_legend(ncol = 1)) +
  labs(x = xlab, y = ylab)

har.subset = seq(1, 37, 2)

data.set <- full %>% dplyr::filter(har %in% har.subset)

xlab = "Harmonic"
ylab = har.label
xmin = min(data.set$har) - 0.1*max(data.set$har)
xmax = max(data.set$har) + 0.1*max(data.set$har)
ymin = min(data.set$amp_mean) - 0.1*max(data.set$amp_mean)
ymax = max(data.set$amp_mean) + 0.1*max(data.set$amp_mean)

p6 = ggplot(data.set) +
  geom_point(aes(x = har, y = amp_mean), shape = 1, size = 3) +
  geom_errorbar(aes(x = har, ymin=amp_mean-amp_sd, ymax=amp_mean+amp_sd), width=.5, position=position_dodge(.9)) +
  scale_x_continuous(breaks = pretty_breaks(n=3), limits = c(xmin,xmax)) +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), limits = c(ymin,ymax)) +
  annotation_logticks(sides = "l") +
  theme_new() +
  guides(col = guide_legend(ncol = 3)) +
  labs(x = xlab, y = ylab)

# Set export directory whether or not saving images is selected
main.directory = getwd()
export.directory = "export"
dir.create(file.path(main.directory, export.directory), showWarnings = FALSE)
setwd(file.path(main.directory, export.directory))

sc = 1    # Set scaling for saved images

if(export.data == "yes"){
  processed.data <- dat %>% select(field.fitted, psf.m3, psf.norm, mh.norm, direction, fwhm, fifthThird)
  write.csv(processed.data, "processed.csv", row.names = FALSE)
  processed.har <- full %>% select(har,amp_mean,amp_sd) %>% dplyr::filter(har %in% har.subset)
  write.csv(processed.har,"odd harmonics.csv", row.names = FALSE)
}

if(export.plots == "yes"){
  # ggsave("field-data.png", p1, width = 4.5 * sc, height = 4 * sc, dpi = "retina")
  # ggsave("raw-data.png", p2, width = 4.5 * sc, height = 4 * sc, dpi = "retina")
  ggsave("psf-m3.png", p3, width = 4.5 * sc, height = 4 * sc, dpi = "retina")
  ggsave("psf-norm.png", p4, width = 4.5 * sc, height = 4. * sc, dpi = "retina")
  ggsave("mh-norm.png", p5, width = 4.5 * sc, height = 4 * sc, dpi = "retina")
  ggsave("harmonics.png", p6, width = 4.5 * sc, height = 4 * sc, dpi = "retina")
}

if(export.grid == "yes"){
  combined <- ggarrange(p2, p3, p5, p6, hjust = -0.25,
                        labels = c("(a)", "(b)", "(c)", "(d)"), font.label = list(size = 24),
                        ncol = 2, nrow = 2)
  ggsave("combined.png", combined, width = 9 * sc, height = 8 * sc, dpi = "retina")
}

if (display.psf == "yes"){
  p3
}

if (display.har == "yes"){
  p6
}