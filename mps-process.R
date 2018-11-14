# Options -----------------------------------------------------------------

display.psf = "no"
export.data = "yes"
export.plots = "no"
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

<<<<<<< HEAD
<<<<<<< HEAD
setwd("/Users/ericteeman/Google Drive/Research/Data/MPS")
=======
setwd("/Users/ericteeman/Google Drive/Research/Data/MPS/Improving in vitro MPS/")
>>>>>>> ee45f614482f02ed28b57a42212f703e18a6fbcf
=======
setwd("/Users/ericteeman/Google Drive/Research/Data/MPS/Improving in vitro MPS/")
>>>>>>> ee45f614482f02ed28b57a42212f703e18a6fbcf
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
vol <- 100 #uL

if (file.exists("conc.txt")) {
  conc <- list.files(pattern = "conc.txt", full.names = T, recursive = F) %>% 
    map_df(~ read.conc(.))
  conc <- conc[, 1] #mgFe/mL
  conc.new = conc / 1000 #gFe/mL
  vol.new = vol / 1000 #mL
  mass = conc.new * vol.new #gFe
  psf.label = expression(paste(chi, " [m" ^ "3", "gFe" ^ "-1", "]"))
  har.label = expression(paste("Amplitude [Am" ^ "2", "gFe" ^ "-1", "]"))
} else {
  conc <- 0 #mgFe/mL
  mass <- 1 #filler value to prevent conversion without known conc
  psf.label = expression(paste(chi, " [m" ^ "3", "]"))
  har.label = expression(paste("Amplitude [Am" ^ "2", "]"))
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
  dplyr::filter(field.data >= 0.99 * min(field.data) &
                  field.data <= 0.99 * max(field.data)) %>%
  mutate(field.data = round(field.data * 1000, 4)) %>%
  mutate(field.fitted = round(field.fitted * 1000, 4)) %>%
  mutate(psf.norm = psf.m3/max(psf.m3)) %>%
  mutate(direction = case_when(row_number() <= (n() / 2) ~ "forward", 
                               row_number() > (n() / 2) ~ "reverse")) %>%
  group_by(direction) %>%
  mutate(fwhm = field.fitted[field.fitted > field.fitted[psf.m3 == max(psf.m3)]][which.min(abs(psf.m3[field.fitted > field.fitted[psf.m3 == max(psf.m3)]] - max(psf.m3) / 2))] - field.fitted[field.fitted < field.fitted[psf.m3 == max(psf.m3)]][which.min(abs(psf.m3[field.fitted < field.fitted[psf.m3 == max(psf.m3)]] - max(psf.m3) / 2))]) %>%
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
  summarise_all(mean)

delf = sample.rate / nrow(full)

full = full %>%
  mutate(fft = fft(sample.data.full)) %>%
  mutate(amp = sqrt(((Re(fft)^2) + (Im(fft)^2)) / n())) %>%
  slice(1:(n()/2)) %>%
  mutate(frequencies = row_number() * delf) %>%
  mutate(amp = (amp * 2) / (sensitivity * mu0 * row_number() * delf * 2 * pi)) %>%
  slice(periods*row_number()+1) %>%
  mutate(har = seq.int(n())) %>%
  mutate(oddeven = case_when(har %% 2 != 0 ~ "odd",
                             har %% 2 == 0 ~ "even"))

dat <- dat %>%
  mutate(fifthThird = full$amp[full$har == 5] / full$amp[full$har == 3])

# mutate(mh.norm = cumsum(sample.data.full)) %>%
#   mutate(amp = abs(fft(mh.norm))) %>%
#   mutate(har = seq.int(n())) %>%
#   dplyr::filter(har %in% har.subset) %>%
#   slice(row_number()*65)


# Plots -------------------------------------------------------------------

data.set <- dat 

xlab = "Time [sec]"
ylab = expression(paste("H [mT", mu[0] ^ -1, "]"))

p1 = ggplot(data.set) +
  geom_point(aes(x = time, y = field.data), size = 3, color = "black") +
  geom_line(aes(x = time, y = field.fitted), size = 1, color = "red") +
  scale_x_continuous(labels = scientific_10_full) +
  scale_y_continuous(breaks = waiver()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14, margin = margin(t = 12)),
    axis.text.y = element_text(size = 14, margin = margin(r = 12)),
    axis.title = element_text(size = 16),
    axis.ticks.length = unit(-8, "pt"),
    panel.grid = element_blank(),
    panel.border = element_rect(size = 1, color = "black"),
    legend.position = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.text = element_text(size = 10)
  ) +
  guides(col = guide_legend(ncol = 1)) +
  labs(x = xlab, y = ylab)

ylab = "Amplitude [V]"

p2 = ggplot(data.set) +
  geom_line(aes(x = time, y = background.data), size = 1, color = "blue") +
  geom_line(aes(x = time, y = sample.data), size = 1, color = "black") +
  geom_line(aes(x = time, y = psf.data, group = direction, color = direction), size = 1) +
  scale_x_continuous(labels = scientific_10_full) +
  scale_y_continuous(breaks = waiver()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14, margin = margin(t = 12)),
    axis.text.y = element_text(size = 14, margin = margin(r = 12)),
    axis.title = element_text(size = 16),
    axis.ticks.length = unit(-8, "pt"),
    panel.grid = element_blank(),
    panel.border = element_rect(size = 1, color = "black"),
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.text = element_text(size = 10)
  ) +
  guides(col = guide_legend(ncol = 1)) +
  labs(x = xlab, y = ylab)


xlab = expression(paste("H [mT", mu[0] ^ -1, "]"))
ylab = psf.label

p3 = ggplot(data.set) +
  geom_line(aes(x = field.fitted, y = psf.m3, group = direction, color = direction), size = 1) +
  # geom_point(aes(x = field.fitted, y = psf.m3, group = direction, color = direction), size = 3) +
  scale_x_continuous(breaks = waiver()) +
  scale_y_continuous(labels = scientific_10_full) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14, margin = margin(t = 12)),
    axis.text.y = element_text(size = 14, margin = margin(r = 12)),
    axis.title = element_text(size = 16),
    axis.ticks.length = unit(-8, "pt"),
    panel.grid = element_blank(),
    panel.border = element_rect(size = 1, color = "black"),
    legend.position = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  ) +
  guides(col = guide_legend(ncol = 1)) +
  labs(x = xlab, y = ylab)


ylab = "Amplitude [a.u.]"

p4 = ggplot(data.set) +
  geom_line(aes(x = field.fitted, y = psf.norm, group = direction, color = direction), size = 1) +
  # geom_point(aes(x = field.fitted, y = psf.norm, group = direction, color = direction), size = 2) +
  scale_x_continuous(breaks = waiver()) +
  scale_y_continuous(breaks = waiver()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14, margin = margin(t = 12)),
    axis.text.y = element_text(size = 14, margin = margin(r = 12)),
    axis.title = element_text(size = 16),
    axis.ticks.length = unit(-8, "pt"),
    panel.grid = element_blank(),
    panel.border = element_rect(size = 1, color = "black"),
    legend.position = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  ) +
  guides(col = guide_legend(ncol = 1)) +
  labs(x = xlab, y = ylab)

p5 = ggplot(data.set) +
  geom_line(aes(x = field.fitted, y = mh.norm, group = direction, color = direction), size = 1) +
  scale_x_continuous(breaks = waiver()) +
  scale_y_continuous(breaks = waiver()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 15, margin = margin(t = 12)),
    axis.text.y = element_text(size = 15, margin = margin(r = 12)),
    axis.title = element_text(size = 20),
    axis.ticks.length = unit(-8, "pt"),
    panel.grid = element_blank(),
    panel.border = element_rect(size = 0.75, color = "black"),
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  ) +
  guides(col = guide_legend(ncol = 1)) +
  labs(x = xlab, y = ylab)

har.subset = seq(1, 37, 2)

data.set <- full %>% dplyr::filter(har %in% har.subset)

xlab = "Harmonic"
ylab = har.label
# , breaks = c(seq(5,10,5) %o% 10^(-20:20))

p6 = ggplot(data.set) +
  geom_point(aes(x = har, y = amp), size = 3) +
  scale_x_continuous(breaks = waiver()) +
  scale_y_log10(labels = scientific_10,breaks = c(10 ^ (-50:50))) +
  annotation_logticks(sides = "l") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 15, margin = margin(t = 12)),
    axis.text.y = element_text(size = 15, margin = margin(r = 12)),
    axis.title = element_text(size = 20),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(-8, "pt"),
    panel.grid = element_blank(),
    panel.border = element_rect(size = 0.75, color = "black"),
    legend.position = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  ) +
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
  processed.har <- full %>% select(har,amp) %>% dplyr::filter(har %in% har.subset)
  write.csv(processed.har,"odd harmonics.csv", row.names = FALSE)
}

if(export.plots == "yes"){
  ggsave("raw-data.png", p2, width = 6 * sc, height = 4.5 * sc, dpi = 300)
  ggsave("psf-m3.png", p3, width = 6 * sc, height = 4.5 * sc, dpi = 300)
  ggsave("psf-norm.png", p4, width = 6 * sc, height = 4.5 * sc, dpi = 300)
  ggsave("mh-norm.png", p5, width = 6 * sc, height = 4.5 * sc, dpi = 300)
  ggsave("harmonics.png", p6, width = 6 * sc, height = 4.5 * sc, dpi = 300)
}

if(export.grid == "yes"){
  grid <- ggarrange(p1, p2, p3, p4, p5, p6, ncol = 2, nrow = 3)
  ggsave("all-plot.png", grid, width = 8.5 * sc, heigh = 11 * sc, dpi = 300)
}

if (display.psf == "yes"){
  p3
}
