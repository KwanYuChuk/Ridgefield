library(ncdf4)   
library(raster) 
library(sf)      
library(ggplot2)
library(lubridate)
library(tidyr)
library(dplyr)
library(purrr)

# Open the NetCDF file
rf <- "ridgefield_2015_2024_L6.nc"
rf_data <- nc_open(rf)
View(rf_data)

# Check variable names
names(rf_data$var)

SKYE_IN_mV_Avg_MODIS1 <- ncvar_get(rf_data, "SKYE_IN_umol_Avg_MODIS1") 
SKYE_IN_mV_Avg_MODIS2 <- ncvar_get(rf_data, "SKYE_IN_umol_Avg_MODIS2") 
SKYE_OUT_mV_Avg_MODIS1 <- ncvar_get(rf_data, "SKYE_OUT_umol_Avg_MODIS1") 
SKYE_OUT_mV_Avg_MODIS2 <- ncvar_get(rf_data, "SKYE_OUT_umol_Avg_MODIS2") 
SKYE_IN_mV_Avg_NDVI1 <- ncvar_get(rf_data, "SKYE_IN_umol_Avg_NDVI1") 
SKYE_IN_mV_Avg_NDVI2 <- ncvar_get(rf_data, "SKYE_IN_umol_Avg_NDVI2") 
SKYE_OUT_mV_Avg_NDVI1 <- ncvar_get(rf_data, "SKYE_OUT_umol_Avg_NDVI1") 
SKYE_OUT_mV_Avg_NDVI2 <- ncvar_get(rf_data, "SKYE_OUT_umol_Avg_NDVI2") 

time <- ncvar_get(rf_data, "time")
time_units <- ncatt_get(rf_data, "time", "units")$value
origin <- as.Date(strsplit(time_units, "since ")[[1]][2])
dates <- origin + time

skye <- data.frame(
  date = dates,
  MODIS1_IN = SKYE_IN_mV_Avg_MODIS1,
  MODIS2_IN = SKYE_IN_mV_Avg_MODIS2,
  MODIS1_OUT = SKYE_OUT_mV_Avg_MODIS1,
  MODIS2_OUT = SKYE_OUT_mV_Avg_MODIS2,
  NDVI1_IN = SKYE_IN_mV_Avg_NDVI1,
  NDVI2_IN = SKYE_IN_mV_Avg_NDVI2,
  NDVI1_OUT = SKYE_OUT_mV_Avg_NDVI1,
  NDVI2_OUT = SKYE_OUT_mV_Avg_NDVI2
)

df_skye <- pivot_longer(skye, cols = -date, names_to = "sensor", values_to = "value")
df_skye

############################################################
#                                                          #
#     \       \      ___I_                                 #
#    / \     / \    /\- --\                                #
#    / \     / \   /  \_-__\                               #
#     |       |    |[]| [] |                               #  
#                                                          #
############################################################   

# 2021 ndvi wide band

df_skye <- df_skye %>%
  mutate(date = as.POSIXct(date, origin = "1970-01-01", tz = "UTC"))

df_skye_noon <- df_skye %>%
  filter(date >= as.POSIXct("2021-01-01", tz = "UTC") & 
           date < as.POSIXct("2022-01-01", tz = "UTC")) %>%
  filter(hour(date) == 12 & minute(date) == 0)

head(df_skye_noon)

# Step 1: Filter only the relevant sensors
modis_noon <- df_skye_noon %>%
  filter(value != -9999) %>%
  filter(sensor %in% c("MODIS1_IN", "MODIS2_IN", "MODIS1_OUT", "MODIS2_OUT"))
modis_noon
# Step 2: Pivot to wide format
modis_noon_wide <- modis_noon %>%
  pivot_wider(names_from = sensor, values_from = value)

ggplot(modis_noon, aes(x = date, y = value, color = sensor)) +
  geom_line(size = 1) +
  labs(
    title = "Wide Sensor Time Series (2021)",
    x = "",
    y = "raw sensor reading"
  ) +
  theme_minimal()

# Step 3: Calculate NDVI
ndvi_noon <- modis_noon_wide %>%
  mutate(
    ndvi = ((0.04694/0.06580)*MODIS2_OUT*MODIS1_IN - MODIS1_OUT*MODIS2_IN) / ((0.04694/0.06580)*MODIS2_OUT*MODIS1_IN + MODIS1_OUT*MODIS2_IN)
  ) %>%
  select(date, ndvi)

# View result
head(ndvi_noon)

ggplot(ndvi_noon, aes(x = date, y = ndvi)) +
  geom_line(color = "forestgreen", size = 1) +
  labs(
    title = "NDVI Time Series at 12:00 Noon (2021)",
    x = "Date",
    y = "NDVI"
  ) +
  theme_minimal()

ndvi_noon$date <- as.Date(ndvi_noon$date)
ndvi_noon$num_date <- as.numeric(ndvi_noon$date)

spline_fit <- smooth.spline(x = ndvi_noon$num_date, y = ndvi_noon$ndvi, spar = 0.6)  # Adjust spar (0–1) for smoothness

# Create smoothed NDVI dataframe
ndvi_spline <- data.frame(
  date = as.Date(spline_fit$x, origin = "1970-01-01"),
  ndvi = spline_fit$y
)
ndvi_spline
write.csv(ndvi_spline, "ndvi_wide.csv", row.names = FALSE)

# Plot
ggplot() +
  geom_point(data = ndvi_noon, aes(x = date, y = ndvi), color = "gray", alpha = 0.5) +
  geom_line(data = ndvi_spline, aes(x = date, y = ndvi), color = "forestgreen", size = 1) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(
    title = "NDVI time series (wide band)",
    x = " ",
    y = "NDVI"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1))

############################################################
#                                                          #
#     \       \      ___I_                                 #
#    / \     / \    /\- --\                                #
#    / \     / \   /  \_-__\                               #
#     |       |    |[]| [] |                               #  
#                                                          #
############################################################   

# 2021 ndvi narrow band

df_skye <- df_skye %>%
  mutate(date = as.POSIXct(date, origin = "1970-01-01", tz = "UTC"))

df_skye_noon <- df_skye %>%
  filter(date >= as.POSIXct("2021-01-01", tz = "UTC") &
           date < as.POSIXct("2022-01-01", tz = "UTC")) %>%
  filter(hour(date) == 12 & minute(date) == 0)

head(df_skye_noon)

# Step 1: Filter only the relevant sensors
narrow_noon <- df_skye_noon %>%
  filter(value != -9999) %>%
  filter(sensor %in% c("NDVI1_IN", "NDVI2_IN", "NDVI1_OUT", "NDVI2_OUT"))
narrow_noon
# Step 2: Pivot to wide format
narrow_noon2 <- narrow_noon %>%
  pivot_wider(names_from = sensor, values_from = value)

ggplot(narrow_noon, aes(x = date, y = value, color = sensor)) +
  geom_line(size = 1) +
  labs(
    title = "Narrow Sensor Time Series (2021)",
    x = "",
    y = "raw sensor reading"
  ) +
  theme_minimal()

# Step 3: Calculate NDVI
ndvi_noon <- narrow_noon2 %>%
  mutate(
    ndvi = ((0.01359/0.01380)*NDVI2_OUT*NDVI1_IN - NDVI1_OUT*NDVI2_IN) / ((0.01359/0.013800)*NDVI2_OUT*NDVI1_IN + NDVI1_OUT*NDVI2_IN)
  ) %>%
  select(date, ndvi)

# View result
head(ndvi_noon)

ggplot(ndvi_noon, aes(x = date, y = ndvi)) +
  geom_line(color = "forestgreen", size = 1) +
  labs(
    title = "NDVI Time Series at 12:00 Noon (2021)",
    x = "Date",
    y = "NDVI"
  ) +
  theme_minimal()

ndvi_noon$date <- as.Date(ndvi_noon$date)
ndvi_noon$num_date <- as.numeric(ndvi_noon$date)

df <- data.frame(x = ndvi_noon$num_date, y = ndvi_noon$ndvi)

df_clean <- df[is.finite(df$x) & is.finite(df$y), ]

spline_fit_narrow <- smooth.spline(x = df_clean$x, y = df_clean$y, spar = 0.6)

# Create smoothed NDVI dataframe
ndvi_spline_narrow <- data.frame(
  date = as.Date(spline_fit_narrow$x, origin = "1970-01-01"),
  ndvi = spline_fit_narrow$y
)
ndvi_spline_narrow
write.csv(ndvi_spline, "ndvi_narrow.csv", row.names = FALSE)

# Plot
ggplot() +
  geom_point(data = ndvi_noon, aes(x = date, y = ndvi), color = "gray", alpha = 0.5) +
  geom_line(data = ndvi_spline_narrow, aes(x = date, y = ndvi), color = "forestgreen", size = 1) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(
    title = "NDVI time series (narrow band)",
    x = " ",
    y = "NDVI"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1))
