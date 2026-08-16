#### ===> Analysis set up <=== ####
  ### ==> Packages
  library(tidyverse)
  library(janitor)
  library(patchwork)
  library(viridis)
  library(ggrepel) # for dbRDA
  ## Analyses
  library(vegan)
  library(pairwiseAdonis) #(https://github.com/pmartinezarbizu/pairwiseAdonis
  library(permute)
  library(suncalc)
  library(mgcv)
  ## PSD
  library(zoo)
  library(psd)        # Welch/multitaper-style adaptive PSD
  library(biwavelet)  # continuous wavelet transform
  ## DHW
  library(fuzzyjoin)
  library(sf)
  ### ==> colour pallete
  pal_cols = c("#8da0cb","#66c2a5","#fc8d59")
                # outer ~ mid-shore ~ inner
  
# ========================================================================================
#
#                       ####  ~~~~  NOAA SST & DHW  ~~~ ####
#
# ========================================================================================
# 
# ~ Note: NOAA CRW data from 1986 to 2024
#  
#### ===> Reading Coral Reef Watch data and site meta data <=== ####
  ### ==> CRW data
  crw<-read.csv("Data/Palmyra_CRW_SST_1986-2024.10.31.csv",stringsAsFactors = TRUE)%>%clean_names()
  
  ### ==> Site meta data
  sites<-read.csv("Data/Pal_terrace_metadata.csv",header=T)%>%clean_names()
  
  ### ==> Read Palmyra shapefile
  pal_map<-st_read("Data/repalmyradhwcode/LI_smooth.shp")
  pal_map$L4_ATTRIB<-ordered(factor(pal_map$L4_ATTRIB,levels=c('forereef','reef flat','land on reef')))
  
  ### ==> Filter data to pixels closest to the sites
    # Note: in this particular instance the sites are all too close together to be matched with independent so isolating to the closest pixel (max_dist=1km)
  sstPix<-sites %>% 
    geo_join(crw,by=c("longitude","latitude"),max_dist=1,unit="km",mode='inner') %>%
    mutate(my = format(as.Date(date), "%Y-%m")) %>%
    mutate(date = as.POSIXct(date, format="%Y-%m-%d")) %>%
    mutate(doy = yday(date),
           month = month(date,label = TRUE),
           year = year(date))%>%
    select(date, doy, month, year, latitude.y, longitude.y, crw_dhw, crw_sst, crw_hotspot, crw_sstanomaly) %>%
    rename(latitude = latitude.y, longitude=longitude.y) %>% as.data.frame()
    # check data
    glimpse(sstPix)
 
#### ===> Long-term climatology with ENSO years <=== ####
  ### ==> Set up
    ## ENSO year
    hl_years <- c(1998, 2009, 2015)  
    
    ## spatial extent
    maxlat=5.95       #needs to be the larger (or less negative if in S. Hemi) number
    minlat=5.85
    minlon=-162.17    #needs to be the more negative number if in the W hemi
    maxlon=-162.0
    
    ## extract maximum DHW in time series for each pixel
    anomaly<-crw %>% group_by(latitude,longitude) %>% slice(which.max(crw_dhw))
    
  ### ==> Pixel coverage visualisation
  plot.crw_pixel<-ggplot()+
    geom_raster(data=anomaly,mapping=aes(longitude,latitude,fill=crw_dhw),interpolate = F)+
    geom_sf(data=pal_map, fill=NA, col="black")+#geom_sf(data=lag,aes(fill=RB_ATTRIB))+
    geom_point(data=sites,aes(x=longitude,y=latitude),pch=21,col='black', size=3, fill="white")+
    ylim(minlat,maxlat)+
    xlim(minlon,maxlon)+
    scale_fill_viridis(option='magma')+
    labs(x="Longitude", y="Latitude")+
    theme_classic()+
    theme(legend.position = "none")
  plot.crw_pixel
  
  # save plot
  #ggsave("Feedback_Outputs/Supp_CRW_pixel.png", plot.crw_pixel, dpi = 600, width = 5, height = 3.5, units = "in", scale = 1)  
#
#
#
# ========================================================================================
#
#                       ####  ~~~~  USWFS Western Terrace Logger  ~~~ ####
#
# ========================================================================================
#
# ~ Note: Defining daytime as time between avg sunrise (06:35) and avg sunset (18:35) time
# ~ Note: Defining night time as time between avg sunset (18:35) and avg sunrise (06:35) time 
#
#### ===> Reading USWFS terrace data <=== ####  
uswf<-read.csv("Data/WesternTerrace_FWS_hourlyTemp_2016_2025.csv",stringsAsFactors = TRUE)%>%
    clean_names()%>%  
    select(date, temp_c)%>%
    #Convert date from factor to date variable  
    mutate(datetime = as.POSIXct(date, format="%Y-%m-%d %H:%M:%S"),
           date = as.POSIXct(date, format="%Y-%m-%d"))
    ## check data
    glimpse(uswf) 
    
  ### ==> Calculate daily mean night time SST
  # Define night time start and end points  
  sunrise <- hms("06:30:00")
  sunset <- hms("18:35:00")
  
  # Filter to night time period
  uswf.night <- uswf %>%
    filter(hms(format(datetime, "%H:%M:%S")) <= sunrise | 
             hms(format(datetime, "%H:%M:%S")) >= sunset)
  
  # Summarise 
  uswf.night_sum<-uswf.night%>%
    group_by(date)%>%
    summarise(temp_c = mean(temp_c))%>%
    mutate(doy = yday(date),
           month = month(date,label = TRUE),
           year = year(date))
#
#
#
# ========================================================================================
#
#                  ####  ~~~~  NOAA Tidal data  ~~~ ####
#
# ========================================================================================
#### ===> Read tidal data for 2023 <=== ####
# Palmyra atoll lat and long
lat <- 5.8833        
lon <- -162.0833   
  ### ==> high-low from subordinate station
  tides_pal<-read.csv("Data/TPT2739_tides_23.csv",stringsAsFactors = TRUE)%>%
    clean_names()%>%
    mutate(datetime = as.POSIXct(date_time, format="%d/%m/%Y %H:%M", tz = "Pacific/Kiritimati"))
    # Cosine (harmonic) interpolation between successive high/low waters, onto a
    # 10-min grid. This is the standard approximation for reconstructing a tidal
    # curve from turning points only.
    tide_pal.interp <- do.call(rbind, lapply(1:(nrow(tides_pal) - 1), function(i) {
      t1 <- tides_pal$datetime[i];   t2 <- tides_pal$datetime[i + 1]
      h1 <- tides_pal$prediction[i]; h2 <- tides_pal$prediction[i + 1]
      tt <- seq(t1, t2, by = "10 min")
      frac <- as.numeric(difftime(tt, t1, units = "secs")) /
        as.numeric(difftime(t2, t1, units = "secs"))
      data.frame(datetime = tt,
                 height = (h1 + h2)/2 - (h2 - h1)/2 * cos(pi * frac))
    }))%>%
      distinct(datetime, .keep_all = TRUE) %>%
      mutate(dh = c(NA, diff(height)) / c(NA, as.numeric(diff(datetime), units = "hours")))
    
    ## extract high water points
    hw <- tide_pal.interp %>%
      mutate(lead_dh = dplyr::lead(dh)) %>%
      filter(dh > 0, lead_dh <= 0) %>%
      pull(datetime)
    
    ## function to extract total hours since high water
    tidal_hour <- function(t, hw) {
      i <- findInterval(t, hw)
      out <- as.numeric(difftime(t, hw[pmax(i, 1)], units = "hours"))
      out[i == 0] <- NA_real_
      out
    }
    

#### ===> from reference station for Palmyra atoll: Honolulu station <=== ####
url <- paste0(
  "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?",
  "begin_date=20230926&end_date=20231007",
  "&station=1612340",          # the harmonic reference station ID
  "&product=predictions",
  "&datum=MLLW",
  "&time_zone=lst_ldt",        # or gmt — match your PAR logger
  #"&interval=h",               # 'h' = hourly; omit entirely for 6-minute
  "&units=metric",
  "&format=csv")

tide_ref <- read.csv(url)%>%
  mutate(datetime = as.POSIXct(Date.Time, format="%Y-%m-%d %H:%M"))
# ========================================================================================
#
#                  ####  ~~~~  Site specific in-situ logger data  ~~~ ####
#
# ========================================================================================
#
# ~ Note: 2023 data: 26/09/23 to 07/10/23, recorded at 5 minute intervals
#       -> Variables: Flow (speed, north, east, heading), PAR, Temperature, DO, DSAT, AOU
# ~ Note: 2024 data: 05/10/24 to 15/10/24, recorded at 1 minute intervals
#       -> Variables: Flow (speed, north, east, heading)
#
# Define night time start and end points  
sunrise <- hms("06:30:00")
sunset <- hms("18:35:00")

#### ===> Read data files for 2023 <=== ####
env23<-read.csv("Data/Env_variables_23.csv",stringsAsFactors = TRUE)%>%
clean_names()%>%
mutate(
# convert date and time to a date-time variable
  datetime = as.POSIXct(paste(date, time), format="%d/%m/%Y %H:%M:%S"),
# convert date to date variable
  date = as.POSIXct(date, format="%d/%m/%Y"),
# calculate instantaneous par 
  # Data is in 5 min intervals -> need to multiply by 300 (5 min x 60s/min) to get instantaneous light measurements as per Morgan et al. (2020)
  par.inst_umol = (par*300),  
  # convert from u mol/m2 s to mol/m2 s -> divide all values by 10^6
  par.inst = (par.inst_umol/10^6),
# adjust site names to reflect distance from shore
  site = factor(site, levels=c("RT1","RT4","RT7"), labels = c("Offshore", "Midshore","Inshore"))
)%>%
# rename north and east and v and u
  rename(u = east, 
         v = north)

### ==> Read data files for 2024 <=== ####
env24<-read.csv("Data/Env_flow_24.csv",stringsAsFactors = TRUE)%>%
  clean_names()%>%
  mutate(
    # convert date and time to a date-time variable
    datetime = as.POSIXct(date_time, format="%d/%m/%Y %H:%M"),
    # convert date to date variable
    date = as.POSIXct(date_time, format="%d/%m/%Y"),
    # adjust site names to reflect distance from shore
    site = factor(site, levels=c("RT1","RT4","RT7"), labels = c("Offshore", "Midshore","Inshore"))
  )%>%
  # rename north and east and v and u
  rename(u = east, 
         v = north)

#### ===> 24 hour cycle calculation <=== ####
  ### ==> 2023 data
  diel23<-env23%>%
      select(datetime, site,par, temp, do, dsat, aou, speed, u, v)%>%
      mutate(tod = hour(datetime) + minute(datetime) / 60)%>%
      pivot_longer(c(par, temp, do, dsat, aou, speed, u, v),
                  names_to = "variable", 
                  values_to = "value")%>%
      group_by(site, variable, tod)%>%
      summarise(mean = mean(value, na.rm = TRUE),
                se   = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))))%>%
      mutate(variable=as.factor(variable))%>%as.data.frame()
  
  ### ==> 2024 data
  diel24<-env24%>%
    select(datetime, site,speed, u, v)%>%
    mutate(tod = hour(datetime) + minute(datetime) / 60)%>%
    pivot_longer(c(speed, u, v),
                 names_to = "variable", 
                 values_to = "value")%>%
    group_by(site, variable, tod)%>%
    summarise(mean = mean(value, na.rm = TRUE),
              se   = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))))%>%
    mutate(variable=as.factor(variable))%>%as.data.frame()
  
#### ===> Isolate temperature data to plot with USWF data <=== ####
  ### ==> Filter to night time SST
  temp_night<-env23%>%
  select(date, datetime, site, temp)%>%
  filter(hms(format(datetime, "%H:%M:%S")) <= sunrise | 
           hms(format(datetime, "%H:%M:%S")) >= sunset)
    ## Summarise
    temp_night_sum<-temp_night%>%
      group_by(date, site)%>%
      summarise(temp_c = mean(temp))%>%
      mutate(doy = yday(date),
             month = month(date,label = TRUE),
             year = year(date))%>%as.data.frame()

#### ===> Calculate the Daily light integral <=== ####    
DLI<-env23%>%
      filter(hms(format(datetime, "%H:%M:%S")) >= sunrise &
               hms(format(datetime, "%H:%M:%S")) <= sunset)%>%
      group_by(site, date)%>%
      summarise(DLI = sum(par.inst))%>%as.data.frame()
#
#
#
# ========================================================================================
#
#                       ####  ~~~~  PLOT: SST and DWH   ~~~ ####
#
# ========================================================================================
#### ===> NOAA SST plot <=== ####
plot.crw_sst<-ggplot()+
  # all years except ENSO years
  geom_line(data = sstPix%>%filter(year != hl_years), 
            aes(x=doy, y=crw_sst, group = year), color="gray80", alpha=0.9)+
  # plot ENSO years
  geom_line(data = sstPix%>%filter(year == hl_years),
            aes(x=doy, y=crw_sst, group = year, color=factor(year)),linewidth = 0.9, alpha = 0.8)+
  # plot  2023 
  geom_line(data = sstPix%>%filter(year == 2023),
            aes(x=doy, y=crw_sst, group = year, color=factor(year)), linewidth = 1 )+
  # highlight sampling period between september and november
  geom_rect(aes(xmin = 244, xmax = 305,ymin =-Inf, ymax = Inf), inherit.aes = FALSE, fill = "lightblue", alpha = 0.50)+
  # change color of ENSO years
  scale_color_manual(values = c(
    "1998" = "#f98e09",
    "2009" = "#bc3754",
    "2015" = "#57106e",
    "2023" = "#000004"), name = "ENSO Years")+
  # convert x scale to months
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+
  labs(y="NOAA SST (\u00B0C)")+
  theme_classic(base_size = 11)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
plot.crw_sst
    
#### ===> NOAA DHW plot <=== ####
plot.crw_dhw<-ggplot()+
  # all years except ENSO years
  geom_line(data = sstPix%>%filter(year != hl_years), 
            aes(x=doy, y=crw_dhw, group = year), color="gray80", alpha=0.9)+
  # plot ENSO years
  geom_line(data = sstPix%>%filter(year == hl_years),
            aes(x=doy, y=crw_dhw, group = year, color=factor(year)),linewidth = 0.9, alpha = 0.8)+
  # plot  2023 
  geom_line(data = sstPix%>%filter(year == 2023),
            aes(x=doy, y=crw_dhw, group = year, color=factor(year)), linewidth = 1 )+
  # highlight sampling period between september and november
  geom_rect(aes(xmin = 244, xmax = 305,ymin =-Inf, ymax = Inf), inherit.aes = FALSE, fill = "lightblue", alpha = 0.50)+
  # change color of ENSO years
  scale_color_manual(values = c(
    "1998" = "#f98e09",
    "2009" = "#bc3754",
    "2015" = "#57106e",
    "2023" = "#000004"), name = "ENSO Years")+
  # convert x scale to months
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+
  labs(y="NOAA DHW (\u00B0C weeks)")+
  theme_classic(base_size = 11)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
plot.crw_dhw

#### ===> Plot in-situ USWF nightly mean SST <=== ####
plot.uswf<-ggplot()+
  # plot UWSF all years: 2016 to 2025
  geom_line(data = uswf.night_sum,
            aes(x=doy, y=temp_c, group=year), color="gray50", alpha=0.9)+
  # isolate USWF 2023
  geom_line(data = uswf.night_sum%>%filter(year == 2023),
            aes(x=doy, y=temp_c, group = year), color="#000004", linewidth = 1)+
  # highlight sampling period between September and November
  geom_rect(aes(xmin = 244, xmax = 305,ymin =-Inf, ymax = Inf), inherit.aes = FALSE, fill = "lightblue", alpha = 0.50)+
  # convert x scale to months
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+
  labs(y="USWFS In situ night time SST (\u00B0C)")+
  theme_classic(base_size = 11)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
plot.uswf

#### ===> plot box plots for period with NOAA SST and insitu USWF <=== ####
  ### ==> Group situ with NOAA SST and USWF in situ data
    ## plot dates to isolate NOAA SST and USWF
    start_doy = 269
    end_doy = 280
    ## subset NOAA SST
    sstPix_sampl<-sstPix%>%
      filter(year == 2023, doy >= start_doy, doy <= end_doy)%>%
      select(date, doy, month, year, crw_sst)%>%
      rename(temp_c = crw_sst)%>%
      mutate(site="NOAA")
    ## subset USWF
    uswf.night_sum_sampl<-uswf.night_sum %>%
      filter(year == 2023, doy >= start_doy, doy <= end_doy)%>%
      mutate(site="USFWS")
    ## join frames
    sst_boxplot<-bind_rows(sstPix_sampl, uswf.night_sum_sampl, temp_night_sum)
    ## adjust site factor
    sst_boxplot$site<-factor(as.factor(sst_boxplot$site), levels = c("NOAA", "USFWS", "Offshore", "Midshore","Inshore"))
    ## check frame
    glimpse(sst_boxplot)
    
  ### ==> Plot to sampling period
  plot.intemp<-ggplot(data = sst_boxplot,
                      aes(x=site, y=temp_c, fill=site))+
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6)+
    scale_fill_manual(values = c("grey60","grey60","#8da0cb","#66c2a5","#fc8d59"))+
    labs(y="Temperature (\u00B0C)")+
    theme_classic(base_size = 11)+
    theme(axis.title.x = element_blank(),
          legend.position = "none")
  plot.intemp
  
#### ===> plot diel in situ temperature variation <=== ####
plot.diel23_temp<-ggplot(data = diel23%>%filter(variable=="temp"),
                         aes(x=tod, y=mean,  fill=site, color=site))+
  geom_ribbon(aes(ymin = mean-se, ymax = mean+se), alpha = 0.2, color=NA)+
  geom_line(linewidth = 0.9)+
  scale_x_continuous(breaks = seq(0, 24, 6),
                     labels = \(x) sprintf("%02d:00", x),
                     limits = c(0, 24), expand = c(0, 0))+
  labs(y="Temperature (\u00B0C)", x="Time of day")+
  scale_fill_manual(values = pal_cols)+
  scale_color_manual(values = pal_cols)+
  theme_classic(base_size = 11)
plot.diel23_temp

#### ===> Stitch plots <=== ####  
  ## plot layout
  layout.temp<-("
               AABB
               CCDD
               #EE#
               ")
  ## stitch plot
  plot.temp<-(
              plot.crw_sst+plot.crw_dhw+
              plot.uswf+plot.intemp+
              plot.diel23_temp
              )+
      plot_layout(design = layout.temp)+
      plot_annotation(tag_levels = "A")   
  plot.temp
  
  ## Save plot
  #ggsave("Feedback_Outputs/Temperature.pdf", 
  #       plot =plot.temp, width = 5, height = 6, units = "in", scale = 1.5, dpi = 600) 

#
#
#
# ========================================================================================
#
#                     ####  ~~~~  Tidal modulation of PAR  ~~~ ####
#
# ========================================================================================
#### ===> Join PAR and normalise data <=== ####
# ~ Note: Raw PAR is ~95% solar geometry, so any tidal signal is invisible until sun angle is divided out.
# ~ Note: par_rel = PAR / sin(solar elevation) is a relative transmission index: 
        # a flat plateau on clear days, with departures reflecting cloud, turbidity or depth
        # The sin_elev > 0.35 cut (~20 deg elevation) is essential:
        # near sunrise and sunset the ratio blows up and swamps everything
tide_par <- env23 %>%
  select(datetime, site, par) %>%
  ## force the same tz as the tide data; drop this line if already correct
  mutate(datetime = as.POSIXct(format(datetime), tz = "Pacific/Kiritimati")) %>%
  mutate(
    tide_h   = approx(tide_pal.interp$datetime, tide_pal.interp$height, xout = datetime)$y,
    dh       = approx(tide_pal.interp$datetime, tide_pal.interp$dh,     xout = datetime)$y,
    phase    = if_else(dh > 0, "Flood", "Ebb"),
    t_hour   = tidal_hour(datetime, hw),
    sin_elev = sin(getSunlightPosition(date = datetime, lat = lat, lon = lon)$altitude),
    ## par > 0 guard: some quantum sensors report a small negative offset in
    ## dim light, which gets amplified by the division.
    par_rel  = if_else(sin_elev > 0.35 & par > 0, par / sin_elev, NA_real_)
  ) %>%
  # tz argument matters - as.Date() on POSIXct converts via UTC by default  
  mutate(day = as.Date(datetime, tz = "Pacific/Kiritimati"))   

  ### ==> Check data for issues 
  # Solar geometry
  range(tide_par$sin_elev, na.rm = TRUE)
  tide_par %>% slice_max(sin_elev, n = 1) %>% select(datetime, sin_elev)
  # par_rel
  summary(tide_par$par_rel)
  # check daylight filter
  tide_par %>% filter(!is.na(par_rel)) %>%
    mutate(h = as.numeric(format(datetime, "%H"))) %>% count(h)
  # identifiably
    # ~ Note: Tidal phase advances ~50 min/day, so a short record risks aliasing tidal hour against time of day. 
    # Every cell should be populated. Empty cells => tide and time-of-day are confounded, 
    # and s(t_hour) would partly be fitting residual diurnal structure.
  dat_chk <- tide_par %>%
    filter(!is.na(par_rel), !is.na(t_hour)) %>%
    mutate(solar_h = as.numeric(format(datetime, "%H")) +
             as.numeric(format(datetime, "%M"))/60)
    with(dat_chk, table(cut(solar_h, seq(8, 18, 2)),
                        cut(t_hour,  seq(0, 12.42, 2.07))))

#### ===> Plot time series for tide and PAR with noon correction <=== ####
# ~ Note: Tide restricted to the same daylight window as PAR, so the eye doesn't
        # pair daytime PAR with night-time tidal peaks.
  ### ==> Isolate tide for day time
  tide_day <- tide_pal.interp %>%
    filter(datetime >= min(tide_par$datetime, na.rm = TRUE),
           datetime <= max(tide_par$datetime, na.rm = TRUE)) %>%
    mutate(sin_elev = sin(getSunlightPosition(date = datetime, lat = lat, lon = lon)$altitude)) %>%
    filter(sin_elev > 0.35) %>%
    mutate(day = as.Date(datetime, tz = "Pacific/Kiritimati"))
  
    ## Insert explicit NA rows at night gaps. geom_line() ALWAYS breaks at NA,
    # 15 min threshold sits safely above the 5-min sampling and far below the ~15 h night gap.
    insert_na <- function(d, yvar, gap_min = 15) {
      d <- arrange(d, datetime)
      gaps <- which(c(0, diff(as.numeric(d$datetime))/60) > gap_min)
      if (!length(gaps)) return(d)
      pad <- d[gaps, ]
      pad$datetime  <- pad$datetime - 60
      pad[[yvar]]   <- NA
      bind_rows(d, pad) %>% arrange(datetime)
    }
    
    tide_day_plot <- insert_na(tide_day, "height")
    
    # group_by(site) matters
    tide_par_plot <- tide_par %>%
      filter(!is.na(par_rel)) %>%
      arrange(site, datetime) %>%
      group_by(site) %>%
      group_modify(~ insert_na(.x, "par_rel")) %>%
      ungroup()
    
  ### ==> Rescale tide onto the par_rel axis for the secondary axis
  par_rng  <- range(tide_par_plot$par_rel, na.rm = TRUE)
  tide_rng <- range(tide_day_plot$height,  na.rm = TRUE)
  sf       <- diff(par_rng) / diff(tide_rng)
  shift    <- par_rng[1] - tide_rng[1] * sf
  
  ### ==> plot
  plot.tidePAR <- ggplot() +
    geom_line(data = tide_par_plot,
              aes(datetime, par_rel, colour = site),
              linewidth = 0.5, show.legend = FALSE) +
    geom_line(data = tide_day_plot,
              aes(datetime, height * sf + shift),
              colour = "grey10", linewidth = 0.7) +
    scale_y_continuous(name = expression(PAR/sin(theta[sun])),
                       sec.axis = sec_axis(~ (. - shift) / sf,
                                           name = "Predicted tide (m)")) +
    xlab("Date")+
    scale_x_datetime(date_breaks = "3 days", date_labels = "%d %b") +
    scale_colour_manual(values = pal_cols) +
    facet_wrap(~ site, ncol = 1) +
    theme_classic(base_size = 11)
  plot.tidePAR
  
    ## Save plot
    #ggsave("Feedback_Outputs/tidePAR_noon.png", 
    #       plot=plot.tidePAR, 
    #       width = 5, height = 5, units = "in", scale = 1, dpi = 600) 
    
#### ===> Plot tidal-hour <=== ####
# ~ Note: Collapses the time axis. If tide mattered you would see a cyclic smooth at
  # the inshore site and a flat one offshore
  
plot.tidePAR_HR<-tide_par %>%
    filter(!is.na(par_rel), !is.na(t_hour)) %>%
    ggplot(aes(t_hour, par_rel, colour = site)) +
    geom_point(alpha = 0.15, size = 0.4) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = 6)) +
    facet_wrap(~ site, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = pal_cols) +
    labs(x = "Hours since high water", y = "Normalised PAR") +
    theme_classic(base_size = 11)+
    theme(legend.position = "none")
plot.tidePAR_HR

  ## Save plot
  #ggsave("Feedback_Outputs/tidePAR_hour.png", 
  #       plot=plot.tidePAR_HR, 
  #       width = 5, height = 5, units = "in", scale = 1, dpi = 600) 

#### ===> Model relationship between tide and PAR <=== ####
  ### ==> Set up data frame
  dat <- tide_par %>%
    filter(!is.na(par_rel), !is.na(t_hour)) %>%
    arrange(site, datetime) %>%
    mutate(t_num = as.numeric(datetime), site = factor(site)) %>%
    ## AR.start flag: resets the AR1 process at each site-day boundary, so the
    ## model does not treat dusk as adjacent to the following dawn.
    group_by(site, day) %>% 
    mutate(new_day = row_number() == 1) %>% 
    ungroup()

  ### ==> Estimate residual autocorrelation from an UNCORRECTED fit first.
  m0 <- bam(par_rel ~ site + s(t_hour, by = site, bs = "cc", k = 6) + s(t_num, k = 15),
            data = dat, knots = list(t_hour = c(0, 12.42)), method = "fREML")
  
  acf(resid(m0), lag.max = 40)
  r <- acf(resid(m0), plot = FALSE)$acf[2]   # lag-1; ~0.70 at 5-min sampling
  r

  ### ==> TIDAL PHASE, site-specific.
  # ~ Note: Check the EDF, not just the p-value: 
          # EDF ~1 means the smooth has been shrunk to a straight line
          # A real semi-diurnal signal needs EDF ~3-4.
  m1 <- bam(par_rel ~ site + s(t_hour, by = site, bs = "cc", k = 6) + s(t_num, k = 15),
            data = dat, knots = list(t_hour = c(0, 12.42)),
            rho = r, AR.start = dat$new_day, method = "fREML")
  summary(m1)
  plot(m1, pages = 1, scale = 0, shade = TRUE)

  ### ==> TIDAL HEIGHT, with a SENSITIVITY SWEEP over the temporal smooth.
  # ~ Npte: s(t_num) and the tidal terms compete for the same slow variance
          # any single k is a choice that drives the answer. Should report the sweep, not one model.
  fit_height <- function(k_time) {
    bam(par_rel ~ site + s(tide_h, by = site, k = 5) + s(t_num, k = k_time),
        data = dat, rho = r, AR.start = dat$new_day, method = "fREML")
  }
  
  m2c <- fit_height(10)    # stiff  -> tide terms significant
  m2  <- fit_height(25)    # medium -> mostly null
  m2b <- fit_height(40)    # flexible, AIC-preferred -> all null
  
  summary(m2c); summary(m2); summary(m2b)
  AIC(m2c, m2, m2b)
   # ~ Note: AIC from bam() with rho is approximate (effective n < n). Treat differences under ~4 units as ties.
  
  ### ==>  NESTED COMPARISON at MATCHED temporal flexibility.
  # ~ Note: All three must share the same s(t_num, k) or the comparison confounds the
          # tidal term with the weather term.
  m1b <- bam(par_rel ~ site + s(t_hour, by = site, bs = "cc", k = 6) + s(t_num, k = 25),
             data = dat, knots = list(t_hour = c(0, 12.42)),
             rho = r, AR.start = dat$new_day, method = "fREML")
  
  m3  <- bam(par_rel ~ site + s(t_hour, bs = "cc", k = 6) + s(t_num, k = 25),
             data = dat, knots = list(t_hour = c(0, 12.42)),
             rho = r, AR.start = dat$new_day, method = "fREML")
  
  m4  <- bam(par_rel ~ site + s(t_num, k = 25),
             data = dat, rho = r, AR.start = dat$new_day, method = "fREML")
  
  AIC(m1b, m3, m4)
  # m3 and m4 identical to 2 dp, df identical to 3 dp -> 
  # the tidal smooth has been penalised to ZERO effective degrees of freedom. 
  # The model with a tidal term IS the model without one.


  ### ==> Basis-dimension check.
  gam.check(m2b)
  # ~ Note: with AR1 correction and 5-min data the k-index is unreliable - it
  # flags residual autocorrelation, not basis dimension. A term with EDF 1.0
  # against k'=4 cannot be under-parameterised.Do not chase these by raising k.
  
  ### ==> check the site gradient
  dat %>%
    filter(day > as.Date("2023-09-27"), day < as.Date("2023-10-06")) %>%
    group_by(site) %>%
    summarise(mean_rel   = mean(par_rel),
              median_rel = median(par_rel),
              n = n())
  # Observed: 798 / 679 / 524, monotonic, equal n, mean ~ median.
  # Close agreement with model coefficients) 
  # confirms the site effect is not an artefact of the smooths
  # but equally means it carries no more protection against confounding variables

# ========================================================================================
#
#                   ####  ~~~~  PLOT: Environmental patterns 2023  ~~~ ####
#
# ========================================================================================
#### ===> Box plots <=== ####
  ### ==> PAR
  plot_box.par<-ggplot(data = env23%>%filter(hms(format(datetime, "%H:%M:%S")) >= sunrise &
                                                hms(format(datetime, "%H:%M:%S")) <= sunset),
                        aes(x=site, y=par, fill=site, color=site))+ 
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6)+
    labs(y=expression(PAR~(µmol~m^-2~s^-1)), x="")+
    scale_fill_manual(values = pal_cols)+
    scale_color_manual(values = pal_cols)+
    theme_classic(base_size = 11)+
    theme(legend.position = "none")
  plot_box.par
  
  ### ==> AOU
  plot_box.aou<-ggplot(data = env23,
                        aes(x=site, y=aou, fill=site, color=site))+ #fct_rev(site)
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6)+
    labs(y=expression(AOU~(mg~l^-1)), x="")+
    scale_fill_manual(values = pal_cols)+
    scale_color_manual(values = pal_cols)+
    theme_classic(base_size = 11)+
    theme(legend.position = "none")
  plot_box.aou
  
  ### ==> Speed
  plot_box.speed<-ggplot(data = env23,
                        aes(x=site, y=speed, fill=site, color=site))+ #fct_rev(site)
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6)+
    labs(y=expression(Speed~(cm~s^-1)), x="")+
    scale_fill_manual(values = pal_cols)+
    scale_color_manual(values = pal_cols)+
    theme_classic(base_size = 11)+
    theme(legend.position = "none")
  plot_box.speed

#### ===> 24 hour variation <=== #### 
  ### ==> PAR
  plot.diel23_par<-ggplot(data = diel23%>%filter(variable=="par"),
                           aes(x=tod, y=mean,  fill=site, color=site))+
    geom_ribbon(aes(ymin = mean-se, ymax = mean+se), alpha = 0.2, color=NA)+
    geom_line(linewidth = 0.9)+
    scale_x_continuous(breaks = seq(0, 24, 6),
                       labels = \(x) sprintf("%02d:00", x),
                       limits = c(0, 24), expand = c(0, 0))+
    labs(y=expression(PAR~(µmol~m^-2~s^-1)), x="Time of day")+
    scale_fill_manual(values = pal_cols)+
    scale_color_manual(values = pal_cols)+
    theme_classic(base_size = 11)+
    theme(legend.position = "none")
  plot.diel23_par
  ### ==> AOU
  plot.diel23_aou<-ggplot(data = diel23%>%filter(variable=="aou"),
                           aes(x=tod, y=mean,  fill=site, color=site))+
    geom_ribbon(aes(ymin = mean-se, ymax = mean+se), alpha = 0.2, color=NA)+
    geom_line(linewidth = 0.9)+
    scale_x_continuous(breaks = seq(0, 24, 6),
                       labels = \(x) sprintf("%02d:00", x),
                       limits = c(0, 24), expand = c(0, 0))+
    labs(y=expression(AOU~(mg~l^-1)), x="Time of day")+
    scale_fill_manual(values = pal_cols)+
    scale_color_manual(values = pal_cols)+
    theme_classic(base_size = 11)+
    theme(legend.position = "none")
  plot.diel23_aou
  ### ==> Speed
  plot.diel23_speed<-ggplot(data = diel23%>%filter(variable=="speed"),
                           aes(x=tod, y=mean,  fill=site, color=site))+
    geom_ribbon(aes(ymin = mean-se, ymax = mean+se), alpha = 0.2, color=NA)+
    geom_line(linewidth = 0.9)+
    scale_x_continuous(breaks = seq(0, 24, 6),
                       labels = \(x) sprintf("%02d:00", x),
                       limits = c(0, 24), expand = c(0, 0))+
    labs(y=expression(Speed~(cm~s^-1)), x="Time of day")+
    scale_fill_manual(values = pal_cols)+
    scale_color_manual(values = pal_cols)+
    theme_classic(base_size = 11)+
    theme(legend.position = "none")
  plot.diel23_speed
  
#### ===> Rose plot for summarised velocity <=== ####   
  ### ==> Create rose
  nbin <- 16                 # 16 compass sectors
  w    <- 360 / nbin         # 22.5° each
  
  rose <- env23 %>%
    filter(!is.na(u), !is.na(v)) %>%
    mutate(
      spd     = speed,
      dir     = (atan2(u, v) * 180 / pi) %% 360,
      dir_bin = (round(dir / w) * w) %% 360,        # 355° -> 360 -> 0
      spd_cls = cut(spd, breaks = c(0, 2, 4, 6, 8, 10, Inf),
                    labels = c("<2", "2–4", "4–6", "6–8", "8–10", "≥10"),
                    right = FALSE)
    ) %>%
    summarise(n = n(), .by = c(site, dir_bin, spd_cls)) %>%
    mutate(pct = 100 * n / sum(n), .by = site)      # % within each site
  
  ### ==> Plot rose plot
  plot.rose23<-ggplot(rose, aes(dir_bin, pct, fill = spd_cls)) +
    geom_col(width = w, colour = "white", linewidth = 0.15) +
    coord_polar(start = -w * pi / 180 / 2) +
    scale_x_continuous(
      limits = c(-w/2, 360 - w/2),
      breaks = seq(0, 315, 45),
      labels = c("N","NE","E","SE","S","SW","W","NW")
    ) +
    scale_fill_viridis_d(option = "mako", direction = -1, name = "cm/s") +
    facet_wrap(~ site) +
    labs(x = NULL, y = "% of observations") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom")+
    guides(fill = guide_legend(nrow = 1))
  plot.rose23

#### ===> Stitch plots together <=== ####  
  ## layout
  layout_env23<- "
                  ABB
                  CDD
                  EFF
                  GGG
                 "

  ## combine plots
  plot.env<-(
     plot_box.par+plot.diel23_par+
     plot_box.speed+plot.diel23_speed+
     plot_box.aou+plot.diel23_aou+
     plot.rose23)+
    plot_layout(design = layout_env23)+
    plot_annotation(tag_levels = "A")
  plot.env
  
  ## Save plot
  ggsave("Feedback_Outputs/Environment2023.pdf", 
         plot=plot.env, 
         width = 5, height = 6, units = "in", scale = 1.5, dpi = 600) 
  
# ========================================================================================
#
#                   ####  ~~~~  PLOT: Environmental patterns 2024  ~~~ ####
#
# ========================================================================================  
#### ===> Box plot <=== ####
plot_box.speed24<-ggplot(data = env24,
                       aes(x=site, y=speed, fill=site, color=site))+ #fct_rev(site)
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6)+
  labs(y=expression(Speed~(cm~s^-1)), x="")+
  scale_fill_manual(values = pal_cols)+
  scale_color_manual(values = pal_cols)+
  theme_classic(base_size = 11)+
  theme(legend.position = "none")
plot_box.speed24

#### ===> 24 hour variation <=== ####
plot.diel24_speed<-ggplot(data = diel24%>%filter(variable=="speed"),
                          aes(x=tod, y=mean,  fill=site, color=site))+
  geom_ribbon(aes(ymin = mean-se, ymax = mean+se), alpha = 0.2, color=NA)+
  geom_line(linewidth = 0.9)+
  scale_x_continuous(breaks = seq(0, 24, 6),
                     labels = \(x) sprintf("%02d:00", x),
                     limits = c(0, 24), expand = c(0, 0))+
  labs(y=expression(Speed~(cm~s^-1)), x="Time of day")+
  scale_fill_manual(values = pal_cols)+
  scale_color_manual(values = pal_cols)+
  theme_classic(base_size = 11)+
  theme(legend.position = "none")
plot.diel24_speed

#### ===> Rose plot for summarised velocity <=== ####
  ### ==> create rose
  nbin <- 16                 # 16 compass sectors
  w    <- 360 / nbin         # 22.5° each

  rose24 <- env24 %>%
    filter(!is.na(u), !is.na(v)) %>%
    mutate(
      spd     = speed,
      dir     = (atan2(u, v) * 180 / pi) %% 360,
      dir_bin = (round(dir / w) * w) %% 360,        # 355° -> 360 -> 0
      spd_cls = cut(spd, breaks = c(0, 2, 4, 6, 8, 10, Inf),
                    labels = c("<2", "2–4", "4–6", "6–8", "8–10", "≥10"),
                    right = FALSE)
    ) %>%
    summarise(n = n(), .by = c(site, dir_bin, spd_cls)) %>%
    mutate(pct = 100 * n / sum(n), .by = site)      # % within each site
  
  ### ==> Plot rose plot
  plot.rose24<-ggplot(rose24, aes(dir_bin, pct, fill = spd_cls)) +
    geom_col(width = w, colour = "white", linewidth = 0.15) +
    coord_polar(start = -w * pi / 180 / 2) +
    scale_x_continuous(
      limits = c(-w/2, 360 - w/2),
      breaks = seq(0, 315, 45),
      labels = c("N","NE","E","SE","S","SW","W","NW")
    ) +
    scale_fill_viridis_d(option = "mako", direction = -1, name = "cm/s") +
    facet_wrap(~ site) +
    labs(x = NULL, y = "% of observations") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom")+
    guides(fill = guide_legend(nrow = 1))
  plot.rose24
  
#### ===> Stitch plots together <=== ####
  ## layout
  layout_env24<- "
                  ABB
                  CCC
                 "
  
  ## combine plots
  plot.env24<-(
      plot_box.speed24+plot.diel24_speed+
      plot.rose24)+
    plot_layout(design = layout_env24)+
    plot_annotation(tag_levels = "A")
  plot.env24
  
  ## Save plot
  #ggsave("Feedback_Outputs/Environment2024.pdf", 
  #       plot=plot.env24, 
  #       width = 5, height = 5, units = "in", scale = 1.5, dpi = 600) 
  
# ========================================================================================
#
#           ####  ~~~~  Testing spatial difference individual parameters  ~~~ ####
#
# ======================================================================================== 
set.seed(1984)
#### ===> Summarise data to daily averages to reduce noise <=== ####
  # run if running section alone -> same as next section for testing 
  ### ===> 2023 
  env23.daily<-env23%>%
    select(date, site, par, temp, aou, speed, u, v)%>%
    pivot_longer(cols = par:v,
                 names_to = "parameter",
                 values_to = "measurement")%>%
    group_by(date, site, parameter)%>%
    summarise(mean = mean(measurement),
              sd = sd(measurement))%>%
    ungroup() %>% 
    pivot_wider(names_from = parameter, 
                values_from = c(mean, sd)) 
  ### ===> 2024 
  env24.daily<-env23%>%
    select(date, site, speed, u, v)%>%
    pivot_longer(cols = speed:v,
                 names_to = "parameter",
                 values_to = "measurement")%>%
    group_by(date, site, parameter)%>%
    summarise(mean = mean(measurement),
              sd = sd(measurement))%>%
    ungroup() %>% 
    pivot_wider(names_from = parameter, 
                values_from = c(mean, sd))  
  
#### ===> set up data for analysis <=== ####
  ### ==> meta data
  meta_DLI<-select(DLI, site)
  meta_23<-select(env23.daily, site)
  meta_24<-select(env24.daily, site)

  ### ==> distance matrix
  dis_DLI<-vegdist(DLI$DLI, method = "euclidean")
  dis_par<-vegdist(env23.daily$mean_par, method = "euclidean")
  dis_temp<-vegdist(env23.daily$mean_temp, method = "euclidean")
  dis_aou<-vegdist(env23.daily$mean_aou, method = "euclidean")
  dis_speed<-vegdist(env23.daily$mean_speed, method = "euclidean")
  dis_u<-vegdist(env23.daily$mean_u, method = "euclidean")
  dis_v<-vegdist(env23.daily$mean_v, method = "euclidean")
  dis_speed24<-vegdist(env24.daily$mean_speed, method = "euclidean")
  dis_u24<-vegdist(env24.daily$mean_u, method = "euclidean")
  dis_v24<-vegdist(env24.daily$mean_v, method = "euclidean")
  
#### ===> DLI <=== ####
per_DLI<-adonis2(dis_DLI~site, data = meta_DLI, permutations = 9999, by = "margin")
per_DLI
  ## Pair-wise
  meta_DLI$sites<-meta_DLI$site
  per_DLI.pw<-pairwise.adonis2(dis_DLI~sites, data = meta_DLI, nperm = 999)
  per_DLI.pw

#### ===> PAR <=== ####
per_par<-adonis2(dis_par~site, data = meta_23, permutations = 9999, by = "margin")
per_par
  ## Pair-wise
  per_par.pw<-pairwise.adonis2(dis_par~site, data = meta_23, nperm = 999)
  per_par.pw
  
#### ===> Temperature <=== ####
per_temp<-adonis2(dis_temp~site, data = meta_23, permutations = 9999, by = "margin")
per_temp
  
#### ===> AOU <=== ####
per_aou<-adonis2(dis_aou~site, data = meta_23, permutations = 9999, by = "margin")
per_aou
  
#### ===> Speed 2023 <=== ####
per_speed<-adonis2(dis_speed~site, data = meta_23, permutations = 9999, by = "margin")
per_speed
  ## Pair-wise
  per_speed.pw<-pairwise.adonis2(dis_speed~site, data = meta_23, nperm = 999)
  per_speed.pw
  
#### ===> u 2023 <=== ####
per_u<-adonis2(dis_u~site, data = meta_23, permutations = 9999, by = "margin")
per_u
  ## Pair-wise
  per_u.pw<-pairwise.adonis2(dis_u~site, data = meta_23, nperm = 999)
  per_u.pw
  
#### ===> v 2023 <=== ####
per_v<-adonis2(dis_v~site, data = meta_23, permutations = 9999, by = "margin")
per_v
  ## Pair-wise
  per_v.pw<-pairwise.adonis2(dis_v~site, data = meta_23, nperm = 999)
  per_v.pw
  
#### ===> Speed 2024 <=== ####
per_speed24<-adonis2(dis_speed24~site, data = meta_24, permutations = 9999, by = "margin")
per_speed24
  ## Pair-wise
  per_speed24.pw<-pairwise.adonis2(dis_speed24~site, data = meta_24, nperm = 999)
  per_speed24.pw
  
#### ===> u 2024 <=== ####
per_u24<-adonis2(dis_u24~site, data = meta_24, permutations = 9999, by = "margin")
per_u24
  ## Pair-wise
  per_u24.pw<-pairwise.adonis2(dis_u24~site, data = meta_24, nperm = 999)
  per_u24.pw
  
#### ===> v 2024 <=== ####
per_v24<-adonis2(dis_v24~site, data = meta_24, permutations = 9999, by = "margin")
per_v24
  ## Pair-wise
  per_v24.pw<-pairwise.adonis2(dis_v24~site, data = meta_24, nperm = 999)
  per_v24.pw
  
# ========================================================================================
#
#                   ####  ~~~~  RDA of 2023 Environmental data  ~~~ ####
#
# ========================================================================================   
set.seed(1984)
#### ===> Summarise data to daily averages to reduce noise <=== ####
# run if running section alone -> same as previous section for testing 
env23.daily<-env23%>%
    select(date, site, par, temp, aou, speed, u, v)%>%
    pivot_longer(cols = par:v,
                 names_to = "parameter",
                 values_to = "measurement")%>%
    group_by(date, site, parameter)%>%
    summarise(mean = mean(measurement),
              sd = sd(measurement))%>%
    ungroup() %>% 
    pivot_wider(names_from = parameter, 
              values_from = c(mean, sd))
  
#### ===> set up data for the RDA <=== ####
  ### ==> Set up data and meta data 
  env.meta<-select(env23.daily, site, date)
  env.meta <- env.meta %>% mutate(date_f = factor(date))
  env.data<-select(env23.daily, -site, -date)
  
  ### ==> check mean-SD coupling
  env.pairs <- c("aou","temp","u","v","speed","par")
  map_dfr(env.pairs, \(p) tibble(
    var = p,
    r = cor(env.data[[paste0("mean_", p)]], env.data[[paste0("sd_", p)]],
            use = "pairwise.complete.obs")))
        # dropping both par sd and mean speed -> both |r| > ~0.7
        env.data_red<-select(env.data, -mean_speed, -sd_par)
    
  ### ===> screen for skew predictors          
  summarise(env.data_red, across(everything(), \(x) sd(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE))))
  # >3 = single-point leverage
  map_dbl(env.data_red, \(x) max(abs(scale(x)), na.rm = TRUE))  

  ### ==> check redundancy
  d <- as.dist(1 - abs(cor(env.data_red, use = "pairwise.complete.obs")))
  plot(hclust(d), main = "predictor redundancy")
  abline(h = 0.3, lty = 2)   # cuts at |r| = 0.7
    ## drop mean aou, aou sd, speed sd
    env.data_red<-select(env.data_red, -mean_aou, -sd_aou, -sd_speed)
    
  ### ==> Standardise data and convert to a dissimilarity matrix
    ## full data set
    env.data.std<-decostand(env.data, method = "standardize")
    ## reduced data set
    env.data_red.std<-decostand(env.data_red, method = "standardize")
    
#### ===> Test within group dispersion <=== ####
  ### ==> Full data set
  beta_full<-betadisper(vegdist(env.data.std, method = "euclidean"), group = env.meta$site, type = "centroid") #using centroid to match with RDA geometry (LS)
  beta_full
    ## test
    permutest(beta_full, permutations = how(nperm = 9999, blocks = env.meta$date))
    tapply(beta_full$distances, env.meta$site, mean)
    
  ### ==> Reduced data set
  beta_red<-betadisper(vegdist(env.data_red.std, method = "euclidean"), group = env.meta$site, type = "centroid") #using centroid to match with RDA geometry (LS)
  beta_red
    ## test
    permutest(beta_red, permutations = how(nperm = 9999, blocks = env.meta$date))
    tapply(beta_red$distances, env.meta$site, mean)
  
#### ===> RDA on the full data set <=== ####
  ### ==> Run RDA     
  RDA_full<-rda(env.data.std~site, data=env.meta)  
    ## check significance of the RDA
    anova(RDA_full, permutations = how(nperm = 9999, blocks = env.meta$date)) # global
    anova(RDA_full, by = "axis", permutations = how(nperm = 9999, blocks = env.meta$date)) # per axis
    #anova(RDA_full, by = "term", permutations = 999) # per term
    ## check adjusted rsquared
    RsquareAdj(RDA_full)
    ## check the axes results of the RDA
    summary(RDA_full)
      # RDA 1 explains: 47.7%
      # RDA 2 explains: 5.1%
    
  ### ==> Test which variable is absorbing the residuals in PC1: temporal dates vs spatial sites
  RDA_full_both <- rda(env.data.std ~ site + Condition(date_f), data = env.meta)
    summary(RDA_full_both)
      # RDA 1 explains: 71.7%
      # RDA 2 explains: 7.7%
  RDA_full_date <- rda(env.data.std ~ date_f + Condition(site), data = env.meta)
    summary(RDA_full_date)
      # PC1 was a day effect -> so use RDA_full_both for plotting
    ## run ANOVAs
    anova(RDA_full_date, permutations = 9999)
    anova(RDA_full_both, by = "axis", permutations = 9999)
    RsquareAdj(RDA_full_both)     # for the site effect adjusted for date
    ## pairwise testing to see which sites differ
    # pairwise RDAs, conditioning on date, with p-adjustment
    pairs <- combn(levels(env.meta$site), 2, simplify = FALSE)
    map_dfr(pairs, \(p) {
      keep <- env.meta$site %in% p
      m <- rda(env.data.std[keep, ] ~ droplevels(site) + Condition(droplevels(date_f)),
               data = droplevels(env.meta[keep, ]))
      a <- anova(m, permutations = 9999)
      tibble(pair = paste(p, collapse = " vs "), F = a$F[1], p = a$`Pr(>F)`[1])
    }) %>% mutate(p_adj = p.adjust(p, "holm"))
    
  ### ==> Check variance attributed to site and date
  vp_full <- varpart(env.data.std, ~ site, ~ date_f, data = env.meta)
  plot(vp_full); vp_full

  ### ==> Isolate constraints for plotting 
    ## isolate constrained (CCA) u scores -> site score: ordination of rows
    rda_full_u.scores<-data.frame(RDA_full_both$CCA$u)
      row.names(rda_full_u.scores)<-c(1:36)
      rda_full_u.scores.j<- inner_join(rownames_to_column(env.meta), 
                                       rownames_to_column(data.frame(rda_full_u.scores)), 
                                       by = "rowname")
    ## isolate constrained (CCA) was scores -> weighted average site score
    rda_full_wa.scores<-data.frame(RDA_full_both$CCA$wa)
      row.names(rda_full_wa.scores)<-c(1:36)
      rda_full_wa.scores.j<- inner_join(rownames_to_column(env.meta), 
                                       rownames_to_column(data.frame(rda_full_wa.scores)), 
                                       by = "rowname")
    ## isolate constrained (CCA) v scores -> species score: coordinates of response
    rda_full_v.scores<-data.frame(RDA_full_both$CCA$v)
      # check the importance of each set of scores
      goodness(RDA_full_both, display = "species", model = "CCA", summ = TRUE)
    ## check site centroids
    RDA_full_both$CCA$centroids
    
  ### ==> plot RDA
  plot.rda_full<-ggplot()+
    geom_segment(data = rda_full_v.scores, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
                 arrow=arrow(length=unit(0.3,"cm")), alpha = 0.5, linewidth = 1, color = 'black')+
    stat_ellipse(data=rda_full_wa.scores.j, aes(x=RDA1, y=RDA2, color=site),
                 linewidth=1, level = 0.68, type="norm")+
    geom_point(data=rda_full_wa.scores.j, aes(x=RDA1, y=RDA2, fill=site),
               shape=21, size=3)+
    geom_text(data = rda_full_v.scores, 
              aes(x = RDA1, y = RDA2, label = rownames(rda_full_v.scores)), col = 'red') +
    scale_fill_manual(values = c(pal_cols))+
    scale_color_manual(values = c(pal_cols))+
    labs(x="RDA1 (71.6%)",
         y="RDA2 (7.7%)",
         color="",fill="")+
    coord_fixed() +
    theme_classic(base_size = 10)+
    theme(legend.position = "bottom")+
    guides(fill = guide_legend(nrow = 1))
  plot.rda_full
    ## Save plot
    #ggsave("Feedback_Outputs/RDA_full.pdf", 
    #       plot=plot.rda_full, 
    #       width = 3.5, height = 3.5, units = "in", scale = 1.5, dpi = 600) 
    
    
#### ===> RDA on the reduced data set <=== ####
  ### ==> Run RDA     
  RDA_red<-rda(env.data_red.std~site, data=env.meta)  
    ## check significance of the RDA
    anova(RDA_red, permutations = how(nperm = 999, blocks = env.meta$date)) # global
    anova(RDA_red, by = "axis", permutations = how(nperm = 999, blocks = env.meta$date)) # per axis
    #anova(RDA_full, by = "term", permutations = 999) # per term
    ## check adjusted RDA
    RsquareAdj(RDA_red)
    ## check the axes results of the RDA
    summary(RDA_red)
    # RDA 1 explains: 45.0%
    # RDA 2 explains: 8.4%%
    
  ### ==> Test which variable is absorbing the residuals in PC1: temporal dates vs spatial sites
  RDA_red_both <- rda(env.data_red.std ~ site + Condition(date_f), data = env.meta)
    summary(RDA_red_both)
    # RDA 1 explains: 66.4%
    # RDA 2 explains: 12.4%
  RDA_red_date <- rda(env.data_red.std ~ date_f + Condition(site), data = env.meta)
    summary(RDA_red_date)
    # PC1 was a day effect -> so use RDA_red_both for plotting
    ## run ANOVAs
    anova(RDA_red_date, permutations = 9999)
    anova(RDA_red_both, permutations = how(nperm = 9999, blocks = env.meta$date_f))
    anova(RDA_red_both, by = "axis", permutations = 9999)
    RsquareAdj(RDA_red_both)     # for the site effect adjusted for date
    ## pairwise testing
    pairs <- combn(levels(env.meta$site), 2, simplify = FALSE)
    map_dfr(pairs, \(p) {
      keep <- env.meta$site %in% p
      m <- rda(env.data_red.std[keep, ] ~ droplevels(site) + Condition(droplevels(date_f)),
               data = droplevels(env.meta[keep, ]))
      a <- anova(m, permutations = 9999)
      tibble(pair = paste(p, collapse = " vs "), F = a$F[1], p = a$`Pr(>F)`[1])
    }) %>% mutate(p_adj = p.adjust(p, "holm"))
    
  ### ==> Check variance attributed to site and date
  vp_red <- varpart(env.data_red.std, ~ site, ~ date_f, data = env.meta)
  plot(vp_red); vp_red
    
  ### ==> Isolate constraints for plotting 
  ## isolate constrained (CCA) u scores -> site score: ordination of rows
  rda_red_u.scores<-data.frame(RDA_red_both$CCA$u)
    row.names(rda_red_u.scores)<-c(1:36)
    rda_red_u.scores.j<- inner_join(rownames_to_column(env.meta), 
                                     rownames_to_column(data.frame(rda_red_u.scores)), 
                                     by = "rowname")
  ## isolate constrained (CCA) was scores -> weighted average site score
  rda_red_wa.scores<-data.frame(RDA_red_both$CCA$wa)
    row.names(rda_red_wa.scores)<-c(1:36)
    rda_red_wa.scores.j<- inner_join(rownames_to_column(env.meta), 
                                      rownames_to_column(data.frame(rda_red_wa.scores)), 
                                      by = "rowname")
  ## isolate constrained (CCA) v scores -> species score: coordinates of response
  rda_red_v.scores<-data.frame(RDA_red_both$CCA$v)
    # check the importance of each set of scores
    goodness(RDA_red_both, display = "species", model = "CCA", summ = TRUE)
  ## check site centroids
  RDA_red_both$CCA$centroids
  
  ### ==> plot RDA
  plot.rda_red<-ggplot()+
    geom_segment(data = rda_red_v.scores, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
                 arrow=arrow(length=unit(0.3,"cm")), alpha = 0.5, linewidth = 1, color = 'black')+
    stat_ellipse(data=rda_red_wa.scores.j, aes(x=RDA1, y=RDA2, color=site),
                 linewidth=1, level = 0.68, type="norm")+
    geom_point(data=rda_red_wa.scores.j, aes(x=RDA1, y=RDA2, fill=site),
               shape=21, size=3)+
    geom_text(data = rda_red_v.scores, 
              aes(x = RDA1, y = RDA2, label = rownames(rda_red_v.scores)), col = 'red') +
    scale_fill_manual(values = c(pal_cols))+
    scale_color_manual(values = c(pal_cols))+
    labs(x="RDA1 (66.4%)",
         y="RDA2 (12.4%)",
         color="",fill="")+
    coord_fixed() +
    theme_classic(base_size = 10)+
    theme(legend.position = "bottom")+
    guides(fill = guide_legend(nrow = 1))
  plot.rda_red
    ## Save plot
    #ggsave("Feedback_Outputs/RDA_red.pdf", 
    #       plot=plot.rda_red, 
    #       width = 3.5, height = 3.5, units = "in", scale = 1.5, dpi = 600) 
    
#### ===> stitch both models together <=== ####
  ### ==> stitch together
  plot.rda<-(plot.rda_full+plot.rda_red)+
    plot_annotation(tag_levels = "A")
  plot.rda
  
  ## Save plot
  #ggsave("Feedback_Outputs/RDA_both.pdf", 
  #       plot=plot.rda, 
  #      width = 5, height = 3.5, units = "in", scale = 1.5, dpi = 600) 
  
# ========================================================================================
#
#       ####  ~~~~  Benthic cover data from 2023 and 2024  ~~~ ####
#
# ========================================================================================
#
# ~ Note: 2023 = pre-bleaching amd 2024 = post-bleaching
#
#### ===> Read percent cover data <=== ####
percent_w <- read.csv("Data/Cover_2324.csv",stringsAsFactors = TRUE) %>%
    clean_names()%>%
    # remove species columns with no counts -> 79 to 48 variables
    select_if(negate(function(col) is.numeric(col) && sum(col) == 0.000))%>%
    mutate(
      # reform site to show inshore, midshore and offshore
      site=factor(site, levels = c("RT1","RT4","RT7"),
                        labels = c("Offshore","Midshore","Inshore")),
      # convert date into a date variable
      date_dtm=as.POSIXct(date, format="%d/%m/%Y"),
      # extract year and convert to bleaching status
      year=factor(as.factor(year(date_dtm)),levels = c("2023","2024"), 
                                            labels = c("Pre", "Post")),
      # form factor combing site and year 
      site_state=as.factor(paste(site, year, sep = "."))
    )
  ## check data
  glimpse(percent_w)
  
#### ===> convert to long format data <=== ####  
percent_l<-percent_w%>%
    gather(benthos, cover, acr_bra:turf_rubbl, factor_key = TRUE)%>%
    mutate(ben_group = case_when(benthos=="acr_bra"|benthos=="acr_corym"|benthos=="acr_dig"|benthos=="acr_sub"|benthos=="acr_tab"|benthos=="ast_plat"|benthos=="astreo_enc"|benthos=="b_monti"|
                                 benthos=="dipmas"|benthos=="fav_mass"|benthos=="fungia"|benthos=="lpsenc"|benthos=="monmas"|benthos=="mon_pla"|benthos=="monti_encr"|
                                 benthos=="monti_fol"|benthos=="mont_sub_m"|benthos=="pavenc"|benthos=="pavfol"|benthos=="pavmas"|benthos=="pav_sub"|benthos=="plamas"|
                                 benthos=="pocill"|benthos=="por_mass"|benthos=="por_sm"|benthos=="potenc"|benthos=="psammo"|benthos=="styloph" ~ "Hard coral",
                                 benthos=="collimorph"|benthos=="ensp" ~ "Other",
                                 benthos=="sarco" ~ "Soft coral",
                                 benthos=="sand" ~ "Soft substrate",
                                 benthos=="bran_calc"|benthos=="peysson" ~ "Calcifying algae",
                                 benthos=="cca"|benthos=="cca_rub" ~ "CCA", 
                                 benthos=="caulerpa"|benthos=="dictyosph"|benthos=="mal" ~"Macro algae",
                                 benthos=="lobophora"|benthos=="lob_r" ~ "Lobophora",
                                 benthos=="halimeda" ~ "Halimeda",
                                 benthos=="turf"|benthos=="turf_rubbl" ~ "Turf"))%>%as.data.frame()
  
  ### ==> collapse to major benthic group
  major_ben<-percent_l%>%
    group_by(name, year, site,site_state, transect, ben_group) %>%
    summarise(cover = sum(cover))%>%
    as.data.frame()
  
#### ===> isolate major coral genera <=== ####
  ### ==> filter to hard coral long format
  hc<-percent_l %>%
    subset(ben_group == "Hard coral") %>%
    mutate(genera = case_when(benthos=="acr_bra"|benthos=="acr_corym"|benthos=="acr_dig"|benthos=="acr_sub"|benthos=="acr_tab" ~ "Acropora",
                              benthos=="ast_plat"|benthos=="astreo_enc" ~ "Astreopora",
                              benthos=="b_monti"|benthos=="monmas"|benthos=="mon_pla"|benthos=="monti_encr"|benthos=="monti_fol"|benthos=="mont_sub_m" ~ "Montipora",
                              benthos=="dipmas" ~ "Dipsastrea",
                              benthos=="fav_mass" ~ "Favites",
                              benthos=="fungia" ~ "Fungia",
                              benthos=="lpsenc" ~ "Leptoseris",
                              benthos=="pavenc"|benthos=="pavfol"|benthos=="pavmas"|benthos=="pav_sub" ~ "Pavona",
                              benthos=="plamas" ~ "Platygyra",
                              benthos=="pocill" ~ "Pocillopora",
                              benthos=="por_mass"|benthos=="por_sm"|benthos=="potenc" ~ "Porites",
                              benthos=="psammo" ~ "Psammocora",
                              benthos=="styloph" ~ "Stylophora"))
  ### ==> collapse to genera
  genera<-hc%>%
    group_by(name, year, site, site_state,transect, genera) %>%
    summarise(cover = sum(cover))%>%
    as.data.frame()

# ========================================================================================
#
#                            ####  ~~~~  dbRDA  ~~~ ####
#
# ========================================================================================
set.seed(1984)
#### ===> Set up data <=== ####
  ### ==> meta data
  com_meta<-select(percent_w, name, site, transect, year, site_state)

  ### ==> isolate living community - dropping sand only
  com_raw<-select(percent_w, acr_bra:sarco, bran_calc:turf_rubbl)
    ## square root transform data
    com_t<-sqrt(com_raw)
    ## convert to bray-curtis dissimilarity matrix
    com_d<-vegdist(com_t, method = "bray")
  ### ==> check unconstrained pattern using nmds
  mds <- metaMDS(com_d, k = 2, trymax = 100)
  mds$stress   # < 0.2 usable, < 0.1 good  -> 0.1683329
  plot(mds)
  ordiellipse(mds, com_meta$site_state, kind = "sd", label = TRUE)
  
#### ===> Assess dispersion of the communities <=== ####
  ### ==> Zero adjusted 
  bd_all<-betadisper(com_d, interaction(com_meta$site, com_meta$year), add = "lingoes")
    # ~ Note: adding "lingoes" because of negative values
    permutest(bd_all, permutations = 9999)  
    # group means — the actual effect size
    round(tapply(bd_all$distances, bd_all$group, mean), 3)
    boxplot(bd_all, las = 2, xlab = "")
    # which pairs
    permutest(bd_all, pairwise = TRUE, permutations = 9999)
    # which factor
    summary(aov(bd_all$distances ~ com_meta$site * com_meta$year))
    # check balance
    table(com_meta$site, com_meta$year)
    
  ### ==> Raw: for reporting
  bd_raw <- betadisper(com_d, interaction(com_meta$site, com_meta$year))
  round(tapply(bd_raw$distances, bd_raw$group, mean), 3)

#### ===> Assess the transect effect <=== ####
# ~ Note: this is to justify running at quadrat level and not transect level
adonis2(com_d ~ site_state/transect, data = com_meta, by = "terms", permutations = 9999)
  # site_state explains 46.2% of variation
  # transect within each site_state explains 2.7% of variation -> effect is minimal even if signficant

#### ===> compute dbRDA <=== #### 
  ### ==> dbRDA  normal model
  com_db<-dbrda(com_d ~ site * year, data = com_meta)
  com_db
    ## adjust R2
    RsquareAdj(com_db)
    
  ### ==> assess imaginary fraction
  ev <- eigenvals(com_db, model = "unconstrained")
  neg <- ev[ev < 0]
  sum(neg)                              # total negative inertia
  abs(sum(neg)) / sum(abs(ev))          # as a proportion of unconstrained
  abs(sum(neg)) / com_db$tot.chi        # as a proportion of total
  length(neg)  

  range(rowSums(com_raw))
  sum(as.matrix(com_d) > 0.95) / length(com_d)   # proportion of saturated pairs
  
    ## dbRDA sensitivity model check because of the imaginary fraction being high
    com_db_sq <- dbrda(sqrt(com_d) ~ site * year, data = com_meta)
    com_db_sq
    RsquareAdj(com_db_sq)
    anova(com_db_sq, by = "margin", permutations = 9999)
    
  ### ==> run significance tests
  anova(com_db, by = "margin", permutations = 9999)
  anova(com_db, by = "terms",  permutations = 9999)   # site, year, site:year
  anova(com_db, by = "axis",   permutations = 9999)
  
#### ===> assessing the dbRDA for spatial and temporal effect <=== ####  
  #### ===> Pair-wise comparisons between sites
    ## check the marginal factor centroids
    scores(com_db, display = "cn")  
    ## check site x year cells
    sc <- scores(com_db, display = "sites", scaling = 2)
    cell_cen <- aggregate(sc[, 1:2], by = list(cell = com_meta$cell), FUN = mean)
    cell_cen
      # pre→post vector per site
      for (s in levels(com_meta$site)) {
        pre  <- cell_cen[cell_cen$cell == paste0(s, ".Pre"),  2:3]
        post <- cell_cen[cell_cen$cell == paste0(s, ".Post"), 2:3]
        cat(s, ": d1 =", round(post[[1]] - pre[[1]], 3),
            " d2 =", round(post[[2]] - pre[[2]], 3), "\n")
      }
    
    ## homogenisation check
    D <- as.matrix(dist(cell_cen[, 2:3]))
    dimnames(D) <- list(cell_cen$cell, cell_cen$cell)
    round(D, 3)
      # did between-site distances shrink?
      pre_pairs  <- c(D["Offshore.Pre","Midshore.Pre"], D["Offshore.Pre","Inshore.Pre"],
                      D["Midshore.Pre","Inshore.Pre"])
      post_pairs <- c(D["Offshore.Post","Midshore.Post"], D["Offshore.Post","Inshore.Post"],
                      D["Midshore.Post","Inshore.Post"])
      data.frame(pair = c("Off-Mid","Off-In","Mid-In"), pre = pre_pairs, post = post_pairs,
                 change = post_pairs - pre_pairs)
      # verify results on the disimilarity space
      bd_cell <- betadisper(com_d, com_meta$cell)
      Dc <- as.matrix(dist(bd_cell$centroids))
      dimnames(Dc) <- list(levels(com_meta$cell), levels(com_meta$cell))
      
      data.frame(
        pair = c("Off-Mid", "Off-In", "Mid-In"),
        pre  = c(Dc["Offshore.Pre","Midshore.Pre"], Dc["Offshore.Pre","Inshore.Pre"],
                 Dc["Midshore.Pre","Inshore.Pre"]),
        post = c(Dc["Offshore.Post","Midshore.Post"], Dc["Offshore.Post","Inshore.Post"],
                 Dc["Midshore.Post","Inshore.Post"])
      )
      # check whether midshore post became more like inshore pre
      Dc["Midshore.Post", "Inshore.Pre"]    # vs Dc["Midshore.Pre","Inshore.Pre"] = 0.512
   
#### ==> pair-wise comparisons within sites using SIMPER and cohens-d <=== ####
  ## calculate shift 
  data.frame(
    site  = c("Offshore", "Midshore", "Inshore"),
    shift = c(Dc["Offshore.Pre","Offshore.Post"],
              Dc["Midshore.Pre","Midshore.Post"],
              Dc["Inshore.Pre","Inshore.Post"]))
  ## test shift
  res_year <- do.call(rbind, lapply(levels(com_meta$site), function(s) {
    i  <- com_meta$site == s
    dm <- as.dist(as.matrix(com_d)[i, i])
    md <- droplevels(com_meta[i, ])
    a  <- adonis2(dm ~ year, data = md, permutations = 9999)
    data.frame(site = s, R2 = a$R2[1], F = a$F[1], p = a$`Pr(>F)`[1])
  }))
  res_year$p_holm <- p.adjust(res_year$p, method = "holm")
  res_year$shift  <- c(Dc["Offshore.Pre","Offshore.Post"],
                       Dc["Midshore.Pre","Midshore.Post"],
                       Dc["Inshore.Pre","Inshore.Post"])
  res_year
  
### ===> Difference in composition
# per-site SIMPER
simp <- lapply(levels(com_meta$site), function(s) {
  i <- com_meta$site == s
  x <- summary(simper(com_t[i, ], com_meta$year[i], permutations = 999))[[1]]
  x$taxon <- rownames(x)
  
  # raw mean cover, keyed by taxon, from the metadata — not from simper
  pre  <- colMeans(com_raw[i & com_meta$year == "Pre",  , drop = FALSE])
  post <- colMeans(com_raw[i & com_meta$year == "Post", , drop = FALSE])
  
  x$pre_pct    <- pre[x$taxon]
  x$post_pct   <- post[x$taxon]
  x$change_pct <- x$post_pct - x$pre_pct
  x$direction  <- ifelse(x$change_pct > 0, "increase", "decrease")
  
  x <- x[x$ratio > 1, c("taxon", "average", "ratio", "p",
                        "pre_pct", "post_pct", "change_pct", "direction")]
  rownames(x) <- NULL
  x
})
  names(simp) <- levels(com_meta$site)

  lapply(simp, function(x) {
    x[, c("average","ratio","p","pre_pct","post_pct","change_pct")] <-
      round(x[, c("average","ratio","p","pre_pct","post_pct","change_pct")], 3)
    head(x, 10)
  })
  ## effect sizes of the taxa
  sel <- function(s, y) com_raw[com_meta$site == s & com_meta$year == y, , drop = FALSE]
  
  cohen_d <- function(s) {
    a <- sel(s, "Pre"); b <- sel(s, "Post")
    sp <- sqrt(((nrow(a)-1)*apply(a,2,var) + (nrow(b)-1)*apply(b,2,var)) /
                 (nrow(a)+nrow(b)-2))
    (colMeans(b) - colMeans(a)) / sp
  }
  round(sapply(c("Offshore","Midshore","Inshore"), cohen_d), 2)

#### ===> plot global dbRDA <=== ####
  ### ==> extract axis scores: dbRDA 1 and dbRDA2
  sit_g <- as.data.frame(scores(com_db, display = "wa", choices = 1:2, scaling = 2))
    # rename the axes
    names(sit_g) <- c("ax1","ax2")
    ## join site and year meta data to axes data
    sit_g$site <- com_meta$site 
    sit_g$year <- com_meta$year
    sit_g$site_state <- com_meta$site_state
    
  ### ==> extract centroids
  cen_g <- aggregate(sit_g[,1:2], by = list(site = sit_g$site, year = sit_g$year), mean)
    ## convert to a wide format 
    arr_g <- merge(subset(cen_g, year == "Pre"),  subset(cen_g, year == "Post"),
                   by = "site", suffixes = c("_pre","_post"))
    
  ### ==> calculate axis variation contribution values
  ev_con <- eigenvals(com_db, model = "constrained")
    ## % of total inertia — matches summary(), use this for axis labels
    pct_tot <- round(100 * ev_con[1:2] / com_db$tot.chi, 1)
    ## % of constrained inertia — how the constrained signal splits across axes
    pct_con <- round(100 * ev_con[1:2] / sum(ev_con), 1)
    
  rbind(pct_tot, pct_con)
  
  ### ==> plot 
  plot.dbRDA_global <- ggplot(sit_g, aes(ax1, ax2)) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_point(aes(colour = site, shape = year), size = 1.5, alpha = 0.40) +
    stat_ellipse(aes(colour = site, linetype = year), type = "t",
                 level = 0.95, linewidth = 0.6) +
    geom_segment(data = arr_g,
                 aes(x = ax1_pre, y = ax2_pre, xend = ax1_post, yend = ax2_post,
                     colour = site),
                 arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
                 linewidth = 1.1, inherit.aes = FALSE, show.legend = FALSE) +
    scale_colour_manual(values = pal_cols, name = "Site") +
    scale_shape_manual(values = c(Pre = 1, Post = 16), name = "Period") +
    scale_linetype_manual(values = c(Pre = "dashed", Post = "solid"), name = "Period") +
    coord_equal() +
    labs(x = paste0("dbRDA1 (", pct_tot[1], "%)"),
         y = paste0("dbRDA2 (", pct_tot[2], "%)")) + #subtitle = "A   All sites: site \u00d7 period structure"
    theme_classic(base_size = 11) +
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          legend.key.height = unit(0.8, "lines"))
  plot.dbRDA_global
  
#### ===> plot site specific dbRDA <=== ####
# One two-level factor per model -> exactly 1 constrained axis -> is the pre/post gradient
# Plotted against the first residual (MDS1) axis -> allows for direct interpretation
  ### ==> Set up extraction for site level
  d_full <- as.matrix(com_d)
  mult   <- 0.5      # species arrow scaling -- tune after a first look
  p_cut  <- 0.05     # SIMPER p threshold for which taxa to draw
  
  ### ==> extract model details per site
  site_out <- lapply(levels(com_meta$site), function(s) {
    
    i  <- com_meta$site == s
    md <- droplevels(com_meta[i, ])
    m  <- dbrda(as.dist(d_full[i, i]) ~ year, data = md)
    sppscores(m) <- com_t[i, ]
    
    # Rank-1 constrained model: take the constrained axis and the first residual
    # axis directly, which avoids scores(choices = 1:2) failing on rank 1.
    sit <- data.frame(ax1 = m$CCA$wa[, 1], ax2 = m$CA$u[, 1],
                      year = md$year, site = s)
    
    sp <- as.data.frame(scores(m, display = "species", scaling = 2))
    sp <- sp[, 1:2, drop = FALSE]
    names(sp) <- c("ax1", "ax2")
    sp$taxon <- rownames(sp)
    sp$site  <- s
    
    # The sign of a constrained axis is arbitrary -- force Post to the right so
    # all three panels read the same way.
    if (mean(sit$ax1[sit$year == "Post"]) < 0) {
      sit$ax1 <- -sit$ax1
      sp$ax1  <- -sp$ax1
    }
    
    # Draw only taxa that discriminate pre/post beyond chance
    keep <- simp[[s]]$taxon[simp[[s]]$p < p_cut]
    keep <- keep[!is.na(keep)]
    sp   <- sp[sp$taxon %in% keep, ]
    
    pv  <- anova(m, permutations = 9999)
    lab <- sprintf("%s   (R\u00b2 = %.3f, p = %.4f)",
                   s, RsquareAdj(m)$r.squared, pv$`Pr(>F)`[1])
    
    list(sit = sit, sp = sp, lab = setNames(lab, s))
  })
    
  ### ==> create data frame for plotting
  sit_b <- bind_rows(lapply(site_out, `[[`, "sit"))
  sp_b  <- bind_rows(lapply(site_out, `[[`, "sp"))
  labs  <- unlist(lapply(site_out, `[[`, "lab"))
    # form factors
    sit_b$site <- factor(sit_b$site, levels = levels(com_meta$site))
    sp_b$site  <- factor(sp_b$site,  levels = levels(com_meta$site))

  ### ==> plot
  plot.dbRDA_site <- ggplot(sit_b, aes(ax1, ax2)) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_point(aes(colour = year, shape = year), size = 1.5, alpha = 0.60) +
    stat_ellipse(aes(colour = year, linetype = year), type = "t", level = 0.95, linewidth = 0.7) +
    geom_segment(data = sp_b,
                 aes(x = 0, y = 0, xend = ax1 * mult, yend = ax2 * mult),
                 arrow = arrow(length = unit(0.15, "cm")),
                 colour = "grey25", linewidth = 0.4, inherit.aes = FALSE) +
    geom_text_repel(data = sp_b,
                    aes(x = ax1 * mult, y = ax2 * mult, label = taxon),
                    size = 2.9, fontface = "italic", colour = "grey15",
                    segment.colour = "grey65", min.segment.length = 0.2,
                    box.padding = 0.35, max.overlaps = Inf, inherit.aes = FALSE) +
    facet_wrap(~ site, nrow = 1, scales = "free",
               labeller = labeller(site = labs)) +
    scale_color_viridis_d(option="magma", begin = 0.2, end=0.6)+
    scale_shape_manual(values = c(Pre = 1, Post = 16), name = NULL) +
    scale_linetype_manual(values = c(Pre = "dashed", Post = "solid"), name = "Period") +
    labs(x = "dbRDA1 (constrained: pre \u2192 post)",
         y = "MDS1 (residual)") + #subtitle = "B   Within each site: pre vs post, with discriminating taxa"
    theme_classic(base_size = 11) +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "grey96", colour = NA),
          strip.text = element_text(size = 9),
          legend.position = "bottom")
    plot.dbRDA_site 
    
#### ===> stitch figures <=== #### 
plot.dbRDA <- plot.dbRDA_global / plot.dbRDA_site + 
      plot_layout(heights = c(1.15, 1))+
      plot_annotation(tag_levels = "A")
plot.dbRDA
  ## Save plot
  #ggsave("Feedback_Outputs/dbRDA.pdf", 
  #       plot=plot.dbRDA, 
  #       width = 5, height = 6, units = "in", scale = 1.5, dpi = 600) 
  
# ========================================================================================
#
#                       ####  ~~~~  PLOT: changes to percent cover  ~~~ ####
#
# ========================================================================================
#### ===> box plot 1: benthic group -> hard coral <=== ####
plot.box_hc<-major_ben%>%
  filter(ben_group=="Hard coral")%>%
  ggplot(aes(x=site, y=cover, fill=year, color=year))+
  geom_point(shape=21, alpha=0.5, size=1,
             position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.3))+
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6,position = position_dodge(width=0.7))+
  labs(y="Hard coral (%)", x="", fill="Period", color="Period")+
  scale_fill_viridis_d(option="magma", begin = 0.2, end= 0.6)+
  scale_color_viridis_d(option="magma", begin = 0.2, end=0.6)+
  theme_classic(base_size = 11)+
  theme(legend.position = "none")
plot.box_hc

#### ===> box plot 2: benthic group -> CCA <=== ####
plot.box_CCA<-major_ben%>%
  filter(ben_group=="CCA")%>%
  ggplot(aes(x=site, y=cover, fill=year, color=year))+
  geom_point(shape=21, alpha=0.5, size=1,
             position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.3))+
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6,position = position_dodge(width=0.7))+
  labs(y="CCA (%)", x="", fill="Period", color="Period")+
  scale_fill_viridis_d(option="magma", begin = 0.2, end= 0.6)+
  scale_color_viridis_d(option="magma", begin = 0.2, end=0.6)+
  theme_classic(base_size = 11)+
  theme(legend.position = "none")
plot.box_CCA

#### ===> box plot 3: benthic group -> Turf <=== ####
plot.box_turf<-major_ben%>%
  filter(ben_group=="Turf")%>%
  ggplot(aes(x=site, y=cover, fill=year, color=year))+
  geom_point(shape=21, alpha=0.5, size=1,
             position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.3))+
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6,position = position_dodge(width=0.7))+
  labs(y="Turf (%)", x="", fill="Period", color="Period")+
  scale_fill_viridis_d(option="magma", begin = 0.2, end= 0.6)+
  scale_color_viridis_d(option="magma", begin = 0.2, end=0.6)+
  theme_classic(base_size = 11)+
  theme(legend.position = "none")
plot.box_turf

#### ===> box plot 4: algae -> Halimeda <=== ####
plot.box_hal<-major_ben%>%
  filter(ben_group=="Halimeda")%>%
  ggplot(aes(x=site, y=cover, fill=year, color=year))+
  geom_point(shape=21, alpha=0.5, size=1,
             position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.3))+
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6,position = position_dodge(width=0.7))+
  labs(y="Halimeda (%)", x="", fill="Period", color="Period")+
  scale_fill_viridis_d(option="magma", begin = 0.2, end= 0.6)+
  scale_color_viridis_d(option="magma", begin = 0.2, end=0.6)+
  theme_classic(base_size = 11)+
  theme(legend.position = "none")
plot.box_hal

#### ===> box plot 5: algae -> Lobophora <=== ####
plot.box_lob<-major_ben%>%
  filter(ben_group=="Lobophora")%>%
  ggplot(aes(x=site, y=cover, fill=year, color=year))+
  geom_point(shape=21, alpha=0.5, size=1,
             position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.3))+
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6,position = position_dodge(width=0.7))+
  labs(y="Lobophora (%)", x="", fill="Period", color="Period")+
  scale_fill_viridis_d(option="magma", begin = 0.2, end= 0.6)+
  scale_color_viridis_d(option="magma", begin = 0.2, end=0.6)+
  theme_classic(base_size = 11)+
  theme(legend.position = "none")
plot.box_lob

#### ===> box plot 6: algae -> Peysonellia <=== ####
plot.box_pey<-percent_w%>%
  ggplot(aes(x=site, y=peysson, fill=year, color=year))+
  geom_point(shape=21, alpha=0.5, size=1,
             position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.3))+
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6,position = position_dodge(width=0.7))+
  labs(y="Peyssonnelia (%)", x="", fill="Period", color="Period")+
  scale_fill_viridis_d(option="magma", begin = 0.2, end= 0.6)+
  scale_color_viridis_d(option="magma", begin = 0.2, end=0.6)+
  theme_classic(base_size = 11)+
  theme(legend.position = "bottom")
plot.box_pey

#### ===> box plot 7: algae -> Caluerpa <=== ####
plot.box_cau<-percent_w%>%
  ggplot(aes(x=site, y=caulerpa, fill=year, color=year))+
  geom_point(shape=21, alpha=0.5, size=1,
             position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.3))+
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6,position = position_dodge(width=0.7))+
  labs(y="Caulerpa (%)", x="", fill="Period", color="Period")+
  scale_fill_viridis_d(option="magma", begin = 0.2, end= 0.6)+
  scale_color_viridis_d(option="magma", begin = 0.2, end=0.6)+
  theme_classic(base_size = 11)+
  theme(legend.position = "none")
plot.box_cau

#### ===> box plot 8: coral genera -> Montipora <=== ####
plot.box_mon<-genera%>%
  filter(genera=="Montipora")%>%
  ggplot(aes(x=site, y=cover, fill=year, color=year))+
  geom_point(shape=21, alpha=0.5, size=1,
             position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.3))+
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6,position = position_dodge(width=0.7))+
  labs(y="Montipora (%)", x="", fill="Period", color="Period")+
  scale_fill_viridis_d(option="magma", begin = 0.2, end= 0.6)+
  scale_color_viridis_d(option="magma", begin = 0.2, end=0.6)+
  theme_classic(base_size = 11)+
  theme(legend.position = "none")
plot.box_mon

#### ===> box plot 9: coral genera -> Acropora <=== ####
plot.box_acr<-genera%>%
  filter(genera=="Acropora")%>%
  ggplot(aes(x=site, y=cover, fill=year, color=year))+
  geom_point(shape=21, alpha=0.5, size=1,
             position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.3))+
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6,position = position_dodge(width=0.7))+
  labs(y="Acropora (%)", x="", fill="Period", color="Period")+
  scale_fill_viridis_d(option="magma", begin = 0.2, end= 0.6)+
  scale_color_viridis_d(option="magma", begin = 0.2, end=0.6)+
  theme_classic(base_size = 11)+
  theme(legend.position = "none")
plot.box_acr

#### ===> stitch plots together <=== ####
  ## layout
  layout_pc<- "
              ABC
              DEF
              GHI
             "
  
  ## combine plots
  plot.pc<-(
           plot.box_hc+plot.box_mon+plot.box_acr+
           plot.box_turf+plot.box_CCA+plot.box_hal+
           plot.box_lob+plot.box_pey+plot.box_cau
           )+
    plot_layout(design = layout_pc)+
    plot_annotation(tag_levels = "A")
  plot.pc
  
  ## Save plot
  ggsave("Feedback_Outputs/Percent_cover.pdf", 
         plot=plot.pc, 
         width = 5, height = 6, units = "in", scale = 1.5, dpi = 600) 

# ========================================================================================
#
#             ####  ~~~~  statistical testing changes to cover  ~~~ ####
#
# ========================================================================================
set.seed(1984)
#### ===> hard coral <== ####
  ### ==> set up
    ## meta data
    meta.hc<-major_ben%>%filter(ben_group=="Hard coral")%>%select(year, site, site_state)
    ## cover data
    col.hc<-major_ben%>%filter(ben_group=="Hard coral")%>%select(cover)
      # convert to distance matrix
      dis.hc<-vegdist(col.hc, method = "euclidean")
      
  ### ==> run PERMANOVA
  perm.hc<-adonis2(dis.hc~site*year, data = meta.hc, permutations = 9999, by = "margin")
  perm.hc
    ## pair-wise testing
    perm.hc.period<-pairwise.adonis2(dis.hc~site_state, data = meta.hc, nperm = 9999)
    perm.hc.period
      
#### ===> CCA <=== ####
  ### ==> set up
    ## meta data
    meta.cca<-major_ben%>%filter(ben_group=="CCA")%>%select(year, site, site_state)
    ## cover data
    col.cca<-major_ben%>%filter(ben_group=="CCA")%>%select(cover)
      # convert to distance matrix
      dis.cca<-vegdist(col.cca, method = "euclidean")
    
  ### ==> run PERMANOVA
  perm.cca<-adonis2(dis.cca~site*year, data = meta.cca, permutations = 9999, by = "margin")
  perm.cca
    ## pair-wise testing
    perm.cca.period<-pairwise.adonis2(dis.cca~site_state, data = meta.cca, nperm = 9999)
    perm.cca.period
    
#### ===> Turf <=== ####
  ### ==> set up
    ## meta data
    meta.turf<-major_ben%>%filter(ben_group=="Turf")%>%select(year, site, site_state)
    ## cover data
    col.turf<-major_ben%>%filter(ben_group=="Turf")%>%select(cover)
      # convert to distance matrix
      dis.turf<-vegdist(col.turf, method = "euclidean")
    
  ### ==> run PERMANOVA
  perm.turf<-adonis2(dis.turf~site*year, data = meta.turf, permutations = 9999, by = "margin")
  perm.turf
    ## pair-wise testing
    perm.turf.period<-pairwise.adonis2(dis.turf~site_state, data = meta.turf, nperm = 9999)
    perm.turf.period
    
#### ===> Halimeda <=== ####
  ### ==> set up
    ## meta data
    meta.hal<-major_ben%>%filter(ben_group=="Halimeda")%>%select(year, site, site_state)
    ## cover data
    col.hal<-major_ben%>%filter(ben_group=="Halimeda")%>%select(cover)
      # convert to distance matrix
      dis.hal<-vegdist(col.hal, method = "euclidean")
    
  ### ==> run PERMANOVA
  perm.hal<-adonis2(dis.hal~site*year, data = meta.hal, permutations = 9999, by = "margin")
  perm.hal
    
#### ===> Lobophora <=== ####
  ### ==> set up
    ## meta data
    meta.lob<-major_ben%>%filter(ben_group=="Lobophora")%>%select(year, site, site_state)
    ## cover data
    col.lob<-major_ben%>%filter(ben_group=="Lobophora")%>%select(cover)
      # convert to distance matrix
      dis.lob<-vegdist(col.lob, method = "euclidean")
    
  ### ==> run PERMANOVA
  perm.lob<-adonis2(dis.lob~site*year, data = meta.lob, permutations = 9999, by = "margin")
  perm.lob
    ## pair-wise testing
    perm.lob.period<-pairwise.adonis2(dis.lob~site_state, data = meta.lob, nperm = 9999)
    perm.lob.period
    
#### ===> Peysonellia <=== ####
  ### ==> set up
    ## meta data
    meta.pey<-percent_w%>%select(year, site, site_state)
    ## cover data
    col.pey<-percent_w%>%select(peysson)
      # convert to distance matrix
      dis.pey<-vegdist(col.pey, method = "euclidean")
    
  ### ==> run PERMANOVA
  perm.pey<-adonis2(dis.pey~site*year, data = meta.pey, permutations = 9999, by = "margin")
  perm.pey
    ## pair-wise testing
    perm.pey.period<-pairwise.adonis2(dis.pey~site_state, data = meta.pey, nperm = 9999)
    perm.pey.period
    
#### ===> Caulerpa <=== ####
  ### ==> set up
    ## meta data
    meta.cau<-percent_w%>%select(year, site, site_state)
    ## cover data
    col.cau<-percent_w%>%select(caulerpa)
      # convert to distance matrix
      dis.cau<-vegdist(col.cau, method = "euclidean")
    
  ### ==> run PERMANOVA
  perm.cau<-adonis2(dis.cau~site*year, data = meta.cau, permutations = 9999, by = "margin")
  perm.cau
    ## pair-wise testing
    perm.cau.period<-pairwise.adonis2(dis.cau~site_state, data = meta.cau, nperm = 9999)
    perm.cau.period
    
#### ===> Montipora <=== ####
  ### ==> set up
    ## meta data
    meta.mon<-genera%>%filter(genera=="Montipora")%>%select(year, site, site_state)
    ## cover data
    col.mon<-genera%>%filter(genera=="Montipora")%>%select(cover)
      # convert to distance matrix
      dis.mon<-vegdist(col.mon, method = "euclidean")
    
  ### ==> run PERMANOVA
  perm.mon<-adonis2(dis.mon~site*year, data = meta.mon, permutations = 9999, by = "margin")
  perm.mon
    ## pair-wise testing
    perm.mon.period<-pairwise.adonis2(dis.mon~site_state, data = meta.mon, nperm = 9999)
    perm.mon.period
    
#### ===> Acropora <=== ####
  ### ==> set up
    ## meta data
    meta.acr<-genera%>%filter(genera=="Acropora")%>%select(year, site, site_state)
    ## cover data
    col.acr<-genera%>%filter(genera=="Acropora")%>%select(cover)
      # convert to distance matrix
      dis.acr<-vegdist(col.acr, method = "euclidean")
    
  ### ==> run PERMANOVA
  perm.acr<-adonis2(dis.acr~site*year, data = meta.acr, permutations = 9999, by = "margin")
  perm.acr
    ## pair-wise testing
    perm.acr.period<-pairwise.adonis2(dis.acr~site_state, data = meta.acr, nperm = 9999)
    perm.acr.period


  