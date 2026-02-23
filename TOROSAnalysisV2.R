TOROSAnalysis <- function(night_folders = NULL,
                          preload_data = F ,
                          preloaded_data = NULL ,
                          write_lcs = F ,
                          run_Var_Analysis = F ,
                          run_RMS_Analysis = F ,
                          lccsv_folder = "" ,
                          plot_lcs = F ,
                          lc_folder = "" ,
                          field_name = "Observed Field" ,
                          target = 0 ,
                          fit_data = F ,
                          fC_guess = NULL ,
                          plot_LSP = F ,
                          plot_plc = F ,
                          plot_lc = F ,
                          plc_use_fitted_P = F ,
                          plc_use_LSP_P = F ,
                          plc_LSP_per_num = 1 ,
                          plc_period = 0.0 ,
                          run_LSST_test = F ,
                          obs_plot_min = 0 ,
                          obs_plot_max = 10 ,
                          LSST_target_name = "" ,
                          LSST_N = 0 ,
                          LSST_sigma = 0 ,
                          LSST_M = 0
                          )
{

  ###Subfunctions to do specific tasks
  print("Loading in subfunctions...")

  loadInNights <- function(night_folders)
  {

    image_photometry_list <- list()
    caught_ids <- list()
    #We are loading in the aperture photometry data of an entire image

    for (folder in night_folders)
    {

      night_file_list <- list.files(path = folder, full.names = T) #Load in file names into a list

      for (file in night_file_list)
      {

        df <- read.csv(file = file) #Read night.csv into a dataframe
        image_photometry_list <- append(image_photometry_list , list(df))

        caught_ids <- append(caught_ids , list(df$alligned_id)) #Generate list of lists of IDs found in fields

      }
    }

    common_ids <- Reduce(f = intersect , caught_ids) #This will generate of list of IDs present among ALL fields

    return(list(CommonIDs = common_ids , Photometry = image_photometry_list))
  }
  print("loadInNights loaded!")

  makeStarStat <- function(data_list , id_list)
  {

    #We now will create individual csvs for light curves of stars
    starlc <- list()

    for (id in id_list)
    {

      temp_star_frame <- data.frame()

      for (df in data_list)
      {

        #Pull data from field with starID = id
        found_star_data <- df[df$alligned_id == id, ]

        #Create a time object
        found_star_data$TimeStamp <- as.POSIXct(
          paste(found_star_data$Date , found_star_data$Time) ,
          format = "%Y-%m-%d %H:%M:%OS", tz = "UTC"
        )


        temp_star_frame <- rbind(temp_star_frame , data.frame(
          StarID = found_star_data$alligned_id,
          Date = found_star_data$Date,
          Time = found_star_data$Time,
          TimeStamp = found_star_data$TimeStamp,
          Flux = found_star_data$f_n,
          Magnitude = (25 - 2.5*log10(found_star_data$f_n)),
          X_Position = found_star_data$x_n,
          Y_Position = found_star_data$y_n,
          RA = found_star_data$RA,
          DEC = found_star_data$DEC
        ))
      }

      temp_star_frame <- temp_star_frame[order(temp_star_frame$TimeStamp) , ] #Sorts by date
      starlc[[as.character(id)]] <- temp_star_frame #Put star lightcurve into a dict

    }

    #We now will account for any offsets per night
    night_medians <- list()
    unique_dates <- unique(starlc[[1]]$Date)#Obtains unique datees

    #Find median magnitude from ALL lcs for each night
    for (date in unique_dates)
    {

      mags_one_night <- lapply(starlc , function(df) { df <- df[df$Date == date , ]}) #Creates list where all LCS are from same night
      median_mag_one_night <- median(sapply(mags_one_night , function(df) { median(df$Magnitude) })) #Gets median mag for this night
      night_medians <- append(night_medians , median_mag_one_night) #Places that night's median into a list

    }

    #Compute median magnitude offsets with night 1 as the baseline
    offset_list <- list()
    night_num <- 1
    for (med in night_medians)
    {

      offset <- night_medians[[1]] - med
      print(paste0("Offset for night " , night_num , " is " , offset, "."))
      offset_list <- append(offset_list , offset)

      #The we subtract off that offset from all magnitudes with that date!
      starlc <- lapply(starlc , function(df)
      {

        df[df$Date == unique_dates[night_num] , ]$Magnitude <- df[df$Date == unique_dates[night_num] , ]$Magnitude + rep(offset , length(df[df$Date == unique_dates[night_num] , ]$Magnitude))
        return(df)
        #The lapply essentially ensures that the rows with the correct date only have the offset removed

      })

      night_num <- night_num + 1

    }





    return(starlc)
  }
  print("makeStarStat loaded!")

  systematicRemoval <- function(lcs)
  {
    new_lcs <- lcs

    library(matrixStats)
    library(magicaxis)
    #We will choose a star, find all stars within a 1500 pixel radius of similar brightness, and build up a systematic to remove from the data.
    for (tar_star in lcs)
    {

      sys_removed <- T

      #We find this star's median magnitude and position
      target_med_mag <- median(tar_star$Magnitude)
      target_x <- mean(tar_star$X_Position)
      target_y <- mean(tar_star$Y_Position)

      #We next will loop through every other star and build up a list of nearby stars
      similar_stars <- list()
      sim_idx <- 1

      for (star in lcs)
      {

        star_x <- mean(star$X_Position)
        star_y <- mean(star$Y_Position)
        star_med_mag <- median(star$Magnitude)

        delta_pos <- sqrt((target_x - star_x)^2 + (target_y - star_y)^2) #Gets distance to target
        delta_mag <- abs(target_med_mag - star_med_mag) #Gets mag diff
        sc_rms <- sd(magclip(star$Magnitude , sigma = 5)$x) #Gets RMS of star mag, only want stable stars

        if (delta_pos > 0 && delta_pos < 1000 && delta_mag < 1.2 && sc_rms < .05)
        {

          #If a star is close enough and similar brightness, add it to the list
          similar_stars[[sim_idx]] <- star$StarID[1]
          sim_idx <- sim_idx + 1

        }
      }

      if (length(similar_stars) < 2)
      {

        sys_removed <- F

      }

      #We will build a systematic by obtaining the medians of all lightcurves per timestep and removing that systematic from the LC
      unique_dates <- sort(unique(tar_star$Date))
      if (sys_removed == T)
      {

        systematics <- c() #Will be a vector of systematics, systematics for each date of obs will be appended to vector

        for (date in unique_dates)
        {

          lcs_from_date <- lapply(lcs , function(df) df[df$Date == date , ]) #Creates lc list of similar stars for one date
          similar_lcs <- lapply(similar_stars , function(id) lcs_from_date[[as.character(id)]]) #Obtain lcs of similar stars
          mag_mat <- sapply(similar_lcs , function(df) df$Magnitude) #Builds matrix where columns are stars and rows are timesteps
          norm_mag_mat <- sweep(mag_mat , 2 , colMedians(mag_mat) , FUN = '-') #Normalize columns by median of star's lc
          systematic <- rowMedians(norm_mag_mat , na.rm = T) #Systematic is the median normalized changes in magnitude among all lightcurves
          systematics <- c(systematics , systematic)

        }


      }


      if (sys_removed == F)
      {

        systematics <- rep(0 , nrow(tar_star))

      }

      #Now we put this data into new_lcs so we can return it with more info about the stars
      new_lcs[[as.character(tar_star$StarID[1])]]$AMag <- tar_star$Magnitude - systematics
      new_lcs[[as.character(tar_star$StarID[1])]]$SMag <- systematics
      new_lcs[[as.character(tar_star$StarID[1])]]$SysRemoved <- rep(sys_removed , nrow(tar_star))



    }


    return(new_lcs)

  }
  print("systematicRemoval loaded!")

  writeLcCsvs <- function(lcs , field_name , lccsv_folder)
  {

    #This function just writes the csvs to disk

    for (star in lcs)
    {

      write.csv(star , paste0(lccsv_folder , field_name , star$StarID[1] , ".csv"))

    }
  }
  print("writeLcCsvs loaded!")

  makeLightCurves <- function(lcs , lc_folder , field_name)
  {
    library(ggplot2)
    library(lubridate)
    library(dplyr)

    for (star in lcs)
    {

      #We need to create a useful time axis
      star <- star %>%
        mutate(
          full_datetime = ymd_hms(paste(Date , Time)),
          night_start = as.POSIXct(Date , format = "%Y-%m-%d" , tz = "UTC"),
          intra_night_time = as.numeric(difftime(full_datetime , night_start , units = "mins"))
        )

      unique_dates <- sort(unique(star$Date))
      date_to_index <- setNames(seq_along(unique_dates) - 1 , unique_dates)
      star <- star %>%
        mutate(
          night_index = date_to_index[Date]
        )

      max_night_duration <- star %>%
        group_by(Date) %>%
        summarize(max_time = max(intra_night_time)) %>%
        pull(max_time) %>%
        max()

      star <- star %>%
        mutate(

          x_time = night_index * (max_night_duration + 1.5) + intra_night_time

        )



      p <- ggplot(data = star , aes(x = x_time , color = Date)) +
        geom_point(aes(y = Magnitude , shape = "Raw Mag") , size = 2 , alpha = .5) +
        geom_point(aes(y = AMag , shape = "Adj Mag") , size = 2 , alpha = .5) +
        scale_y_reverse() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs( color = "Viewing Night" ,
              x = "Run Time for Individual Night [min]" ,
              y = "Instrumental Magnitude" ,
              title = paste0("Lightcurve for " , field_name , star$StarID[1]) ,
              subtitle = paste0("Systematic Removal: " , star$SysRemoved[1]) ,
              shape = "Legend" ) +
        scale_color_discrete() +
        scale_shape_manual(values = c(
          "Raw Mag" = 4,
          "Adj Mag" = 16
        )) +
        theme(
          plot.title = element_text(size = 25 , face = "bold") ,
          plot.subtitle = element_text(size = 20) ,
          axis.title.x = element_text(size = 20 , face = "bold") ,
          axis.title.y = element_text(size = 20 , face = 'bold') ,
          axis.text.x = element_text(size = 15) ,
          axis.text.y = element_text(size = 15) ,
          legend.title = element_text(size = 12 , face = "bold") ,
          legend.text = element_text(size = 12)
      )

      ggsave(
        filename = paste0(field_name , star$StarID[1] , "_lightcurve.pdf") ,
        plot = p ,
        path = lc_folder ,
        width = 8 , height = 6 , units = 'in'
      )

      # print(p)
      if (star$SysRemoved[1] == T)
      {

        pp <- ggplot(data = star , aes(x = x_time , color = Date)) +
          geom_line(aes( y = SMag ) , group = 1) +
          scale_y_reverse() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          labs( color = 'Viewing Night' ,
                x = 'Time' ,
                y = 'Systematic Magnitude' ,
                title = paste0("Systematic for " , field_name , star$StarID[1])) +
          scale_color_discrete() +
          theme(
            plot.title = element_text(size = 25 , face = "bold") ,
            #plot.subtitle = element_text(size = 20) ,
            axis.title.x = element_text(size = 20 , face = "bold") ,
            axis.title.y = element_text(size = 20 , face = 'bold') ,
            axis.text.x = element_text(size = 15) ,
            axis.text.y = element_text(size = 15) ,
            legend.title = element_text(size = 12 , face = "bold") ,
            legend.text = element_text(size = 12)
      )

        ggsave(
          filename = paste0(field_name , star$StarID[1] , "systematic_lightcurve.pdf") ,
          plot = pp ,
          path = lc_folder ,
          width = 8 , height = 6 , units = 'in'
        )

      }





















    }




  }
  print("makeLightCurves loaded!")

  fitCurveToData <- function(lcs , fC_target , field_name , guess)
  {

    library(ggplot2)

    #Lets get the running day number of our target
    tar_star <- lcs[[as.character(fC_target)]]
    tar_star$NumericTime <- as.numeric(
      difftime(tar_star$TimeStamp , min(tar_star$TimeStamp) , units = 'days')
    )

    #We want to fit model to our running day vs AMag measurements
    what_to_fit_to <- data.frame(Running_Day = tar_star$NumericTime , Mag = tar_star$AMag)

    #Do the fitting
    fitting_data <- nls(Mag ~ A_fit * sin(2 * pi * (Running_Day - H_fit) / P_fit) + V_fit ,
                        data = what_to_fit_to ,
                        start = list(A_fit = guess[1] , P_fit = guess[2] , H_fit = guess[3] , V_fit = guess[4]))

    #Save the params
    A <- coef(fitting_data)[[1]]
    P <- coef(fitting_data)[[2]]
    H <- coef(fitting_data)[[3]]
    V <- coef(fitting_data)[[4]]

    #Create sinusoidal function
    model <- function(t , A , P , H , V)
    {

      return (A * sin(2 * pi * (t - H) / P) + V)

    }

    #Plot the data against the fitted curve
    fitting_time <- seq(min(tar_star$NumericTime) , max(tar_star$NumericTime) , length.out = 1000)
    curve_mag <- model(t = fitting_time , A = A , P = P , H = H , V = V)

    p <- ggplot() +
      geom_point(aes( x = tar_star$NumericTime , y = tar_star$AMag , color = "Observations") , size = 1) +
      geom_line(aes( x = fitting_time , y = curve_mag , color = "Fitted Function")) +
      labs(title = paste0("Fitted Curve to " , field_name , tar_star$StarID) ,
           subtitle = paste0("Fitted Amplitude: " , A , " | Fitted Period: " , P , " days") ,
           x = "Running Time [days]" ,
           y = "Instrumental Magnitude"
           ) +
      scale_color_manual(
        name = "Legend" ,
        values = c(
          "Observations" = "black" ,
          "Fitted Function" = "red"
        )
      ) +
      scale_y_reverse() +
      theme(
        plot.title = element_text(size = 25 , face = "bold") ,
        plot.subtitle = element_text(size = 20) ,
        axis.title.x = element_text(size = 20 , face = "bold") ,
        axis.title.y = element_text(size = 20 , face = 'bold') ,
        axis.text.x = element_text(size = 15) ,
        axis.text.y = element_text(size = 15) ,
        legend.title = element_text(size = 12 , face = "bold") ,
        legend.text = element_text(size = 12)
      )


    print(p)

    fitting_list <- list(Fit_P = P , Fit_A = A , Fit_H = H , Fit_V = V)
    return(fitting_list)




  }
  print("fitCurveToData loaded!")

  lombScarg <- function(lcs , lsp_target , field_name , plc_lsp_num)
  {

    library(lomb)

    #Run a LSP on the data of a specific star
    tar_star <- lcs[[as.character(lsp_target)]]

    #Lets get the running day number of our target
    tar_star <- lcs[[as.character(lsp_target)]]
    tar_star$NumericTime <- as.numeric(
      difftime(tar_star$TimeStamp , min(tar_star$TimeStamp) , units = 'days')
    )

    #Create LSP object
    lomb <- lsp(x = tar_star$AMag , times = tar_star$NumericTime , type = "period" , ofac = 3 , plot = F)

    #Find the top ranked periods
    ord <- order(lomb$power , decreasing = T)
    top_periods <- data.frame(Period = lomb$scanned[ord] , Power = lomb$power[ord])

    #Clean up harmonics
    clean_periods <- top_periods

    for (i in 1:(nrow(clean_periods) - 1))
    {

      ratios <- clean_periods$Period[i] / clean_periods$Period[(i + 1):nrow(clean_periods)] #Obtain ratios of current period vs all other periods worse than it
      bad_ratios <- which(abs(ratios - round(ratios)) < .01) #Get indices of periods with nearly integer ratios

      #If there are bad ratios, kill them
      if (length(bad_ratios) != 0)
      {

        #Take these indices out
        clean_periods <- clean_periods[-(i + bad_ratios)]

      }
    }

    #Plot the periodogram with extra data
    p <- ggplot() +
      geom_line(aes( x = clean_periods$Period , y = clean_periods$Power , color = "LSP") , size = 1) +
      geom_hline(aes(yintercept = lomb$sig.level , color = "P < .01") , linetype = 'dotted' , size = 1.5) +
      geom_vline(aes(xintercept = clean_periods$Period[1] , color = '1st Period' ) , linetype = 'dashed' , size = 1.8) +
      geom_vline(aes(xintercept = clean_periods$Period[2] , color = '2nd Period' ) , linetype = 'dashed' , size = 1.8) +
      geom_vline(aes(xintercept = clean_periods$Period[3] , color = '3rd Period' ) , linetype = 'dashed' , size = 1.8) +
      labs(title = paste0("Lomb-Scargle Periodogram for " , field_name , " " , lsp_target) ,
           subtitle = paste0("1st Period: " , round(clean_periods$Period[1] , digits = 3) , " days | 1st Power: " , round(clean_periods$Power[1] , digits = 3) , " | 1st Period P-Value: " , lomb$p.value ,
                             "\n2nd Period: " , round(clean_periods$Period[2] , digits = 3) , " days | 2nd Power: " , round(clean_periods$Power[2] , digits = 3) ,
                             "\n3rd Period: " , round(clean_periods$Period[3] , digits = 3) , " days | 3rd Power: " , round(clean_periods$Power[3] , digits = 3)

                             ) ,
           x = "Period [days]" ,
           y = "Power"
           ) +
      scale_color_manual(
        name = "Legend" ,
        values = c(
          "LSP" = "black" ,
          "P < .01" = "red" ,
          "1st Period" = "purple" ,
          "2nd Period" = "blue" ,
          "3rd Period" = "cyan"
        )
      ) +
      theme(
        plot.title = element_text(size = 25 , face = "bold") ,
        plot.subtitle = element_text(size = 20) ,
        axis.title.x = element_text(size = 20 , face = "bold") ,
        axis.title.y = element_text(size = 20 , face = 'bold') ,
        axis.text.x = element_text(size = 15) ,
        axis.text.y = element_text(size = 15) ,
        legend.title = element_text(size = 12 , face = "bold") ,
        legend.text = element_text(size = 12)
      )

    print(p)
    plc_periods <- NULL
    for (i in 1:plc_lsp_num)
    {

      plc_periods <- c(plc_periods , clean_periods$Period[i])

    }

    return(list(Periods = plc_periods , LSP_Data = lomb))


  }
  print("lombScarg loaded!")

  viewLC <- function(lcs , lc_target , field_name)
  {

    library(ggplot2)

    tar_star <- lcs[[as.character(lc_target)]] #Get target

    tar_star$NumericTime <- as.numeric(
      difftime(tar_star$TimeStamp , min(tar_star$TimeStamp) , units = 'days')
    )

    #First we will plot the raw LC and systematic
    raw <- ggplot(data = tar_star , aes(color = Date)) +
      geom_point(aes(x = NumericTime , y = Magnitude , shape = "Raw Magnitude")) +
      geom_point(aes(x = NumericTime , y = AMag , shape = "Adjusted Magnitude")) +
      scale_y_reverse() +
      scale_color_discrete() +
      scale_shape_manual(
        values = c(
          "Raw Magnitude" = 4 ,
          "Adjusted Magnitude" = 16
        )
      ) +
      labs(title = "Running Time Light Curve" ,
           subtitle = paste0(field_name , " " , tar_star$StarID[1]) ,
           x = "Running Time [days]" ,
           y = "Instrumental Magnitude" ,
           color = "Observation Date" ,
           shape = "Magnitude Type"
           ) +
      theme(
            plot.title = element_text(size = 25 , face = "bold") ,
            plot.subtitle = element_text(size = 20) ,
            axis.title.x = element_text(size = 20 , face = "bold") ,
            axis.title.y = element_text(size = 20 , face = 'bold') ,
            axis.text.x = element_text(size = 15) ,
            axis.text.y = element_text(size = 15) ,
            legend.title = element_text(size = 12 , face = "bold") ,
            legend.text = element_text(size = 12)
          )

    print(raw)

    if (tar_star$SysRemoved[1] == T)
    {
      observation_number <- seq(from = 1 , to = length(tar_star$SMag))

      sys <- ggplot(data = tar_star , aes(color = Date)) +
        geom_line(aes(x = observation_number , y = SMag)) +
        scale_y_reverse() +
        scale_color_discrete() +
        geom_hline(aes(yintercept = 0)) +
        labs(title = "Systematic Light Curve" ,
             subtitle = paste0(field_name , " " , tar_star$StarID[1]) ,
             x = "Observation Number" ,
             y = "Instrumental Magnitude" ,
             color = "Observation Date"
             ) +
        theme(
              plot.title = element_text(size = 25 , face = "bold") ,
              plot.subtitle = element_text(size = 20) ,
              axis.title.x = element_text(size = 20 , face = "bold") ,
              axis.title.y = element_text(size = 20 , face = 'bold') ,
              axis.text.x = element_text(size = 15) ,
              axis.text.y = element_text(size = 15) ,
              legend.title = element_text(size = 12 , face = "bold") ,
              legend.text = element_text(size = 12)
            )

      print(sys)

    }

    return(tar_star)

  }
  print("viewLC loaded!")

  makePhasesCurve <- function(lcs , plc_target , plc_period , field_name , use_fitted , use_lsp , fitted_P , lsp_P , was_fitting_run , was_lsp_run)
  {

    library(ggplot2)

    tar_star <- lcs[[as.character(plc_target)]] #Get target

    tar_star$NumericTime <- as.numeric(
      difftime(tar_star$TimeStamp , min(tar_star$TimeStamp) , units = 'days')
    ) #t - T for all data points of tar_star


    if ((use_fitted == T && was_fitting_run == F) || (use_lsp == T && was_lsp_run == F))
    {

      print("Error! Either the fitting function or LSP function was not run and no calculated period can be used!")
      print("plc_period will be used instead!!!")
      cur_period <- plc_period
      prefix <- NULL

    }
    else if (use_fitted == T & was_fitting_run == T)
    {

      cur_period <- fitted_P
      prefix <- "Fitted"

    }
    else if (use_lsp == T && was_lsp_run == T)
    {

      cur_period <- lsp_P

      #If the LSP option is chosen, the PLC function will plot for the top 3 periods found by the LSP
      for (i in 1:length(cur_period))
      {

        tar_star$Phase <- (tar_star$NumericTime / cur_period[i]) %% 1 #Create phase column
        p <- ggplot(data = tar_star, aes(x = Phase , y = AMag , color = Date)) +
          geom_point(alpha = .6) +
          scale_y_reverse() +
          labs(
            title = paste0("Phased Light Curve for " , field_name , " " , plc_target),
            subtitle = paste0("LSP Period " , i , ": " , cur_period[i] , " days") ,
            x = "Phase",
            y = "Instrumental Magnitude",
            color = "Observation Night"
          ) +
          scale_color_discrete() +
          theme(
            plot.title = element_text(size = 25 , face = "bold") ,
            plot.subtitle = element_text(size = 20) ,
            axis.title.x = element_text(size = 20 , face = "bold") ,
            axis.title.y = element_text(size = 20 , face = 'bold') ,
            axis.text.x = element_text(size = 15) ,
            axis.text.y = element_text(size = 15) ,
            legend.title = element_text(size = 12 , face = "bold") ,
            legend.text = element_text(size = 12)
          )

        print(p)





      }

      return(tar_star)

    }
    else if (use_lsp == F && use_fitted == F)
    {

      cur_period <- plc_period
      prefix <- NULL


    }

    #Create phase column
    tar_star$Phase <- (tar_star$NumericTime / cur_period) %% 1 #Create phase column

    p <- ggplot(data = tar_star, aes(x = Phase , y = AMag , color = Date)) +
      geom_point(alpha = .6) +
      scale_y_reverse() +
      labs(
        title = paste0("Phased Light Curve for " , field_name , " " , plc_target),
        subtitle = paste0(prefix , " Period: " , cur_period , " days") ,
        x = "Phase",
        y = "Instrumental Magnitude",
        color = "Observation Night"
      ) +
      scale_color_discrete() +
      theme(
        plot.title = element_text(size = 25 , face = "bold") ,
        plot.subtitle = element_text(size = 20) ,
        axis.title.x = element_text(size = 20 , face = "bold") ,
        axis.title.y = element_text(size = 20 , face = 'bold') ,
        axis.text.x = element_text(size = 15) ,
        axis.text.y = element_text(size = 15) ,
        legend.title = element_text(size = 12 , face = "bold") ,
        legend.text = element_text(size = 12)
      )

    print(p)
    return(tar_star)

  }
  print("makePhasesCurve loaded!")

  jstetcalc <- function(df)
  {

    mags <- df$AMag
    mags <- sort(mags)
    counts <- 10^((-(mags - 25))/2.5)
    er <- 1 / sqrt(counts)
    rms <- sd(mags)
    nms <- length(mags)

    cut <- as.integer(.1 * nms)
    d90 <- abs(mags[cut + 1] - mags[length(mags) - cut])

    Jt <- rep(0.0 , nms)
    Jb <- rep(0.0 , nms)

    wk <- 1
    MeanMag <- mean(mags)

    for (i in seq(from = 1 ,to = nms - 2 , by = 2))
      {

        Sigi <- (mags[i] - MeanMag) / (er[i])*(sqrt(nms / (nms - 1)))
        Sigj <- (mags[i+1] - MeanMag) / (er[i+1])*(sqrt(nms / (nms - 1)))

        Pk <- Sigi*Sigj

        if (Pk > 0.0)
        {
          sgnPk <- 1.0
        }
        if (Pk == 0.0)
        {
          sgnPk <- 0.0
        }
        if (Pk < 0.0)
        {
          sgnPk <- -1.0
        }

        Jt[i] <- wk*sgnPk*(sqrt(abs(Pk)))
        Jb[i] <- (wk)

      }

    jstetson <- sum(Jt) / sum(Jb)
    result <- list(Jstetson = jstetson , D90 = d90)

    # print("Returning JStet result")
    # print(result)
    return(result)



  }
  print("jstetcalc loaded!")

  runVarAnalysis <- function(lcs , field_name)
  {

    library(magicaxis)
    library(ggplot2)

    #We will run calculations to see if stars are variable
    id_list <- NULL
    j_list <- NULL
    d90_list <- NULL


    for (star in lcs)
    {

      #Obtain variability results
      var_results <- jstetcalc(star)

      #Save needed data
      id_list <- c(id_list , star$StarID[1])
      j_list <- c(j_list , var_results$Jstetson)
      d90_list <- c(d90_list , var_results$D90)

    }

    #We now calculate medians and thresholds
    sigma_clipped_j_list <- magclip(j_list , sigma = 5)$x
    sigma_clipped_D90_list <- magclip(d90_list , sigma = 5)$x

    median_J <- median(j_list)
    median_d90 <- median(d90_list)

    sd_J <- sd(sigma_clipped_j_list)
    #sd_d90 <- sd(sigma_clipped_D90_list)

    threshJ <- median_J + 2 * sd_J
    #threshd90 <- median_d90 + 2 * sd_d90

    #We construct our variability dataframe
    var_df <- data.frame(StarID = id_list , Stetson_J = j_list , D90 = d90_list)


    var_df$J_Flag <- var_df$Stetson_J > threshJ #Flag stars that are above the threshold J Value
    var_df$D90_Flag <- F #Intialize this column as booleans

    ###To flag the stars that pass the D90 requirements, we will bin stars based on magnitudes, find the mean D90, and then flag anything above 2sigma
    mean_mags <- sapply(lcs , function(df) mean(df$AMag)) #Get mean mags of the stars
    # mag_order <- order(mean_mags , decreasing = F)
    # ordered_by_mag_lcs <- lcs[mag_order]

    #Go through lcs based on each 1% of mags (dimmest to brightest)
    percentiles <- quantile(mean_mags , probs = seq(0 , 1 , .01))
    bins <- cut(mean_mags , breaks = percentiles , include.lowest = T , labels = F)
    d90_quartile_threshes <- list()
    for (bin_index in sort(unique(bins)))
    {

      #Get this percentile of lcs
      lc_in_this_bin <- lcs[bins == bin_index]

      #Calculate the threshold D90 for this grouping of LCs
      ids_this_bin <- sapply(lc_in_this_bin , function(df) return(df$StarID[1]))
      var_df_this_bin <- var_df[var_df$StarID %in% ids_this_bin , ]
      mediand90 <- median(var_df_this_bin$D90)
      sdd90 <- sd(var_df_this_bin$D90)
      threshd90 <- mediand90 + 2 * sdd90

      #Update the D90 flag status of the stars in this bin based on threshd90
      var_df[var_df$StarID %in% ids_this_bin ,]$D90_Flag <- var_df[var_df$StarID %in% ids_this_bin , ]$D90 > threshd90
    }





    #We create a histogram of J values and D90 values
    Jplot <- ggplot(data = var_df , aes(x = Stetson_J)) +
      geom_histogram(binwidth = .5) +
      geom_vline(xintercept = threshJ , color = 'red' , linetype = 'solid' , size = .7) +
      labs(title = paste0("Stetson J Values for " , field_name) ,
           subtitle = paste0("Mean J: " , mean(j_list) , " | Median J: " , median(j_list) , " | StDev: " , sd(j_list) ,
                             "\n Threshold J: " , threshJ) ,
           x = "J Value" ,
           y = "Count"
      ) +
      xlim(0 , 40) +
      theme(
        plot.title = element_text(size = 25 , face = "bold") ,
        plot.subtitle = element_text(size = 20) ,
        axis.title.x = element_text(size = 20 , face = "bold") ,
        axis.title.y = element_text(size = 20 , face = 'bold') ,
        axis.text.x = element_text(size = 15) ,
        axis.text.y = element_text(size = 15) ,
        legend.title = element_text(size = 12 , face = "bold") ,
        legend.text = element_text(size = 12)
      )
    d90plot <- ggplot(data = var_df , aes(x = D90)) +
      geom_histogram(binwidth = .005) +
      labs(title = paste0("D90 Values for " , field_name) ,
           subtitle = paste0("Mean D90: " , mean(d90_list) , " | Median D90: " , median(d90_list) , " | StDev: " , sd(d90_list)) ,
           x = "D90 Value" ,
           y = "Count"
      ) +
      theme(
        plot.title = element_text(size = 25 , face = "bold") ,
        plot.subtitle = element_text(size = 20) ,
        axis.title.x = element_text(size = 20 , face = "bold") ,
        axis.title.y = element_text(size = 20 , face = 'bold') ,
        axis.text.x = element_text(size = 15) ,
        axis.text.y = element_text(size = 15) ,
        legend.title = element_text(size = 12 , face = "bold") ,
        legend.text = element_text(size = 12)
      )

    print(Jplot)
    print(d90plot)
    return(var_df)

  }
  print("runVarAnalysis loaded!")

  makeRMSplot <- function(lcs , var_data , field_name)
  {
    library(magicaxis)
    library(ggplot2)

    #We will plot RMS vs. magnitude
    id_list <- NULL
    rms_list <- NULL
    med_mag_list <- NULL
    status_list <- NULL

    for (star in lcs)
    {

      median_mag <- median(star$AMag)
      sigma_clipped_mags <- magclip(star$AMag , sigma = 5)$x
      mag_RMS <- sd(sigma_clipped_mags)

      id_list <- c(id_list , star$StarID[1])
      med_mag_list <- c(med_mag_list , median_mag)
      rms_list <- c(rms_list , mag_RMS)

      #Create a variability status if the star has a J flag, D90 flag, or both
      status <- "No Var Flag"

      if (var_data[var_data$StarID == star$StarID[1] , ]$J_Flag == T && var_data[var_data$StarID == star$StarID[1] , ]$D90_Flag == T)
      {
        status <- "J + D90 Flag"
      }

      else if (var_data[var_data$StarID == star$StarID[1] , ]$J_Flag == T)
      {
        status <- "J Flag"
      }

      else if (var_data[var_data$StarID == star$StarID[1] , ]$D90_Flag == T)
      {
        status <- "D90 Flag"
      }

      status_list <- c(status_list , status)



    }

    #Create a df to hold the info
    rms_df <- data.frame(
      StarID = id_list ,
      Median_Mag = med_mag_list ,
      RMS = rms_list ,
      Variability_Status = status_list ,
      stringsAsFactors = F
    )

    #Plot RMS vs. Mag
    p <- ggplot(data = rms_df , aes( x = Median_Mag , y = RMS , color = Variability_Status)) +
      geom_point() +
      scale_color_manual( values = c('No Var Flag' = 'black' , 'J + D90 Flag' = 'purple' , 'J Flag' = 'red' , 'D90 Flag' = 'blue' )) +
      scale_y_log10() +
      labs(
        title = paste0("RMS vs. Median Instrumental Magnitude for " , field_name) ,
        x = "Instrumental Magnitude" ,
        y = "RMS" ,
        color = "Possible Variable"
      ) +
      theme(
        plot.title = element_text(size = 25 , face = "bold") ,
        axis.title.x = element_text(size = 20 , face = "bold") ,
        axis.title.y = element_text(size = 20 , face = 'bold') ,
        axis.text.x = element_text(size = 15) ,
        axis.text.y = element_text(size = 15) ,
        legend.title = element_text(size = 12 , face = "bold") ,
        legend.text = element_text(size = 12)
      )

    print(p)

    return(rms_df)






  }
  print("makeRMSplot loaded!")

  TOROSObsPredict <- function(lcs , target , N , sigma_level , target_name , minp , maxp)
  {

    library(magicaxis)
    library(ggplot2)

    #This sim will find minimum number of TOROS observations required to recover the target star's variability
    target_star <- lcs[[as.character(target)]]
    lcs[[as.character(target)]] <- NULL #Removes target from LCs to avoid contamination

    thresh_list <- NULL
    J_Stet_Comp_Table <- data.frame()

    for (i in 3:nrow(target_star))
    {

      print(paste0("Sample size i = " , i))

      #We loop through lightcurves to compute the threshold J level at this sample size i
      print("Computing threshold J")
      comp_J_list <- NULL

      #Compute threshold J value for this i sample size
      for (lc in lcs)
      {

        print(paste0("ThreshJ, analyzing Star " , lc$StarID[1]))

        #Sample each star with sample size i
        sample_lc_ind <- sample(seq_len(nrow(target_star)) , size = i)
        cut_lc <- lc[sample_lc_ind , ]

        #Calculate the J value of this sample
        JD90 <- jstetcalc(cut_lc)
        J <- JD90[[1]] #1 gives J, 2 gives D90
        comp_J_list <- c(comp_J_list , J)

      }

      #Compute threshold for this i
      cur_med_J <- median(comp_J_list)
      cur_sd_J <- sd(magclip(comp_J_list , sigma = 5)$x)

      cur_thresh_J <- cur_med_J + (cur_sd_J * sigma_level)
      thresh_list <- c(thresh_list , cur_thresh_J)


      #Now we do the sampling from our target star N times
      print("Computing target Js")
      tar_J_list <- numeric(N)

      for (j in 1:N)
      {

        #Same process as before
        sample_tar_ind <- sample(seq_len(nrow(target_star)) , size = i)
        cut_tar_lc <- target_star[sample_tar_ind , ]

        tarJD90 <- jstetcalc(cut_tar_lc)
        tarJ <- tarJD90[[1]]
        tar_J_list[j] <- tarJ

      }

      J_Stet_Comp_Table <- rbind(J_Stet_Comp_Table , data.frame(t(as.matrix(tar_J_list))))

    }

    #We will bind the threshold values to the table
    J_Stet_Comp_Table <- cbind(J_Stet_Comp_Table , thresh_list)

    #Lets name our columns and rows
    colnames(J_Stet_Comp_Table) <- paste0("Sample Number " , 1:N)
    colnames(J_Stet_Comp_Table)[ncol(J_Stet_Comp_Table)] <- "Threshold J"
    rownames(J_Stet_Comp_Table) <- paste0("Sample Size " , 3:nrow(target_star))

    #We now will create our plot, first making a count of samples passed
    counts_passed <- NULL
    sample_df <- J_Stet_Comp_Table

    for (ii in 1:nrow(sample_df))
    {

      sample_df[ ii , ] <- as.numeric(sample_df[ ii , 1:ncol(sample_df)]) > sample_df[ ii , ncol(sample_df)] #Checks if J val in each column is greater/less than threshold in last column.
      counts_passed <- c(counts_passed , sum(sample_df[ ii , ]))

    }

    sample_sizes <- seq(from = 3 , to = nrow(sample_df) + 2)
    new_df <- data.frame(Sample_Sizes = sample_sizes , Passed = counts_passed , Percent_Passed = (counts_passed) / N)

    #Now we plot percent passing vs. sample size, effectively plotting the number of TOROS observations needed to confirm the variability of your target star
    p <- ggplot(data = new_df , aes(x = Sample_Sizes)) +
      geom_line(aes(y = Percent_Passed) , size = .8) +
      geom_hline(aes(yintercept = 1) , color = 'red' , linetype = 'dotted' , size = .9) +
      labs(
        title = paste0("% of Passing Samples vs. Sample Size for " , target_name) ,
        x = "Sample Size (Number of Observations)" ,
        y = "Percentage Passeed"
      ) +
      xlim(minp , maxp) +
      theme(
        plot.title = element_text(size = 25 , face = "bold") ,
        #plot.subtitle = element_text(size = 20) ,
        axis.title.x = element_text(size = 20 , face = "bold") ,
        axis.title.y = element_text(size = 20 , face = 'bold') ,
        axis.text.x = element_text(size = 15) ,
        axis.text.y = element_text(size = 15) ,
        legend.title = element_text(size = 12 , face = "bold") ,
        legend.text = element_text(size = 12)
      )

    print(p)

    LSSTPred1 <- list(J_Comp = J_Stet_Comp_Table , Sample_Data = new_df)

    return (LSSTPred1)



  }
  print("TOROSObsPredict loaded!")

  TOROSTimePredict <- function(lcs , target , target_name , I , M)
  {
    library(dplyr)
    library(ggplot2)

    #Mark target star and remove from lcs
    target_star <- lcs[[as.character(target)]]
    lcs[[as.character(target)]] <- NULL

    J_list <- NULL
    tint_list <- NULL

    #We employ the same method with a known sample size and record J vs. Time Interval
    for (j in 1:M)
    {

      print(paste0("Iteration " , j))

      #We take our sample
      sample_tar_ind <- sample(seq_len(nrow(target_star)) , size = I)
      cut_tar_lc <- target_star[sample_tar_ind , ]
      cut_tar_lc <- cut_tar_lc %>%
        arrange(TimeStamp)

      #We calculate J
      JD90 <- jstetcalc(cut_tar_lc)
      J <- JD90[[1]]
      J_list <- c(J_list , J)

      #We find the time interval between earliest and latest randomly selected lc point
      time_int <- cut_tar_lc$TimeStamp[nrow(cut_tar_lc)] - cut_tar_lc$TimeStamp[1]
      time_int <- as.numeric(time_int , units = 'days')
      tint_list <- c(tint_list , time_int)

    }

    #Build all the data into a frame
    J_timeint_df <- data.frame(Time_Interval = tint_list , J = J_list)

    #We will truncate times since we only care about if observation need to stretch over multiple nights of observation
    J_timeint_df$Time_Interval <- trunc(J_timeint_df$Time_Interval , 2)

    #Now we will obtain needed statistical measures for our plot
    uniq_tint <- unique(J_timeint_df$Time_Interval)
    Jm_list <- NULL
    Jsd_list <- NULL
    Jtime_list <- NULL

    for (time in uniq_tint)
    {

      cut_J_df <- J_timeint_df[J_timeint_df$Time_Interval == time , ]
      if (nrow(cut_J_df) < 2) {next}

      J_mean <- mean(cut_J_df$J)
      J_sd <- sd(cut_J_df$J)
      J_time <- cut_J_df$Time_Interval[1]

      Jm_list <- c(Jm_list , J_mean)
      Jsd_list <- c(Jsd_list , J_sd)
      Jtime_list <- c(Jtime_list , J_time)

    }

    new_Jtime_df <- data.frame(Time_Interval = Jtime_list , J_mean = Jm_list , J_sd = Jsd_list)

    p <- ggplot(data = new_Jtime_df , aes(x = Time_Interval , y = J_mean)) +
      geom_point(color = 'red') +
      geom_errorbar(aes(ymin = J_mean - J_sd , ymax = J_mean + J_sd) , width = .01) +
      labs(
        title = paste0("J vs. Time Interval of Observations for " , target_name) ,
        x = "Time Interval [days]" ,
        y = "Calulcated J Value" ,
        subtitle = paste0("Sample Size: " , I , " | Iterations: " , M)
      ) +
      theme(
        plot.title = element_text(size = 25 , face = "bold") ,
        plot.subtitle = element_text(size = 20) ,
        axis.title.x = element_text(size = 20 , face = "bold") ,
        axis.title.y = element_text(size = 20 , face = 'bold') ,
        axis.text.x = element_text(size = 15) ,
        axis.text.y = element_text(size = 15) ,
        legend.title = element_text(size = 12 , face = "bold") ,
        legend.text = element_text(size = 12)
      )

    print(p)

    LSSTPred2 <- list(Time_Int_Data = new_Jtime_df)
    return(LSSTPred2)


  }
  print("TOROSTimePredict loaded!")

  #Initialize these so program doesnt crash
  fit_P <- 1
  lsp_P <- 1
  Extra_Info <- list()

  #Actually run the functions now based on user requests and inputs
  if (preload_data == F)
  {

    print("Running loadInNights...")
    field_data <- loadInNights(night_folders = night_folders)
    Extra_Info$Night_Data <- field_data

    print("Running makeStarStat...")
    star_stats <- makeStarStat(data_list = field_data$Photometry , id_list = field_data$CommonIDs)
    print("Running systematicRemoval...")
    sys_Rem_LCs <- systematicRemoval(lcs = star_stats)
    ExtraInfo$Star_LCs <- sys_Rem_LCs

    print("Running runVarAnalysis...")
    var_data <- runVarAnalysis(lcs = sys_Rem_LCs , field_name = field_name)
    print("Running makeRMSPlot...")
    rms_data <- makeRMSplot(lcs = sys_Rem_LCs , field_name = field_name , var_data = var_data)

    if (!(target %in% field_data$CommonIDs) && target != 0)
    {

      print("Error! Your selected star does not have sufficient data with which to do analysis!")
      return()

    }

    if (plot_lcs == T)
    {

      print("Running makeLightCurves...")
      makeLightCurves(lcs = sys_Rem_LCs , lc_folder = lc_folder ,field_name = field_name)

    }

    if (plot_lc == T)
    {

      print("Running viewLC...")
      viewLC(lcs = sys_Rem_LCs , field_name = field_name , lc_target = target)

    }

    if (fit_data == T)
    {

      print("Running fitCurveToData...")
      fit_list <- fitCurveToData(lcs = sys_Rem_LCs , fC_target = target , field_name = field_name , guess = fC_guess)
      Extra_Info$Fitting_Data <- fit_list

    }

    if (plot_LSP == T)
    {

      print("Running lombScarg")
      lsp <- lombScarg(lcs = sys_Rem_LCs , lsp_target = target , field_name = field_name , plc_lsp_num = plc_LSP_per_num)
      Extra_Info$LSP_Data <- lsp$LSP_Data

    }

    if(run_Var_Analysis == T)
    {

      print("Running runVarAnalysis...")
      var_data <- runVarAnalysis(lcs = sys_Rem_LCs , field_name = field_name)
      Extra_Info$Variability <- var_data

    }

    if (run_RMS_Analysis == T)
    {

      print("Running makeRMSplot...")
      rms_data <- makeRMSplot(lcs = sys_Rem_LCs , field_name = field_name , var_data = var_data)
      Extra_Info$Error_Analysis <- rms_data

    }

    if(plot_plc == T)
    {

      print("Running makePhasesCurve...")
      makePhasesCurve(lcs = sys_Rem_LCs , plc_target = target , plc_period = plc_period , field_name = field_name , use_fitted = plc_use_fitted_P , use_lsp = plc_use_LSP_P ,fitted_P = fit_list$Fit_P , lsp_P = lsp$Periods , was_fitting_run = fit_data , was_lsp_run = plot_LSP)

    }

    if(write_lcs == T)
    {

      print("Running writeLcCsvs...")
      writeLcCsvs(lcs = sys_Rem_LCs , field_name = field_name , lccsv_folder = lccsv_folder)

    }

    if(run_LSST_test == T)
    {

      print("Running TOROSObsPredict...")
      obs_pred <- TOROSObsPredict(lcs = sys_Rem_LCs , target = target , target_name = LSST_target_name , N = LSST_N , sigma_level = LSST_sigma , minp = obs_plot_min , maxp = obs_plot_max)
      Extra_Info$Observation_Test <- obs_pred

      #Get first N to achieve 100%
      first_index_N <- which(obs_pred$Sample_Data$Percent_Passed == 1)[1]
      opt_I <- obs_pred$Sample_Data$Sample_Sizes[first_index_N]

      print("Running TOROSTimePredict...")
      time_pred <- TOROSTimePredict(lcs = sys_Rem_LCs , target = target , target_name = LSST_target_name , M = LSST_M , I = opt_I)
      Extra_Info$Interval_Test <- time_pred

    }

    return (Extra_Info)

  }

  if (preload_data == T)
  {

    if (!(target %in% preloaded_data$Night_Data$CommonIDs) && target !=0)
    {

      print("Error! Your selected star does not have sufficient data with which to do analysis!")
      return()

    }

    preloaded_lcs <- preloaded_data$Star_LCs
    Extra_Info$Night_Data <- preloaded_data$Night_Data
    Extra_Info$Star_LCs <- preloaded_data$Star_LCs

    if (plot_lcs == T)
    {

      print("Running makeLightCurves...")
      makeLightCurves(lcs = preloaded_lcs , lc_folder = lc_folder ,field_name = field_name)

    }

    if (plot_lc == T)
    {

      print("Running viewLC...")
      viewLC(lcs = preloaded_lcs , field_name = field_name , lc_target = target)

    }

    if (fit_data == T)
    {

      print("Running fitCurveToData...")
      fit_list <- fitCurveToData(lcs = preloaded_lcs , fC_target = target , field_name = field_name , guess = fC_guess)
      Extra_Info$Fitting_Data <- fit_list

    }

    if (plot_LSP == T)
    {

      print("Running lombScarg")
      lsp <- lombScarg(lcs = preloaded_lcs , lsp_target = target , field_name = field_name , plc_lsp_num = plc_LSP_per_num)
      Extra_Info$LSP_Data <- lsp$LSP_Data


    }

    if(run_Var_Analysis == T)
    {

      print("Running runVarAnalysis...")
      var_data <- runVarAnalysis(lcs = preloaded_lcs , field_name = field_name)
      Extra_Info$Variability <- var_data

    }

    if (run_RMS_Analysis == T)
    {

      print("Running makeRMSplot...")
      rms_data <- makeRMSplot(lcs = preloaded_lcs , field_name = field_name , var_data = var_data)
      Extra_Info$Error_Analysis <- rms_data

    }

    if (plot_plc == T)
    {

      print("Running makePhasesCurve...")
      makePhasesCurve(lcs = preloaded_lcs , plc_target = target , plc_period = plc_period , field_name = field_name , use_fitted = plc_use_fitted_P , use_lsp = plc_use_LSP_P ,fitted_P = fit_list$Fit_P , lsp_P = lsp$Periods , was_fitting_run = fit_data , was_lsp_run = plot_LSP)

    }

    if(write_lcs == T)
    {

      print("Running writeLcCsvs...")
      writeLcCsvs(lcs = preloaded_lcs , field_name = field_name , lccsv_folder = lccsv_folder)

    }

    if(run_LSST_test == T)
    {

      print("Running TOROSObsPredict...")
      obs_pred <- TOROSObsPredict(lcs = preloaded_lcs , target = target , target_name = LSST_target_name , N = LSST_N , sigma_level = LSST_sigma , minp = obs_plot_min , maxp = obs_plot_max)
      Extra_Info$Observation_Test <- obs_pred

      #Get first N to achieve 100%
      first_index_N <- which(obs_pred$Sample_Data$Percent_Passed == 1)[1]
      opt_I <- obs_pred$Sample_Data$Sample_Sizes[first_index_N]

      print("Running TOROSTimePredict...")
      int_pred <- TOROSTimePredict(lcs = preloaded_lcs , target = target , target_name = LSST_target_name , M = LSST_M , I = opt_I)
      Extra_Info$Interval_Test <- int_pred

    }



    return (list(Extra_Info))

  }

}
