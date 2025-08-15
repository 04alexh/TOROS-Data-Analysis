TOROSAnalysis <- function(night_folders = NULL,
                          preload_lcs = F ,
                          preloaded_lcs = NULL ,
                          write_lcs = F ,
                          lccsv_folder = "" ,
                          plot_lcs = F ,
                          lc_folder = "" ,
                          field_name = "Observed Field" ,
                          plot_plc = F ,
                          plc_target = 0 ,
                          plc_period = 0.0 ,
                          run_obs_pred = F ,
                          run_time_pred = F ,
                          LSST_target = 0 ,
                          LSST_target_name = "" ,
                          LSST_N = 0 ,
                          LSST_sigma = 0 ,
                          LSST_M = 0,
                          LSST_I = 0
                          ) #field_name , file_prefix = "Alignened")
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

  systematicRemoval <- function(lcs) #, method)
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

        if (delta_pos > 0 && delta_pos < 2500 && delta_mag < 1.5 && sc_rms < .05)
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
              x = "Time" ,
              y = "Instrumental Magnitude" ,
              title = paste0("Lightcurve for " , field_name , star$StarID[1]) ,
              subtitle = paste0("Systematic Removal: " , star$SysRemoved[1]) ,
              shape = "Legend" ) +
        scale_color_discrete() +
        scale_shape_manual(values = c(
          "Raw Mag" = 4,
          "Adj Mag" = 16
        ))

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
          scale_color_discrete()

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

  makePhasesCurve <- function(lcs , plc_target , plc_period , field_name)
  {

    library(ggplot2)

    tar_star <- lcs[[as.character(plc_target)]]

    tar_star$NumericTime <- as.numeric(
      difftime(tar_star$TimeStamp , min(tar_star$TimeStamp) , units = 'days')
    )

    tar_star$Phase <- (tar_star$NumericTime / plc_period) %% 1


    p <- ggplot(data = tar_star, aes(x = Phase , y = AMag , color = Date)) +
      geom_point(alpha = .6) +
      scale_y_reverse() +
      labs(
        title = paste0("Phased Light Curve for " , field_name , plc_target),
        x = "Phase",
        y = "Instrumental Magnitude",
        color = "Observation Night"
      ) +
      scale_color_discrete()

    print(p)

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
    # return(result)



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
    sd_d90 <- sd(sigma_clipped_D90_list)

    threshJ <- median_J + 2 * sd_J
    threshd90 <- median_d90 + 2 * sd_d90

    #We construct our variability dataframe
    var_df <- data.frame(StarID = id_list , Stetson_J = j_list , D90 = d90_list)


    var_df$J_Flag <- var_df$Stetson_J > threshJ
    var_df$D90_Flag <- var_df$D90 > threshd90


    #We create a histogram of J values and D90 values
    Jplot <- ggplot(data = var_df , aes(x = Stetson_J)) +
      geom_histogram(binwidth = .5) +
      geom_vline(xintercept = threshJ , color = 'red' , linetype = 'solid' , size = .7) +
      labs(title = paste0("Stetson J Values for " , field_name) ,
           x = "J Value" ,
           y = "Count"
      ) +
      xlim(0 , 40)
    d90plot <- ggplot(data = var_df , aes(x = D90)) +
      geom_histogram(binwidth = .005) +
      geom_vline(xintercept = threshd90 , color = 'red' , linetype = 'solid' , size = .7) +
      labs(title = paste0("D90 Values for " , field_name) ,
           x = "D90 Value" ,
           y = "Count"
      )

    print(Jplot)
    print(d90plot)
    return(var_df)

  }
  print("runVarAnalysis loaded!")

  makeRMSplot <- function(lcs , var_data ,field_name)
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
      Variability_Status = status_list
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
      )

    print(p)

    return(rms_df)






  }
  print("makeRMSplot loaded!")

  TOROSObsPredict <- function(lcs , target , N , sigma_level , target_name)
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
      xlim(0 , 20)

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
      )

    print(p)

    LSSTPred2 <- list(Time_Int_Data = new_Jtime_df)
    return(LSSTPred2)


  }
  print("TOROSTImePredict loaded!")




  #Actually run the functions now based on user requests and inputs
  if (preload_lcs == F)
  {

    print("Running loadInNights...")
    field_data <- loadInNights(night_folders = night_folders)
    print("Running makeStarStat...")
    star_stats <- makeStarStat(data_list = field_data$Photometry , id_list = field_data$CommonIDs)
    # return (star_stats)
    print("Running systematicRemoval...")
    sys_Rem_LCs <- systematicRemoval(lcs = star_stats)
    print("Running runVarAnalysis...")
    var_data <- runVarAnalysis(lcs = sys_Rem_LCs , field_name = field_name)
    print("Running makeRMSPlot...")
    rms_data <- makeRMSplot(lcs = sys_Rem_LCs , field_name = field_name , var_data = var_data)

    if (plot_lcs == T)
    {

      print("Running makeLightCurves...")
      makeLightCurves(lcs = sys_Rem_LCs , lc_folder = lc_folder ,field_name = field_name)

    }

    if(plot_plc == T)
    {

      print("Running makePhasesCurve...")
      makePhasesCurve(lcs = sys_Rem_LCs , plc_target = plc_target , plc_period = plc_period , field_name = field_name)

    }

    if(write_lcs == T)
    {

      print("Running writeLcCsvs...")
      writeLcCsvs(lcs = sys_Rem_LCs , field_name = field_name , lccsv_folder = lccsv_folder)

    }

    if(run_obs_pred == T)
    {

      print("Running TOROSObsPredict...")
      TOROSObsPredict(lcs = sys_Rem_LCs , target = LSST_target , target_name = LSST_target_name , N = LSST_N , sigma_level = LSST_sigma)

    }

    if(run_time_pred == T)
    {

      print("Running TOROSTimePredict...")
      TOROSTimePredict(lcs = sys_Rem_LCs , target = LSST_target , target_name = LSST_target_name , M = LSST_M , I = LSST_I)

    }

    return (list(Night_Data = field_data , Star_LCs = sys_Rem_LCs , Variability = var_data , Error_Analysis = rms_data))

  }

  if (preload_lcs == T)
  {

    if (plot_lcs == T)
    {

      print("Running makeLightCurves...")
      makeLightCurves(lcs = preloaded_lcs , lc_folder = lc_folder ,field_name = field_name)

    }

    if (plot_plc == T)
    {

      print("Running makePhasesCurve...")
      makePhasesCurve(lcs = preloaded_lcs , plc_target = plc_target , plc_period = plc_period , field_name = field_name)

    }

    if(write_lcs == T)
    {

      print("Running writeLcCsvs...")
      writeLcCsvs(lcs = preloaded_lcs , field_name = field_name , lccsv_folder = lccsv_folder)

    }

    if(run_obs_pred == T)
    {

      print("Running TOROSObsPredict...")
      TOROSObsPredict(lcs = preloaded_lcs , target = LSST_target , target_name = LSST_target_name , N = LSST_N , sigma_level = LSST_sigma)

    }

    if(run_time_pred == T)
    {

      print("Running TOROSTimePredict...")
      TOROSTimePredict(lcs = preloaded_lcs , target = LSST_target , target_name = LSST_target_name , M = LSST_M , I = LSST_I)

    }

    print("Running runVarAnalysis...")
    var_data <- runVarAnalysis(lcs = preloaded_lcs , field_name = field_name)
    print("Running makeRMSplot...")
    rms_data <- makeRMSplot(lcs = preloaded_lcs , field_name = field_name , var_data = var_data)
    return (list(Variability = var_data , Error_Analysis = rms_data))

  }

}
