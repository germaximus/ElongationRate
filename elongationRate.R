#############################################################################################################
###########  Prepare mRNA coverage data for analysis  #######################################################
#############################################################################################################

library(magrittr)
setwd("f:/Ribosome_profiling/Mammals/Mouse/Injections/Translation_rate/")

      #Instuction how to load Riso-seq coverage data for every sample, separated by individual reads. Takes a lot of time, therefore I saved the result as Rdata file.
      # setwd("f:/Ribosome_profiling/Mammals/Mouse/Injections/Translation_rate/")
      # library("R.utils")
      # coverage_samples <- vector("list", length=length(list.files(path="./sample_coverage", pattern='coverage')))
      # names(coverage_samples) <- list.files(path="./sample_coverage", pattern='coverage')
      # for(i in 1:length(coverage_samples)) {
      #   handle <- file(paste0(getwd(), "/sample_coverage/", names(coverage_samples)[i]), open = "rb"); linecount <- countLines(handle);  close(handle);
      #   handle <- file(paste0(getwd(), "/sample_coverage/", names(coverage_samples)[i]), open = "r");
      #   file <- vector("list", length=linecount)
      #   for(a in 1:linecount) {
      #     line <- readLines(handle, n=1)
      #     line <- unlist(strsplit(line, "\t"))
      #     names(file)[a] <- line[1]
      #     file[[a]] <- as.integer(line[-1])
      #   }
      #   close(handle)
      #   # in early version of the Coverage.pl genes were written in a random order (as unsorted hash keys). To fix this, sort file names:
      #   file <- file[sort(names(file))]
      #   coverage_samples[[i]] <- file
      # }
      # save(coverage_samples, file = "./Rdata/coverage_samples2.Rdata")
      
load(file="./Rdata/coverage_samples.Rdata")
matching.table <- read.table("./sample_coverage/sample_table.txt", header = TRUE, stringsAsFactors = FALSE)

#Function to replicate perl script Coverage_processor.pl length - mRNA length cut-off (includes 5`-UTR), save - should the file be saved
#accepts coverage_samples Rdata as the input. length - mRNA length cut-off (includes 5`-UTR), save - should the file be saved
coverage_processor <- function(input, length = 2000, save = FALSE, report = FALSE) {
  coverage_filtered <- lapply(input, function(x) { lapply(x, function(y)  {   if(length(y) >= length) {return(y[1:length])}})}) #leaves only mRNA longer than specified by length
  coverage_filtered <- lapply(coverage_filtered, function(x) {x[!sapply(x, is.null)] })
  coverage_normalized <- lapply(coverage_filtered, function(x) {
    lapply(x, function(y) {
      average <- mean(y)  
      if(average >= 0.25){  return(y/average)   }
      else {return(NULL)}
    })
  }); rm(coverage_filtered)
  coverage_normalized <- lapply(coverage_normalized, function(x) {x[!sapply(x, is.null)] })
  coverage_matrix     <- lapply(coverage_normalized, function(x) { matrix(unlist(x), byrow=TRUE, nrow=length(x))  })
  densities_sample    <- lapply(coverage_matrix, function(x) {colSums(x) / nrow(x)})
  if(isTRUE(save)) { save(densities_sample, file = "./Rdata/density_samples.Rdata")  }
  if(isTRUE(report)) { print(paste0("processing is completed for cut-off ",length))  }
  return(densities_sample)
} #example: coverage_processor(input = coverage_samples, length = 2000, save = FALSE, report = TRUE)

#function to retrieve translation rate (no statistics) and intermediate analysis, Warning: function does not adjust all parameters by itself, see the spline fitting that requires manual selection of df parameter
#accepts output from coverage_processor() as input. 
translation_rate <- function(input, match, simplify = FALSE, method = "mean", report = NA) {
    #method can be "derivative" or "mean". It reflects how the location of the center of the ribosome coverage curve is determined.
    data0 <- vector("list", length = length(unique(match[,3])));  names(data0) <- sort(unique(match[,3]))
    for(i in 1:length(data0)) {  
      data0[[i]] <- vector("list", length = length(unique(match[,2]))); names(data0[[i]]) <- as.character(sort(unique(match[,2])))
      for(j in 1:length(data0[[i]])) {
        data0[[i]][[j]] <- vector("list", length = length(match[(match[3] == names(data0)[i]) & (match[2] == names(data0[[i]])[j]), 2]))
        names(data0[[i]][[j]]) <- match[(match[3] == names(data0)[i]) & (match[2] == names(data0[[i]])[j]), 1]
        for(k in 1:length(data0[[i]][[j]])) {    data0[[i]][[j]][k] <- input[names(data0[[i]][[j]])[k]]   }
      }
    }
      
   if(isTRUE(simplify)) {
        library(matrixStats)
        #this route simplifies data by starting with the mean coverage track between replicates
        #trim 100 nt of 5`-UTR` and 20 nt of ORF to get rid of the peak at the start codon
        data1 <- lapply(data0, function(a) { lapply(a, function(b) { tail(rowMeans(simplify2array(b)), -120)  })  })
        #smoothen coverage tracks with splines
        data2 <- lapply(data1, function(a) {
          temp <- vector("list", length = length(a)); names(temp) <- names(a)
          for(name in c("0","30","45")) {  temp[[name]] <- smooth.spline(a[[name]], df=12)      }
          temp[["15"]] <- smooth.spline(a[["15"]], df=20) 
          return(temp)
        })
        #calculate the first derivative
        data3 <- lapply(data2, function(a) { lapply(a, function(b) { predict(b, deriv = 1)    }) })
        #extract peaks of the first derivative
        data4 <- lapply(data3, function(a) { lapply(a, function(b) { d <- approx(x=b$y, y=b$x, xout=max(b$y)); return(d$y)  }) })
        #linear model with means only
        data5 <- sapply(data4, function(a) {
          temp <- data.frame("distance" = unlist(a[2:4], use.names = FALSE)/3, "time"=as.numeric(names(a)[2:4]))
          model <- lm(data=temp, formula = distance ~ time)
          return(model$coefficients[[2]])
        })
        output <- vector("list", length = 6);
        output[["raw"]]        <- data0
        output[["trim"]]       <- data1
        output[["spline"]]     <- data2
        output[["derivative"]] <- data3
        output[["peaks"]]      <- data4
        output[["rate"]]       <- data5
   }
   else {
        #this route preserves individuality of replicates for subsequent statistical analysis
        #trim 100 nt of 5`-UTR` and 20 nt of ORF to get rid of the peak at the start codon
        data1 <- lapply(data0, function(a) { lapply(a, function(b) { lapply(b, function(c) {   tail(c, -120)   }) }) })
        #smoothen coverage tracks with splines
        data2 <- lapply(data1, function(a) {
          temp <- vector("list", length = length(a)); names(temp) <- names(a)
          for(name in c("0","30","45")) {  temp[[name]] <- lapply(a[[name]], function(b) {  smooth.spline(b, df=12) })     }
          temp[["15"]] <- lapply(a[["15"]], function(b) {  smooth.spline(b, df=20) })
          return(temp)
        })
        
        #first approach - finding the mean of the coverage slopes (replicates are still individual)
        if(method == "mean") {
                    data3 <- lapply(data2, function(a) { lapply(a[c('15','30','45')], function(b) { lapply(b, function(c) {
                      upper <- which.max(c$y[100:length(c$y)]) #100 nucleotides is an extra trim in addition to 20nt that were already trimmed
                      lower <- which.min(head(c$y, n = upper + 100))
                      interval <- c$y[lower:(upper + 100)]
                      position_of_mean <- approx(x=interval,y=c$x[1:length(interval)], xout=((max(interval, na.rm = TRUE) - min(interval, na.rm=TRUE)) / 2) + min(interval, na.rm=TRUE))
                      output <- position_of_mean$y + lower
                      return(output)
                    }) }) })
                    
                    #linear models with means only
                    data4 <- lapply(data3, function(a) { lapply(a, function(b) {  mean(unlist(b))   })  })
                    data5 <- lapply(data4, function(a) {
                      means <- data.frame("distance" = unlist(a, use.names = FALSE)/3, "time"=as.numeric(names(a)))
                      model <- lm(data=means, formula = distance ~ time)
                      output <- model
                      return(output)
                    })
                    
                    #full linear model with replicates
                    library(reshape2)
                    data6 <- lapply(data3, function(a) {
                      max_replicates_per_timepoint <- max(sapply(a, function(b) length(b)))
                      df <- vector("list", length = length(names(a))); names(df) <- names(a);
                      for(i in names(df)) { 
                        temp <- unlist(a[[i]], use.names = FALSE); length(temp) <- max_replicates_per_timepoint
                        df[[i]] <- temp
                      }
                      df <- data.frame(df, check.names = FALSE)
                      df <- melt(df, value.name="distance", id.var=c())
                      df <- df[!is.na(df$distance),]
                      df$distance <- df$distance / 3
                      row.names(df) <- NULL
                      colnames(df)[1] <- "time"
                      df$time <- as.numeric(as.character(df$time))
                      model <- lm(data=df, formula = distance ~ time)
                      return(list(model=model, data=df))
                    })
                    output <- vector("list");
                    output[["raw"]]              <- data0
                    output[["trim"]]             <- data1
                    output[["spline"]]           <- data2
                    output[["positions"]]        <- data3
                    output[["mean_positions"]]   <- data4
                    output[["mean_lm"]]          <- data5
                    output[["full_lm"]]          <- data6
        }

        #second approach - derivative
        else if(method == "derivative") { 
              #calculate the first derivative
              data3 <- lapply(data2, function(a) { lapply(a, function(b) { lapply(b, function(c) {  predict(c, deriv = 1)  })  }) })
              #extract peaks of the first derivative
              data4 <- lapply(data3, function(a) { lapply(a, function(b) { lapply(b, function(c) {  d <- approx(x=c$y, y=c$x, xout=max(c$y)); return(d$y) }) }) })
              #calculate means across peaks
              data5 <- lapply(data4, function(a) { lapply(a, function(b) {  mean(unlist(b))   })  })
              #linear model with means only
              data6 <- sapply(data5, function(a) {
                temp <- data.frame("distance" = unlist(a[2:4], use.names = FALSE)/3, "time"=as.numeric(names(a)[2:4]))
                model <- lm(data=temp, formula = distance ~ time)
                return(model)
              })
              #proper linear model with replicates
              library(reshape2)
              data7 <- lapply(data4, function(a) {
                length <- max(sapply(a, function(b) length(b)))
                df <- vector("list", length = length(names(a))); names(df) <- names(a);
                for(i in names(df)) { 
                  temp <- unlist(a[[i]], use.names = FALSE); length(temp) <- length
                  df[[i]] <- temp
                }
                df <- data.frame(df, check.names = FALSE)
                df <- melt(df[c(-1)], value.name="distance", id.var=c())
                df <- df[!is.na(df$distance),]
                df$distance <- df$distance / 3
                row.names(df) <- NULL
                colnames(df)[1] <- "time"
                df$time <- as.numeric(as.character(df$time))
                model <- lm(data=df, formula = distance ~ time)
                return(list(model=model, data=df))
              })
              output <- vector("list");
              output[["raw"]]        <- data0
              output[["trim"]]       <- data1
              output[["spline"]]     <- data2
              output[["derivative"]] <- data3
              output[["peaks"]]      <- data4
              output[["mean_peaks"]] <- data5
              output[["mean_lm"]]    <- data6
              output[["full_lm"]]    <- data7
        }
        else { stop(paste0("unsupported method ", method))}
        

   }
  if(is.na(report)) {return(output)}
  else{return(output[c(report)])}
}# example: translation_rate(input = density_samples, match = matching.table, method = "mean")

#function to draw a nice ggplot for coverage tracks with means and sd
#Accepts output of translation_rate() is the input. Order - the desired order of tissues on the plot
rates_plot <- function(input, filename="test.pdf", order = NULL) {
  
  library(matrixStats)
  #mean coverage and standard deviation
  mean <- lapply(input[["raw"]], function(x) { lapply(x, function(y) { rowMeans(simplify2array(y))[101:2000]  })  })
  stdev <- lapply(input[["raw"]], function(x) { lapply(x, function(y) { rowSds(simplify2array(y))[101:2000]  })  })
  
  if(!is.null(order)) { mean <- mean[order]; stdev <- stdev[order]; }

  #smoothen +/- SD for aesthetic purposes
  mean_plus_stdev <- relist(unlist(mean) + unlist(stdev), skeleton = mean)
  mean_minus_stdev <- relist(unlist(mean) - unlist(stdev), skeleton = mean)
  
  upper_loess <- lapply(mean_plus_stdev,  function(x) { lapply(x, function(y) { loess(y ~c(1:length(y)), span = 0.05)  })  })
  lower_loess <- lapply(mean_minus_stdev, function(x) { lapply(x, function(y) { loess(y ~c(1:length(y)), span = 0.05)  })  })
  
  #Plot
  pdf(file = filename, width = 3, height = length(names(mean)))
  par(mfrow=c(length(names(mean)),1), mar=c(3,5,1,3), mgp=c(3,2,0), cex=0.2, cex.lab=2, cex.axis=2, cex.main=2, lend = 1)
  
  for(i in names(mean)) {
    plot(NULL, xlab = "", ylab= "", type="l", ylim=c(0,3), xlim=c(0,1500))
    if(length(upper_loess[[i]][["0"]]$fitted) > 1) {
      polygon(y=c(upper_loess[[i]][["0"]]$fitted,rev(lower_loess[[i]][["0"]]$fitted)),   x=c(c(1:length(upper_loess[[i]][["0"]]$fitted)),rev(c(1:length(upper_loess[[i]][["0"]]$fitted)))), col="grey75", border=NA)
    }
   
    polygon(y=c(upper_loess[[i]][["45"]]$fitted,rev(lower_loess[[i]][["45"]]$fitted)), x=c(c(1:length(upper_loess[[i]][["45"]]$fitted)),rev(c(1:length(upper_loess[[i]][["45"]]$fitted)))), col="#fcbba1", border=NA)
    polygon(y=c(upper_loess[[i]][["30"]]$fitted,rev(lower_loess[[i]][["30"]]$fitted)), x=c(c(1:length(upper_loess[[i]][["30"]]$fitted)),rev(c(1:length(upper_loess[[i]][["30"]]$fitted)))), col="#addd8e", border=NA)
    polygon(y=c(upper_loess[[i]][["15"]]$fitted,rev(lower_loess[[i]][["15"]]$fitted)), x=c(c(1:length(upper_loess[[i]][["15"]]$fitted)),rev(c(1:length(upper_loess[[i]][["15"]]$fitted)))), col="#9ecae1", border=NA)
    
    if(length(upper_loess[[i]][["0"]]$fitted) > 1) {  lines(mean[[i]][["0"]],  col="grey25")  }
    lines(mean[[i]][["15"]], col="#023858")
    lines(mean[[i]][["30"]], col="#004529")
    lines(mean[[i]][["45"]], col="#cb181d") 
  }
  dev.off()
}

#Visualize translation elongation for mRNAs longer than 2000 (includes 100 nt from 5`-UTR)
density <- coverage_processor(input = coverage_samples, length = 2000, save = FALSE, report = FALSE)

analysis_3mo  <- translation_rate(input = density, match = matching.table[matching.table$age == 3 & matching.table$timepoint < 60,], method = "mean")
analysis_18mo <- translation_rate(input = density, match = matching.table[matching.table$age == 18 & matching.table$timepoint < 60,], method = "mean")

rates_plot(input = analysis_3mo, filename = "coverage_rates_3mo.pdf", order=c("liver","kidney","skeletal"))
rates_plot(input = analysis_18mo, filename = "coverage_rates_18mo.pdf")


      #plot linear models for 3 organs from young mice on the same plot
      cairo_pdf(filename = "rates_3mo_organs.pdf", width = 3, height = 3)
      par(mgp=c(3,1,0), cex=0.6, cex.lab=1, cex.axis=1, cex.main=1)
      plot(NULL, xlab = "", ylab= "", xlim=c(0,60), ylim=c(25,300), las=1)
      clip(14.5,45.5,0,300)
      abline(analysis_3mo[["full_lm"]][["liver"]][["model"]], col="#e31a1c", lwd=2.5)
      abline(analysis_3mo[["full_lm"]][["kidney"]][["model"]], col="#33a02c", lwd=2.5)
      abline(analysis_3mo[["full_lm"]][["skeletal"]][["model"]], col="#1f78b4", lwd=2.5)
      clip(0,60,0,300)
      if(exists("mean_peaks", analysis_3mo)) {
        lines(x=c(15,30,45), y=unlist(analysis_3mo$mean_peaks$liver)[-1]/3, type="p", pch=16, col="#e31a1c")
        lines(x=c(15,30,45), y=unlist(analysis_3mo$mean_peaks$kidney)[-1]/3, type="p", pch=16, col="#33a02c")
        lines(x=c(15,30,45), y=unlist(analysis_3mo$mean_peaks$skeletal)[-1]/3, type="p", pch=16, col="#1f78b4")
      } else {
        lines(x=c(15,30,45), y=unlist(analysis_3mo$mean_positions$liver)/3, type="p", pch=16, col="#e31a1c")
        lines(x=c(15,30,45), y=unlist(analysis_3mo$mean_positions$kidney)/3, type="p", pch=16, col="#33a02c")
        lines(x=c(15,30,45), y=unlist(analysis_3mo$mean_positions$skeletal)/3, type="p", pch=16, col="#1f78b4")
      }
      dev.off()
 
      #plot elongation rate slow-down after 45 seconds
      cairo_pdf(filename = "rates_with_60-90sec.pdf", width = 3, height = 3)
      par(mgp=c(3,1,0), cex=0.6, cex.lab=1, cex.axis=1, cex.main=1)
      plot(NULL, xlab = "", ylab= "", xlim=c(0,100), ylim=c(25,400), las=1)
      abline(analysis[["full_lm"]][["liver"]][["model"]], col="#e31a1c", lwd=2.5)
      abline(analysis[["full_lm"]][["kidney"]][["model"]], col="#33a02c", lwd=2.5)
      abline(analysis[["full_lm"]][["skeletal"]][["model"]], col="#1f78b4", lwd=2.5)
      lines(x=c(15,30,45,60), y=c(unlist(analysis$mean_positions$liver)/3, 326), type="p", pch=16, col="#e31a1c")
      lines(x=c(15,30,45,60), y=c(unlist(analysis$mean_positions$kidney)/3, 260), type="p", pch=16, col="#33a02c")
      lines(x=c(15,30,45,60), y=c(unlist(analysis$mean_positions$skeletal)/3, 215), type="p", pch=16, col="#1f78b4")
      lines(x= 90, y= 340, type="p", pch=16, col="#e31a1c")
      x <- c(15, 30, 45, 60, 90, 120)
      y <- c(unlist(analysis$mean_positions$liver)/3, 326, 340, 340)
      f <- function(x,a,b,c,d) {(a*x^3) + (b*x^2) + (c*x) + d}
      fit <- nls(y ~ f(x,a,b,c,d), start = c(a=1, b=1, c=1, d=1)) 
      co <- coef(fit)
      curve(f(x, a=co[1], b=co[2], c=co[3], d=co[4]), add = TRUE, col="pink", lwd=2.5) 
      dev.off()
      
      
      #  Optional: plot linear model for a single organ, study coefficients
      #  plot(x=analysis[["full_lm"]][["liver"]][["data"]]$time, y= analysis[["full_lm"]][["liver"]][["data"]]$distance, type="p")
      #  abline(analysis[["full_lm"]][["liver"]][["model"]])
      #  summary(analysis[["full_lm"]][["liver"]][["model"]])
      
      #barplot with p-values for comparisons between organ rates
      cairo_pdf(filename="barplot_rates.pdf", height=5, width=5)
      opar <- par()
      par(lwd=2)
      bar_data <- c(analysis_3mo$full_lm$liver$model$coefficients[2], analysis_3mo$full_lm$kidney$model$coefficients[2], analysis_3mo$full_lm$skeletal$model$coefficients[2])
      bar_names <- c("liver","kidney","skeletal")
      barCenters <- barplot(bar_data, names.arg=bar_names, ylim=c(0,10), space=0.5)
      errors <- c(coef(summary(analysis_3mo$full_lm$liver$model))[2,2], coef(summary(analysis_3mo$full_lm$kidney$model))[2,2], coef(summary(analysis_3mo$full_lm$skeletal$model))[2,2])
      barplot(bar_data, names.arg=bar_names, ylim=c(0,10), space=0.5, yaxt="n")
      arrows(barCenters, bar_data-errors, barCenters, bar_data +errors, lwd=2, angle=90, code=3, length=0.1)
      axis(2, las=1)
      par(opar)
      dev.off()
      
      #plot the figure that describes the data analysis for translation rates
      cairo_pdf(filename="analysis_example.pdf", height=3, width=3)
      par(mfrow=c(3,1), mar=c(3,5,1,3), mgp=c(3,2,0), cex=0.2, cex.lab=2, cex.axis=2, cex.main=2, lend = 1)
      
      plot(NULL, xlab = "", ylab= "", ylim=c(0,3), xlim=c(0,1500))
      lines(analysis$trim$kidney[["0"]][[1]], type="l", col="grey25")
      lines(analysis$trim$kidney[["15"]][[1]], type ="l", col="#023858")
      lines(analysis$trim$kidney[["30"]][[1]], type ="l", col="#004529")
      lines(analysis$trim$kidney[["45"]][[1]], type ="l", col="#cb181d")
      
      plot(NULL, xlab = "", ylab= "", ylim=c(0,3), xlim=c(0,1500))
      lines(analysis$spline$kidney[["0"]][[1]], type="l", col="grey25")
      lines(analysis$spline$kidney[["15"]][[1]], type ="l", col="#023858")
      lines(analysis$spline$kidney[["30"]][[1]], type ="l", col="#004529")
      lines(analysis$spline$kidney[["45"]][[1]], type ="l", col="#cb181d")
      
      plot(NULL,xlab = "", ylab= "", ylim=c(-0.006, 0.007), xlim=c(0,1500))
      lines(analysis$derivative$kidney[["0"]][[1]], type="l", col="grey25")
      lines(analysis$derivative$kidney[["15"]][[1]], type ="l", col="#023858")
      lines(analysis$derivative$kidney[["30"]][[1]], type ="l", col="#004529")
      lines(analysis$derivative$kidney[["45"]][[1]], type ="l", col="#cb181d")
      dev.off()

#Compare translation rates between organs, use linear model      
#Accepts output from translation_rate().
prepare_data <- function(input) {
                      df <- vector("list")
                      for(i in names(input$full_lm)) { df[[i]] <- data.frame(cbind(input$full_lm[[i]]$data, "organ" = rep(as.character(i),nrow(input$full_lm[[i]]$data)))) }
                      df <- do.call("rbind", df)
                      df$organ <- as.character(df$organ)
                      row.names(df) <- NULL
                      return(df)
}
                                      
organs_df <- prepare_data(input = analysis_3mo)
organs_df <- rbind(organs_df[organs_df$organ == "liver",], organs_df[organs_df$organ != "liver",])
organs_df$organ <- as.factor(organs_df$organ)
organs_df$organ <- factor(organs_df$organ , levels = c("liver", "kidney", "skeletal"))
row.names(organs_df) <- NULL

#linear model with all organs (better p-values)
full_model <- lm(data=organs_df, formula = distance ~ time*organ)

#pair-wise organ comparison (worse p-values)
l_k_model <- lm(data=organs_df[organs_df$organ != "skeletal",], formula = distance ~ time*organ)
l_s_model <- lm(data=organs_df[organs_df$organ != "kidney",],   formula = distance ~ time*organ)
k_s_model <- lm(data=organs_df[organs_df$organ != "liver",],    formula = distance ~ time*organ)


#Compare translation rates between young and old liver, use linear model      
#Accepts output from translation_rate().
organs_df_3mo  <- prepare_data(input = analysis_3mo)   %>% .[.$organ == 'liver',c('time','distance')] %>% cbind(., data.frame('age'=rep('3', nrow(.))))
organs_df_18mo <- prepare_data(input = analysis_18mo)  %>% .[.$organ == 'liver',c('time','distance')] %>% cbind(., data.frame('age'=rep('18', nrow(.))))
organs_df <- rbind(organs_df_3mo, organs_df_18mo)
organs_df$age <- as.factor(organs_df$age)
organs_df$age <- factor(organs_df$age , levels = c("3", "18"))

full_model <- lm(data=organs_df, formula = distance ~ time*age)



#barplot with p-values for comparisons or rates between ages
cairo_pdf(filename="barplot_liver_rates_ages.pdf", height=5, width=4)
opar <- par()
par(lwd=2)
bar_data <- c(analysis_3mo$full_lm$liver$model$coefficients[2], analysis_18mo$full_lm$liver$model$coefficients[2])
bar_names <- c("3 months","18 months")
barCenters <- barplot(bar_data, names.arg=bar_names, ylim=c(0,8), space=0.5)
errors <- c(coef(summary(analysis_3mo$full_lm$liver$model))[2,2], coef(summary(analysis_18mo$full_lm$liver$model))[2,2])
barplot(bar_data, names.arg=bar_names, ylim=c(0,8), space=0.5, yaxt="n")
arrows(barCenters, bar_data-errors, barCenters, bar_data +errors, lwd=2, angle=90, code=3, length=0.1)
axis(2, las=1)
par(opar)
dev.off()




#translation elongation rate of ribosomal proteins
# extract names of ribosomal proteins
gene_names <- names(coverage_samples$MI26Li.coverage)
names_RP <- gene_names[grep("(^Rpl)|(^Rps)", gene_names)]

#calculate length distribution of RPs (they are fairly short, around 500 nt)
length_RP <- sapply(coverage_samples$MI26Li.coverage[names(coverage_samples$MI26Li.coverage) %in% names_RP], function(x) {length(x)})
plot(density(length_RP))

#extract coverage for RPs only
coverage_RP <- lapply(coverage_samples, function(x) { x[names(x) %in% names_RP] })

#calculate density for RPs
density_RP  <- coverage_processor(input = coverage_RP, length = 500, save = FALSE, report = FALSE)
#calculate density for all proteins with same length settings as RPs
density_all <- coverage_processor(input = coverage_samples, length = 500, save = FALSE, report = FALSE)

analysis_3mo_RP   <- translation_rate(input = density_RP, match = matching.table[matching.table$age == 3 & matching.table$timepoint < 60,], method = "mean")
rates_plot(input = analysis_3mo_RP, filename = "test.pdf", order=c("liver","kidney","skeletal"))












#Plot coverage of harringtonine injected mice
harr <- read.delim("F:/Ribosome_profiling/Mammals/Mouse/Injections/Translation_rate/Harringtonine.txt", sep = "\t", header=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
cairo_pdf(filename = "harringtonine.pdf", height=5, width=4)
par(mfrow=c(ncol(harr),1), mar=c(0.5,1,0,1), oma=c(1,0,0,0), mgp=c(3,0.5,0))
for(i in 1:(ncol(harr)-1)) {
  plot(x=c(-100:(nrow(harr)-101)), y=harr[,i], type="l", ylab="", xaxt="n", yaxt="n")
}
plot(x=c(-100:(nrow(harr)-101)), y=harr[,ncol(harr)], type="l", ylab="", yaxt="n")
dev.off()



######################################################################################
#### Dependence of translation rate estimation on gene length ########################
######################################################################################
load(file="./Rdata/coverage_samples.Rdata")

      # Optional: plot number of mRNAs vs length cut-off (excluding 5'UTR)
      # passed <- c()
      # for(i in 101:3100) {
      #   temp <- lapply(coverage_samples[[1]], function(y)  {   if(length(y) > i) {return(y)}})
      #   temp <- temp[!sapply(temp, is.null)]
      #   passed <- c(passed, length(temp))
      # }
      # cairo_pdf(filename = "mRNA_vs_length_cutoff.pdf", height=4, width=4)
      # par(mar=c(2.5,2.5,1,1), cex.axis=0.75)
      # plot(y=passed, x=c(1:3000), ylim=c(0, max(passed)), xaxt="n",yaxt="n", type="l", lwd=1.5)
      # axis(1, mgp=c(1, 0.3, 0))
      # axis(2, mgp=c(1, 0.75, 0), at=c(0,2000,4000,6000,8000,10000,12000,14000), labels=c(0,2,4,6,8,10,12,14), las=1)
      # mtext(side=1,text="length cut-off [nucleotides]", line=1.25, cex=0.85)
      # mtext(side=2,text="number of mRNAs passed   x1000", line=1.4, cex=0.85)
      # dev.off()

#Dependence of translation rate estimation on gene length cut-off (therefore number of passed genes)
#Computationally exensive, run on a server with multithreading

          #Server code with multithreading 
          # library("parallel")
          # library("pbmcapply")
          # density_samples <- list()
          # density_samples <- pbmclapply(rev(1000:2000), function(i) {
          #   result <- coverage_processor(input = coverage_samples, length = i, save = FALSE, report = FALSE) #copy function code above
          #   return(result)
          # }, mc.cores = 15) #warning: each core allocates ~10 GB of memory
          # save(density_samples, file = "density_samples_bylength.Rdata")
          
load(file = "./Rdata/density_samples_bylength.Rdata")
matching.table <- read.table("./sample_coverage/sample_table.txt", header = TRUE, stringsAsFactors = FALSE)


#Plot dependence of translation rate estimation on CDS length cut-off 
rate_vs_length <- lapply(density_samples[1:1000], function(x){ data <- translation_rate(input = x, match = matching.table, method = "mean", report=c("mean_lm")); return(data[[1]]) })

liver <- sapply(rate_vs_length, function(x) {  return(x$liver$coefficients[2]) })
kidney <- sapply(rate_vs_length, function(x) {  return(x$kidney$coefficients[2]) })
skeletal <- sapply(rate_vs_length, function(x) {  return(x$skeletal$coefficients[2]) })

gene_num <- c()
for(i in 1101:2100) {
  temp <- lapply(coverage_samples[[1]], function(y)  {   if(length(y) > i) {return(y)}})
  temp <- temp[!sapply(temp, is.null)]
  gene_num <- c(gene_num, length(temp))
}
gene_num_scaled <- (7.5 - 3) / (max(gene_num) - min(gene_num)) * (gene_num - min(gene_num)) + 3

cairo_pdf(filename="rate_vs_length.pdf", height=4, width=6)
par(mar=c(3,3,1,3))
plot(NULL, yaxt= "n", xlab = "", ylim=c(3,8), xlim=rev(range(1001,2000)))
polygon(y=c(rev(gene_num_scaled), rep(3,length(gene_num_scaled))), x=c(rev(1001:2000),c(1001:2000)), col="#ffeda0", border=NA)
lines(liver, type="l", x=rev(1001:2000), col="#e31a1c", lwd=2)
lines(kidney, type="l", x=rev(1001:2000), col="#33a02c",lwd=2)
lines(skeletal, type="l", x=rev(1001:2000), col="#1f78b4",lwd=2)
axis(side=2, las = 1)
axis(side=4, at=c(3,4,5,6,7,8), labels=c("3700","5100","6500","7900","9300", "10700"))
dev.off()  



# Optional: function to compare coverage tracks and spline/derivative analysis between two rate_vs_length[[i]] or rate_vs_expression[[i]] (make sure "raw" is reported with translation_rate())
# rates_plot2 <- function(input1, input2, mean.mode = TRUE,
#                         filename1="test_raw.pdf", 
#                         filename2="test_trim.pdf", 
#                         filename3="test_spline.pdf", 
#                         filename4="test_derivative.pdf") {
#   
#                     if(isTRUE(mean.mode)) {
#                           library(matrixStats)
#                           #mean coverage and standard deviation
#                           mean1 <- lapply(input1[["raw"]], function(x) { lapply(x, function(y) { rowMeans(simplify2array(y))[1:2000]  })  })
#                           mean2 <- lapply(input2[["raw"]], function(x) { lapply(x, function(y) { rowMeans(simplify2array(y))[1:2000]  })  })
#           
#                           #Plot1
#                           pdf(file = filename1, width = 3, height = length(names(mean1)))
#                           par(mfrow=c(length(names(mean1)),1), mar=c(3,5,1,3), mgp=c(3,2,0), cex=0.2, cex.lab=2, cex.axis=2, cex.main=2, lend = 1)
#                           for(i in names(mean1)) {
#                             plot(NULL, xlab = "", ylab= "", type="l", ylim=c(0,3), xlim=c(0,1500))
#                             lines(mean1[[i]][["0"]],  col="#023858")
#                             lines(mean1[[i]][["15"]], col="#023858")
#                             lines(mean1[[i]][["30"]], col="#023858")
#                             lines(mean1[[i]][["45"]], col="#023858")
#           
#                             lines(mean2[[i]][["0"]],  col="#cb181d")
#                             lines(mean2[[i]][["15"]], col="#cb181d")
#                             lines(mean2[[i]][["30"]], col="#cb181d")
#                             lines(mean2[[i]][["45"]], col="#cb181d")
#                             
#                             text(x = 750, y = 2.8, labels = i, cex = 2.5)
#                           }; dev.off()
#                           
#                           #trim 100 nt 5'-UTR and 20 nt from CDS start.
#                           mean1 <- lapply(mean1, function(x) { lapply(x, function(y) { return(tail(y, -120))  })  })
#                           mean2 <- lapply(mean2, function(x) { lapply(x, function(y) { return(tail(y, -120))  })  })
#           
#                           #Plot2
#                           pdf(file = filename2, width = 3, height = length(names(mean1)))
#                           par(mfrow=c(length(names(mean1)),1), mar=c(3,5,1,3), mgp=c(3,2,0), cex=0.2, cex.lab=2, cex.axis=2, cex.main=2, lend = 1)
#                           for(i in names(mean1)) {
#                             plot(NULL, xlab = "", ylab= "", type="l", ylim=c(0,3), xlim=c(0,1500))
#                             lines(mean1[[i]][["0"]],  col="#023858")
#                             lines(mean1[[i]][["15"]], col="#023858")
#                             lines(mean1[[i]][["30"]], col="#023858")
#                             lines(mean1[[i]][["45"]], col="#023858")
#           
#                             lines(mean2[[i]][["0"]],  col="#cb181d")
#                             lines(mean2[[i]][["15"]], col="#cb181d")
#                             lines(mean2[[i]][["30"]], col="#cb181d")
#                             lines(mean2[[i]][["45"]], col="#cb181d")
#           
#                             text(x = 750, y = 2.8, labels = i, cex = 2.5)
#                           }; dev.off()
#           
#                           structure <- rep(NA ,length(unlist(mean1)))
#                           spline1 <- relist(structure, skeleton = mean1)
#                           spline2 <- relist(structure, skeleton = mean2)
#           
#                           for(i in names(mean1)) {
#                             spline1[[i]][c(1,3,4)] <- lapply(mean1[[i]][c(1,3,4)], function(x) { smooth.spline(x[!is.na(x)], df=12)    })
#                             spline2[[i]][c(1,3,4)] <- lapply(mean2[[i]][c(1,3,4)], function(x) { smooth.spline(x[!is.na(x)], df=12)    })
#                             spline1[[i]][c(2)] <- lapply(mean1[[i]][c(2)], function(x) { smooth.spline(x[!is.na(x)], df=20)    })
#                             spline2[[i]][c(2)] <- lapply(mean2[[i]][c(2)], function(x) { smooth.spline(x[!is.na(x)], df=20)    })
#                           }
#           
#                           derivative1 <- lapply(spline1, function(x) {  lapply(x, function(y) { predict(y, deriv = 1)     })   })
#                           derivative2 <- lapply(spline2, function(x) {  lapply(x, function(y) { predict(y, deriv = 1)     })   })
#           
#                           #Plot3
#                           pdf(file = filename3, width = 3, height = length(names(mean1)))
#                           par(mfrow=c(length(names(mean1)),1), mar=c(3,5,1,3), mgp=c(3,2,0), cex=0.2, cex.lab=2, cex.axis=2, cex.main=2, lend = 1)
#                           for(i in names(mean1)) {
#                               plot(NULL, xlab = "", ylab= "", type="l", ylim=c(0,3), xlim=c(0,1500))
#                               lines(spline1[[i]][["0"]],  col="#023858")
#                               lines(spline1[[i]][["15"]], col="#023858")
#                               lines(spline1[[i]][["30"]], col="#023858")
#                               lines(spline1[[i]][["45"]], col="#023858")
#           
#                               lines(spline2[[i]][["0"]],  col="#cb181d")
#                               lines(spline2[[i]][["15"]], col="#cb181d")
#                               lines(spline2[[i]][["30"]], col="#cb181d")
#                               lines(spline2[[i]][["45"]], col="#cb181d")
#           
#                               text(x = 750, y = 2.8, labels = i, cex = 2.5)
#                           }; dev.off()
#                           
#                           #Plot4
#                           pdf(file = filename4, width = 3, height = length(names(mean1)))
#                           par(mfrow=c(length(names(mean1)),1), mar=c(3,5,1,3), mgp=c(3,2,0), cex=0.2, cex.lab=2, cex.axis=2, cex.main=2, lend = 1)
#                           for(i in names(mean1)) {
#                             plot(NULL, xlab = "", ylab= "", type="l", ylim=c(-0.006, 0.007), xlim=c(0,1500))
#                             lines(derivative1[[i]][["15"]], col="#023858")
#                             lines(derivative1[[i]][["30"]], col="#023858")
#                             lines(derivative1[[i]][["45"]], col="#023858")
#                             
#                             lines(derivative2[[i]][["15"]], col="#cb181d")
#                             lines(derivative2[[i]][["30"]], col="#cb181d")
#                             lines(derivative2[[i]][["45"]], col="#cb181d")
#                             
#                             text(x = 750, y = 0.006, labels = i, cex = 2.5)
#                           }; dev.off()
#                       }
#                       else{
#                         #all replicates are treated individually
#                         mean1 <- lapply(input1[["raw"]], function(x) { lapply(x, function(y) { lapply(y, function(z) {  z   })  })  })
#                         mean2 <- lapply(input2[["raw"]], function(x) { lapply(x, function(y) { lapply(y, function(z) {  z   })  })  })
#                         #Plot1
#                         pdf(file = filename1, width = 3, height = length(names(mean1)))
#                         par(mfrow=c(length(names(mean1)),1), mar=c(3,5,1,3), mgp=c(3,2,0), cex=0.2, cex.lab=2, cex.axis=2, cex.main=2, lend = 1)
#                         for(i in names(mean1)) {
#                           plot(NULL, xlab = "", ylab= "", type="l", ylim=c(0,3), xlim=c(0,1500))
#                           for(j in 1:length(mean1[[i]][["0"]]))  { lines(unlist(mean1[[i]][["0"]][[j]]),  col="#023858", lwd=0.5)   }
#                           for(j in 1:length(mean1[[i]][["15"]])) { lines(unlist(mean1[[i]][["15"]][[j]]),  col="#023858", lwd=0.5)   }
#                           for(j in 1:length(mean1[[i]][["30"]])) { lines(unlist(mean1[[i]][["30"]][[j]]),  col="#023858", lwd=0.5)   }
#                           for(j in 1:length(mean1[[i]][["45"]])) { lines(unlist(mean1[[i]][["45"]][[j]]),  col="#023858", lwd=0.5)   }
#                           
#                           for(j in 1:length(mean1[[i]][["0"]])) { lines(unlist(mean2[[i]][["0"]][[j]]),  col="#cb181d", lwd=0.5)   }
#                           for(j in 1:length(mean1[[i]][["15"]])) { lines(unlist(mean2[[i]][["15"]][[j]]),  col="#cb181d", lwd=0.5)   }
#                           for(j in 1:length(mean1[[i]][["30"]])) { lines(unlist(mean2[[i]][["30"]][[j]]),  col="#cb181d", lwd=0.5)   }
#                           for(j in 1:length(mean1[[i]][["45"]])) { lines(unlist(mean2[[i]][["45"]][[j]]),  col="#cb181d", lwd=0.5)   }
#                           
#                           text(x = 750, y = 2.8, labels = i, cex = 2.5)
#                         }; dev.off()
#                         
#                         #trim 100 nt 5'-UTR and 20 nt from CDS start.
#                         mean1 <- lapply(mean1, function(x) { lapply(x, function(y) { lapply(y, function(z) {return(tail(z, -120))})  })  })
#                         mean2 <- lapply(mean2, function(x) { lapply(x, function(y) { lapply(y, function(z) {return(tail(z, -120))})  })  })
#                         #Plot2
#                         pdf(file = filename2, width = 3, height = length(names(mean1)))
#                         par(mfrow=c(length(names(mean1)),1), mar=c(3,5,1,3), mgp=c(3,2,0), cex=0.2, cex.lab=2, cex.axis=2, cex.main=2, lend = 1)
#                         for(i in names(mean1)) {
#                           plot(NULL, xlab = "", ylab= "", type="l", ylim=c(0,3), xlim=c(0,1500))
#                           for(j in 1:length(mean1[[i]][["0"]]))  { lines(unlist(mean1[[i]][["0"]][[j]]),  col="#023858", lwd=0.5)   }
#                           for(j in 1:length(mean1[[i]][["15"]])) { lines(unlist(mean1[[i]][["15"]][[j]]),  col="#023858", lwd=0.5)   }
#                           for(j in 1:length(mean1[[i]][["30"]])) { lines(unlist(mean1[[i]][["30"]][[j]]),  col="#023858", lwd=0.5)   }
#                           for(j in 1:length(mean1[[i]][["45"]])) { lines(unlist(mean1[[i]][["45"]][[j]]),  col="#023858", lwd=0.5)   }
#                           
#                           for(j in 1:length(mean2[[i]][["0"]])) { lines(unlist(mean2[[i]][["0"]][[j]]),  col="#cb181d", lwd=0.5)   }
#                           for(j in 1:length(mean2[[i]][["15"]])) { lines(unlist(mean2[[i]][["15"]][[j]]),  col="#cb181d", lwd=0.5)   }
#                           for(j in 1:length(mean2[[i]][["30"]])) { lines(unlist(mean2[[i]][["30"]][[j]]),  col="#cb181d", lwd=0.5)   }
#                           for(j in 1:length(mean2[[i]][["45"]])) { lines(unlist(mean2[[i]][["45"]][[j]]),  col="#cb181d", lwd=0.5)   }
#                           
#                           text(x = 750, y = 2.8, labels = i, cex = 2.5)
#                         }; dev.off()
#                         
#                         structure <- rep(NA ,length(unlist(mean1)))  
#                         spline1 <- relist(structure, skeleton = mean1)
#                         spline2 <- relist(structure, skeleton = mean2)
#                         
#                         for(i in names(mean1)) {
#                           spline1[[i]][c(1,3,4)] <- lapply(mean1[[i]][c(1,3,4)], function(x) { lapply(x, function(y)  { smooth.spline(y[!is.na(y)], df=12)   })   })
#                           spline2[[i]][c(1,3,4)] <- lapply(mean2[[i]][c(1,3,4)], function(x) { lapply(x, function(y)  { smooth.spline(y[!is.na(y)], df=12)   })   })
#                           spline1[[i]][c(2)] <- lapply(mean1[[i]][c(2)], function(x) { lapply(x, function(y)  { smooth.spline(y[!is.na(y)], df=20)   })   })
#                           spline2[[i]][c(2)] <- lapply(mean2[[i]][c(2)], function(x) { lapply(x, function(y)  { smooth.spline(y[!is.na(y)], df=20)   })   })
#                         }
#                         
#                         derivative1 <- lapply(spline1, function(x) {  lapply(x, function(y) { lapply(y, function(z) { predict(z, deriv = 1)  })    })   })
#                         derivative2 <- lapply(spline2, function(x) {  lapply(x, function(y) { lapply(y, function(z) { predict(z, deriv = 1)  })    })   })
#                         
#                         #Plot3
#                         pdf(file = filename3, width = 3, height = length(names(mean1)))
#                         par(mfrow=c(length(names(mean1)),1), mar=c(3,5,1,3), mgp=c(3,2,0), cex=0.2, cex.lab=2, cex.axis=2, cex.main=2, lend = 1)
#                         for(i in names(mean1)) {
#                           plot(NULL, xlab = "", ylab= "", type="l", ylim=c(0,3), xlim=c(0,1500))
#                           for(j in 1:length(spline1[[i]][["0"]]))  { lines(spline1[[i]][["0"]][[j]],  col="#023858", lwd=0.5)    }
#                           for(j in 1:length(spline1[[i]][["15"]])) { lines(spline1[[i]][["15"]][[j]],  col="#023858", lwd=0.5)   }
#                           for(j in 1:length(spline1[[i]][["30"]])) { lines(spline1[[i]][["30"]][[j]],  col="#023858", lwd=0.5)   }
#                           for(j in 1:length(spline1[[i]][["45"]])) { lines(spline1[[i]][["45"]][[j]],  col="#023858", lwd=0.5)   }
#                           
#                           for(j in 1:length(spline2[[i]][["0"]]))  { lines(spline2[[i]][["0"]][[j]],  col="#cb181d", lwd=0.5)    }
#                           for(j in 1:length(spline2[[i]][["15"]])) { lines(spline2[[i]][["15"]][[j]],  col="#cb181d", lwd=0.5)   }
#                           for(j in 1:length(spline2[[i]][["30"]])) { lines(spline2[[i]][["30"]][[j]],  col="#cb181d", lwd=0.5)   }
#                           for(j in 1:length(spline2[[i]][["45"]])) { lines(spline2[[i]][["45"]][[j]],  col="#cb181d", lwd=0.5)   }
#                           
#                           text(x = 750, y = 2.8, labels = i, cex = 2.5)
#                         }; dev.off()
#                         
#                         #Plot4
#                         pdf(file = filename4, width = 3, height = length(names(mean1)))
#                         par(mfrow=c(length(names(mean1)),1), mar=c(3,5,1,3), mgp=c(3,2,0), cex=0.2, cex.lab=2, cex.axis=2, cex.main=2, lend = 1)
#                         for(i in names(mean1)) {
#                           plot(NULL, xlab = "", ylab= "", type="l", ylim=c(-0.006, 0.007), xlim=c(0,1500))
#                           for(j in 1:length(derivative1[[i]][["15"]])) { lines(derivative1[[i]][["15"]][[j]],  col="#023858", lwd=0.5)   }
#                           for(j in 1:length(derivative1[[i]][["30"]])) { lines(derivative1[[i]][["30"]][[j]],  col="#023858", lwd=0.5)   }
#                           for(j in 1:length(derivative1[[i]][["45"]])) { lines(derivative1[[i]][["45"]][[j]],  col="#023858", lwd=0.5)   }
#       
#                           for(j in 1:length(derivative2[[i]][["15"]])) { lines(derivative2[[i]][["15"]][[j]],  col="#cb181d", lwd=0.5)   }
#                           for(j in 1:length(derivative2[[i]][["30"]])) { lines(derivative2[[i]][["30"]][[j]],  col="#cb181d", lwd=0.5)   }
#                           for(j in 1:length(derivative2[[i]][["45"]])) { lines(derivative2[[i]][["45"]][[j]],  col="#cb181d", lwd=0.5)   }
#       
#                           text(x = 750, y = 0.006, labels = i, cex = 2.5)
#                         }; dev.off()
#                     }
# 
# } #rates_plot2(input1 = rate_vs_expression2[[6484]], input2 = rate_vs_expression2[[6485]], mean.mode = FALSE)





######################################################################################
#### Dependence of translation rate estimation on gene expression level ############## 
######################################################################################
 
#Dependence of translation rate estimation on gene expression
#Computationally exensive, run on a server with multithreading
      
      # #Server code with multithreading
      # setwd("/Red6TB/Maxim/Rtest/")
      # load(file="coverage_samples.Rdata")
      # untreated_organs <- list("liver" = coverage_samples[c("MI26Li.coverage", "MI27Li.coverage", "MI28Li.coverage", "MI29Li.coverage")], 
      #                          "kidney" = coverage_samples[c("MI26K.coverage", "MI27K.coverage")], 
      #                          "skeletal" = coverage_samples[c("MI26SKM.coverage", "MI27SKM.coverage")])
      # 
      # library("data.table")
      # ranked_expression <- lapply(untreated_organs, function(x) { 
      #   data <- lapply(x, function(y) {  lapply(y, function(z) { mean(z) })  })
      #   data <- rbindlist(data)
      #   data <- apply(data, 2, mean)
      #   data <- sort(data)
      #   data <- as.list(data)
      #   return(data)
      # })
      # library("parallel")
      # library("pbmcapply")
      # density_biglist <- pbmclapply(1:13000, function(i) {
      # genes_kept <- lapply(ranked_expression, function(x) {   tail(unlist(names(x)), n= -i)    })
      # coverage_kept <- vector("list", length=length(coverage_samples)); names(coverage_kept) <- names(coverage_samples)
      # for(j in names(coverage_samples)){
      #   organ <- matching.table[matching.table$filename == j, 3]
      #   genes <- coverage_samples[[j]][genes_kept[[organ]]] 
      #   coverage_kept[[j]] <- genes
      # }
      # result <- coverage_processor(input = coverage_kept, length = 1200, save = FALSE, report = FALSE)
      # return(result)
      # }, mc.cores = 10)

load(file = "./Rdata/density_samples_expr.Rdata")
matching.table <- read.table("./sample_coverage/sample_table.txt", header = TRUE, stringsAsFactors = FALSE)

#Plot dependence of translation rate estimation on gene expression 
rate_vs_expression <- lapply(density_samples_expr, function(x){ data <- translation_rate(input = x, match = matching.table, method = "mean", report=c("mean_lm"))
                                                                        output <- vector("list")
                                                                        output[["liver"]] <- data[[1]]$liver$coefficients[[2]]
                                                                        output[["kidney"]] <- data[[1]]$kidney$coefficients[[2]]
                                                                        output[["skeletal"]] <- data[[1]]$skeletal$coefficients[[2]]
                                                                        return(output) })

cairo_pdf(filename="rate_vs_expression.pdf", height=4, width=6)
par(mar=c(3,3,1,3))
plot(NULL, yaxt= "n", xlab = "", ylim=c(3,8), xlim=c(1,length(sapply(rate_vs_expression, function(x) x$liver    ))))
lines(sapply(rate_vs_expression, function(x) x$liver    ), type="l", x=c(1:length(sapply(rate_vs_expression, function(x) x$liver    ))), col="#e31a1c", lwd=2)
lines(sapply(rate_vs_expression, function(x) x$kidney   ), type="l", x=c(1:length(sapply(rate_vs_expression, function(x) x$kidney   ))), col="#33a02c",lwd=2)
lines(sapply(rate_vs_expression, function(x) x$skeletal ), type="l", x=c(1:length(sapply(rate_vs_expression, function(x) x$skeletal ))), col="#1f78b4",lwd=2)
axis(side=2, las = 1)
dev.off()  


#########################################################################################################################################################
################################# Study translation rates of individual genes ###########################################################################
#########################################################################################################################################################

load(file="./Rdata/coverage_samples2.Rdata")
matching.table <- read.table("./sample_coverage/sample_table2.txt", header = TRUE, stringsAsFactors = FALSE)

# extract ribo-seq coverage for control, untreated organs
untreated_organs <- list("liver" = coverage_samples[c("MI26Li.coverage", "MI27Li.coverage", "MI28Li.coverage", "MI29Li.coverage")], 
                         "kidney" = coverage_samples[c("MI26K.coverage", "MI27K.coverage")], 
                         "skeletal" = coverage_samples[c("MI26SKM.coverage", "MI27SKM.coverage")])

#Rank genes based on the expression
library("data.table")
ranked_expression <- lapply(untreated_organs, function(x) { 
  data <- lapply(x, function(y) {  lapply(y, function(z) { mean(z) })  })
  data <- rbindlist(data)
  data <- apply(data, 2, mean)
  data <- sort(data)
  data <- as.list(data)
  return(data)
})
full_gene_set <- names(ranked_expression[[1]])

#pick 500 top expressed genes
liver_ranked    <-  tail(unlist(ranked_expression$liver), n = 500)
kidney_ranked   <-  tail(unlist(ranked_expression$kidney), n = 500)
skeletal_ranked <-  tail(unlist(ranked_expression$skeletal), n = 500)

#select genes that are highly expressed in all three organs
overlap <- liver_ranked[names(liver_ranked) %in% names(kidney_ranked)]
overlap <- overlap[names(overlap) %in% names(skeletal_ranked)]
overlap <- names(overlap)

#Function to access a coverage info for a single gene in all organs (by name)
gene_coverage <- function(gene_name, data, match) {
      #prepare data structure,fill with data
      structure <- vector("list", length = length(unique(match[,3])));  names(structure) <- sort(unique(match[,3]))
      for(i in names(structure)) { 
            structure[[i]] <- vector("list", length = length(unique(match[,2]))) 
            names(structure[[i]]) <- as.character(sort(unique(match[,2]))) 
            for(j in names(structure[[i]])) {    
                  structure[[i]][[j]] <-  vector("list", length = nrow(match[(match[3] == i) & (match[2] == j),]))
                  names(structure[[i]][[j]]) <- match[(match[3] == i) & (match[2] == j), 1]
                  #fill structure with data
                  for(k in names(structure[[i]][[j]])) {   
                        structure[[i]][[j]][[k]] <- data[[k]][[gene_name]]
                  }
            }
      }
      #normalize replicates by total number of reads
      data1 <- lapply(structure, function(a) {lapply(a, function(b) { lapply(b, function(c) { output <- c / sum(c, na.rm = TRUE); return(output)  })})})
      
      #take average across replicate
      library(matrixStats)
      data2 <- lapply(data1, function(a) { lapply(a, function(b) { if(length(b) > 1) {rowMeans(simplify2array(b))} else {return(as.numeric(unlist(b)))}  })  })
      return(data2)
}# Example: gene_coverage(gene_name = "Apoe", data = coverage_samles, match = matching.table)


#Interactive plotting of gene coverage
library(shiny)
library(plotly)
library(DT)
library(webshot)

    server <- function(input, output) {
      #initialize observers
      gene_name <- reactiveValues(value=NULL)
      data_set <- reactiveValues(value=NULL)
      
      observeEvent(input$data_set, {  data_set$value  <- input$data_set
                                      print(data_set$value)
      })
      
      observeEvent(input$gene_table_rows_selected, {  gene_name$value  <- get(data_set$value)[input$gene_table_rows_selected]  })
      
      #prepare data for plotting
      plot_data <- eventReactive (gene_name$value, {
                   query_gene <- gene_coverage(gene_name = gene_name$value, data = coverage_samples, match = matching.table)
                   output <- list("liver" = data.frame(cbind("position" = c(1:length(query_gene[["liver"]][[1]])), as.data.frame(query_gene[["liver"]], col.names=names(query_gene[["liver"]])))),
                                  "kidney" = data.frame(cbind("position" = c(1:length(query_gene[["kidney"]][[1]])), as.data.frame(query_gene[["kidney"]], col.names=names(query_gene[["kidney"]])))),
                                  "skeletal" = data.frame(cbind("position" = c(1:length(query_gene[["skeletal"]][[1]])), as.data.frame(query_gene[["skeletal"]], col.names=names(query_gene[["skeletal"]])))) )
                   return(output)

      })
      
      output$gene_table <- renderDataTable({  datatable( data.frame("gene_table" = get(data_set$value)), 
                                                         caption = 'Table 1: Top common highly expressed genes in three organs.',
                                                         colnames = c("Gene name"),
                                                         filter = "top",
                                                         escape = FALSE,
                                                         options = list(dom = 'ftip',
                                                                        pageLength = 10,
                                                                        columnDefs = list(  list(className = 'dt-center', targets = 1),
                                                                                            list(targets = 1, render = JS(  "function(data, type, row, meta) {",
                                                                                                                            "return type === 'display' && data.length > 30 ?",
                                                                                                                            "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
                                                                                                                            "}"
                                                                                                                          )
                                                                                            )
                                                                        )
                                                         ),
                                                         selection = list(mode = 'single', target = 'row', selected = c(1))
      )})
      
      plotInput <- reactive({
        p1 <- plot_ly(plot_data()[["liver"]], hoverinfo='none') %>%
             add_trace(x = ~position, y = ~X0, name = 'coverage', type = "bar", marker = list(color="#9ecae1"), showlegend = FALSE) %>%
             layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                     xaxis = list(range=c(-nrow(plot_data()[["liver"]])*0.05, nrow(plot_data()[["liver"]])),
                                  title = 'nucleotide position',
                                  autorange=FALSE,
                                  ticks='outside',
                                  tickwidth=2),
                 bargap = 0,
                 barmode = "overlay",
                 plot_bgcolor = 'rgba(235,235,235,1)',
                 legend = list(orientation = "h", bgcolor = "#08519c"),
                 margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p2 <- plot_ly(plot_data()[["liver"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X15, name = 'coverage', type = "bar", marker = list(color="#6baed6"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data()[["liver"]])*0.05, nrow(plot_data()[["liver"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p3 <- plot_ly(plot_data()[["liver"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X30, name = 'coverage', type = "bar", marker = list(color="#4292c6"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data()[["liver"]])*0.05, nrow(plot_data()[["liver"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p4 <- plot_ly(plot_data()[["liver"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X45, name = 'coverage', type = "bar", marker = list(color="#2171b5"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data()[["liver"]])*0.05, nrow(plot_data()[["liver"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p5 <- plot_ly(plot_data()[["liver"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X300, name = 'coverage', type = "bar", marker = list(color="#084594"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data()[["liver"]])*0.05, nrow(plot_data()[["liver"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 40)
          )
        p6 <- plot_ly(plot_data()[["kidney"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X0, name = 'coverage', type = "bar", marker = list(color="#fc9272"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data()[["liver"]])*0.05, nrow(plot_data()[["kidney"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p7 <- plot_ly(plot_data()[["kidney"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X15, name = 'coverage', type = "bar", marker = list(color="#fb6a4a"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data()[["liver"]])*0.05, nrow(plot_data()[["kidney"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p8 <- plot_ly(plot_data()[["kidney"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X30, name = 'coverage', type = "bar", marker = list(color="#ef3b2c"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data()[["liver"]])*0.05, nrow(plot_data()[["kidney"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p9 <- plot_ly(plot_data()[["kidney"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X45, name = 'coverage', type = "bar", marker = list(color="#cb181d"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data()[["liver"]])*0.05, nrow(plot_data()[["kidney"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p10 <- plot_ly(plot_data()[["kidney"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X300, name = 'coverage', type = "bar", marker = list(color="#99000d"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data()[["liver"]])*0.05, nrow(plot_data()[["kidney"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 40)
          )
        p11 <- plot_ly(plot_data()[["skeletal"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X0, name = 'coverage', type = "bar", marker = list(color="#a1d99b"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data()[["liver"]])*0.05, nrow(plot_data()[["skeletal"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p12 <- plot_ly(plot_data()[["skeletal"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X15, name = 'coverage', type = "bar", marker = list(color="#74c476"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data()[["liver"]])*0.05, nrow(plot_data()[["skeletal"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p13 <- plot_ly(plot_data()[["skeletal"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X30, name = 'coverage', type = "bar", marker = list(color="#41ab5d"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data()[["liver"]])*0.05, nrow(plot_data()[["skeletal"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p14 <- plot_ly(plot_data()[["skeletal"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X45, name = 'coverage', type = "bar", marker = list(color="#238b45"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data()[["liver"]])*0.05, nrow(plot_data()[["skeletal"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p15 <- plot_ly(plot_data()[["skeletal"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X300, name = 'coverage', type = "bar", marker = list(color="#005a32"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data()[["liver"]])*0.05, nrow(plot_data()[["skeletal"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 40)
          )

        p <- subplot(p1, p6, p11, p2, p7, p12, p3, p8, p13, p4, p9, p14, p5, p10, p15, nrows = 5, shareY = FALSE, shareX = TRUE, margin = c(0, 0.01, 0, 0.02))
      })
      
    output$coverage_plot <- renderPlotly({ plotInput()  })
    
    #export image locally
    output$imagesave1 <- downloadHandler(
                            filename = function() { paste0(gene_name$value, ".pdf")},
                            content  = function(file) { export(plotInput(), file = file)  })
    
    output$imagesave2 <- downloadHandler(
                            filename = function() { paste0(gene_name$value, ".png")},
                            content  = function(file) { export(plotInput(), file = file, vheight = 400, vwidth = 800, zoom = 2)  })


    }
    ui <- fluidPage(   fluidRow(div(column(width=12, plotlyOutput("coverage_plot", height = 400))), style = "padding:25px;"),
                       fluidRow(column(width=4, tags$head(tags$style("#gene_table  {white-space: nowrap;  }")), div(dataTableOutput("gene_table"), style = "font-size: 80%")),
                                column(width=4, selectInput(inputId = "data_set", label = "Select the data set:", choices = c("overlap", "full_gene_set"), selected = "overlap"),
                                                downloadButton("imagesave1", "Save the image as pdf", style = "background-color:#084594; color:white"),
                                                downloadButton("imagesave2", "Save the image as png", style = "background-color:#005a32; color:white")))
                       
          )

shinyApp(ui = ui, server = server)




##### Atf-4 gene uORFs study case #####
Atf4_sequence <- strsplit(read.csv("Atf4.txt", header = FALSE, stringsAsFactors = FALSE)[,1], split = "")
Atf4_length <- length(Atf4_sequence[[1]])

    library("R.utils")
    coverage_samples <- vector("list", length=length(list.files(path="./Atf4", pattern='coverage')))
    names(coverage_samples) <- list.files(path="./Atf4", pattern='coverage')
    for(i in names(coverage_samples)) {
      file <- read.csv(file = paste0(getwd(), "/Atf4/", i), header = FALSE)
      file <- as.numeric(file[,1])
      if(length(file) < Atf4_length) {
         file <- c(file, rep(0, Atf4_length - length(file)))
      }
      coverage_samples[[i]] <- file
    }

df <- data.frame(coverage_samples)
matching.table <- read.table("./Atf4/sample_table.txt", header = TRUE, stringsAsFactors = FALSE)

Atf4_coverage <- function(data = df, match = matching.table) {
  #prepare data structure,fill with data
  structure <- vector("list", length = length(unique(match[,3])));  names(structure) <- sort(unique(match[,3]))
  for(i in names(structure)) { 
    structure[[i]] <- vector("list", length = length(unique(match[,2]))) 
    names(structure[[i]]) <- as.character(sort(unique(match[,2]))) 
    for(j in names(structure[[i]])) {    
      structure[[i]][[j]] <-  vector("list", length = nrow(match[(match[3] == i) & (match[2] == j),]))
      names(structure[[i]][[j]]) <- match[(match[3] == i) & (match[2] == j), 1]
      #fill structure with data
      for(k in names(structure[[i]][[j]])) {  structure[[i]][[j]][[k]] <- data[[k]] }
    }
  }
  #take sum across replicate
  library(matrixStats)
  output <- lapply(structure, function(a) { lapply(a, function(b) { if(length(b) > 1) {rowSums(simplify2array(b))} else {return(as.numeric(unlist(b)))}  })  })
  return(output)
}

Atf4 <- Atf4_coverage()

server <- function(input, output) {

  #prepare data for plotting
  data <- list("liver"     = data.frame(cbind("position" = c(1:length(Atf4[["liver"]][[1]])), as.data.frame(Atf4[["liver"]][1:4], col.names=names(Atf4[["liver"]][1:4])))),
                "kidney"   = data.frame(cbind("position" = c(1:length(Atf4[["kidney"]][[1]])), as.data.frame(Atf4[["kidney"]][1:4], col.names=names(Atf4[["kidney"]][1:4])))),
                "skeletal" = data.frame(cbind("position" = c(1:length(Atf4[["skeletal"]][[1]])), as.data.frame(Atf4[["skeletal"]][1:4], col.names=names(Atf4[["skeletal"]][1:4])))) )

    
  
  plotInput <- reactive({
    p1 <- plot_ly(data[["liver"]], hoverinfo='none') %>%
      add_trace(x = ~position, y = ~X0, name = 'coverage', type = "bar", marker = list(color="#2171b5"), showlegend = FALSE) %>%
      layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
              xaxis = list(range=c(-nrow(data[["liver"]])*0.05, nrow(data[["liver"]])),
                           title = 'nucleotide position',
                           autorange=FALSE,
                           ticks='outside',
                           tickwidth=2),
              bargap = 0.1,
              barmode = "overlay",
              plot_bgcolor = 'rgba(235,235,235,1)',
              legend = list(orientation = "h", bgcolor = "#08519c")
      )
    p2 <- plot_ly(data[["liver"]], hoverinfo='none') %>%
      add_trace(x = ~position, y = ~X15, name = 'coverage', type = "bar", marker = list(color="#2171b5"), showlegend = FALSE) %>%
      layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
              xaxis = list(range=c(-nrow(data[["liver"]])*0.05, nrow(data[["liver"]])),
                           title = 'nucleotide position',
                           autorange=FALSE,
                           ticks='outside',
                           tickwidth=2),
              bargap = 0.1,
              barmode = "overlay",
              plot_bgcolor = 'rgba(235,235,235,1)',
              legend = list(orientation = "h", bgcolor = "#08519c")
      )
    p3 <- plot_ly(data[["liver"]], hoverinfo='none') %>%
      add_trace(x = ~position, y = ~X30, name = 'coverage', type = "bar", marker = list(color="#2171b5"), showlegend = FALSE) %>%
      layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
              xaxis = list(range=c(-nrow(data[["liver"]])*0.05, nrow(data[["liver"]])),
                           title = 'nucleotide position',
                           autorange=FALSE,
                           ticks='outside',
                           tickwidth=2),
              bargap = 0.1,
              barmode = "overlay",
              plot_bgcolor = 'rgba(235,235,235,1)',
              legend = list(orientation = "h", bgcolor = "#08519c")
      )
    p4 <- plot_ly(data[["liver"]], hoverinfo='none') %>%
      add_trace(x = ~position, y = ~X45, name = 'coverage', type = "bar", marker = list(color="#2171b5"), showlegend = FALSE) %>%
      layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
              xaxis = list(range=c(-nrow(data[["liver"]])*0.05, nrow(data[["liver"]])),
                           title = 'nucleotide position',
                           autorange=FALSE,
                           ticks='outside',
                           tickwidth=2),
              bargap = 0.1,
              barmode = "overlay",
              plot_bgcolor = 'rgba(235,235,235,1)',
              legend = list(orientation = "h", bgcolor = "#08519c")
      )
    p5 <- plot_ly(data[["kidney"]], hoverinfo='none') %>%
      add_trace(x = ~position, y = ~X0, name = 'coverage', type = "bar", marker = list(color="#cb181d"), showlegend = FALSE) %>%
      layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
              xaxis = list(range=c(-nrow(data[["liver"]])*0.05, nrow(data[["kidney"]])),
                           title = 'nucleotide position',
                           autorange=FALSE,
                           ticks='outside',
                           tickwidth=2),
              bargap = 0.1,
              barmode = "overlay",
              plot_bgcolor = 'rgba(235,235,235,1)',
              legend = list(orientation = "h", bgcolor = "#08519c")
      )
    p6 <- plot_ly(data[["kidney"]], hoverinfo='none') %>%
      add_trace(x = ~position, y = ~X15, name = 'coverage', type = "bar", marker = list(color="#cb181d"), showlegend = FALSE) %>%
      layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
              xaxis = list(range=c(-nrow(data[["liver"]])*0.05, nrow(data[["kidney"]])),
                           title = 'nucleotide position',
                           autorange=FALSE,
                           ticks='outside',
                           tickwidth=2),
              bargap = 0.1,
              barmode = "overlay",
              plot_bgcolor = 'rgba(235,235,235,1)',
              legend = list(orientation = "h", bgcolor = "#08519c")
      )
    p7 <- plot_ly(data[["kidney"]], hoverinfo='none') %>%
      add_trace(x = ~position, y = ~X30, name = 'coverage', type = "bar", marker = list(color="#cb181d"), showlegend = FALSE) %>%
      layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
              xaxis = list(range=c(-nrow(data[["liver"]])*0.05, nrow(data[["kidney"]])),
                           title = 'nucleotide position',
                           autorange=FALSE,
                           ticks='outside',
                           tickwidth=2),
              bargap = 0.1,
              barmode = "overlay",
              plot_bgcolor = 'rgba(235,235,235,1)',
              legend = list(orientation = "h", bgcolor = "#08519c")
      )
    p8 <- plot_ly(data[["kidney"]], hoverinfo='none') %>%
      add_trace(x = ~position, y = ~X45, name = 'coverage', type = "bar", marker = list(color="#cb181d"), showlegend = FALSE) %>%
      layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
              xaxis = list(range=c(-nrow(data[["liver"]])*0.05, nrow(data[["kidney"]])),
                           title = 'nucleotide position',
                           autorange=FALSE,
                           ticks='outside',
                           tickwidth=2),
              bargap = 0.1,
              barmode = "overlay",
              plot_bgcolor = 'rgba(235,235,235,1)',
              legend = list(orientation = "h", bgcolor = "#08519c")
      )

    p9 <- plot_ly(data[["skeletal"]], hoverinfo='none') %>%
      add_trace(x = ~position, y = ~X0, name = 'coverage', type = "bar", marker = list(color="#238b45"), showlegend = FALSE) %>%
      layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
              xaxis = list(range=c(-nrow(data[["liver"]])*0.05, nrow(data[["skeletal"]])),
                           title = 'nucleotide position',
                           autorange=FALSE,
                           ticks='outside',
                           tickwidth=2),
              bargap = 0.1,
              barmode = "overlay",
              plot_bgcolor = 'rgba(235,235,235,1)',
              legend = list(orientation = "h", bgcolor = "#08519c")
      )
    p10 <- plot_ly(data[["skeletal"]], hoverinfo='none') %>%
      add_trace(x = ~position, y = ~X15, name = 'coverage', type = "bar", marker = list(color="#238b45"), showlegend = FALSE) %>%
      layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
              xaxis = list(range=c(-nrow(data[["liver"]])*0.05, nrow(data[["skeletal"]])),
                           title = 'nucleotide position',
                           autorange=FALSE,
                           ticks='outside',
                           tickwidth=2),
              bargap = 0.1,
              barmode = "overlay",
              plot_bgcolor = 'rgba(235,235,235,1)',
              legend = list(orientation = "h", bgcolor = "#08519c")
      )
    p11 <- plot_ly(data[["skeletal"]], hoverinfo='none') %>%
      add_trace(x = ~position, y = ~X30, name = 'coverage', type = "bar", marker = list(color="#238b45"), showlegend = FALSE) %>%
      layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
              xaxis = list(range=c(-nrow(data[["liver"]])*0.05, nrow(data[["skeletal"]])),
                           title = 'nucleotide position',
                           autorange=FALSE,
                           ticks='outside',
                           tickwidth=2),
              bargap = 0.1,
              barmode = "overlay",
              plot_bgcolor = 'rgba(235,235,235,1)',
              legend = list(orientation = "h", bgcolor = "#08519c")
      )
    p12 <- plot_ly(data[["skeletal"]], hoverinfo='none') %>%
      add_trace(x = ~position, y = ~X45, name = 'coverage', type = "bar", marker = list(color="#238b45"), showlegend = FALSE) %>%
      layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
              xaxis = list(range=c(-nrow(data[["liver"]])*0.05, nrow(data[["skeletal"]])),
                           title = 'nucleotide position',
                           autorange=FALSE,
                           ticks='outside',
                           tickwidth=2),
              bargap = 0.1,
              barmode = "overlay",
              plot_bgcolor = 'rgba(235,235,235,1)',
              legend = list(orientation = "h", bgcolor = "#08519c")
      )

    p <- subplot(p1, p5, p9, p2, p6, p10, p3, p7, p11, p4, p8, p12, nrows = 4, shareY = FALSE, shareX = TRUE, margin = c(0.01, 0.01, 0.02, 0.02))
  })
  
  output$coverage_plot <- renderPlotly({ plotInput()  })
  
  #export image locally
  output$imagesave1 <- downloadHandler(
    filename = function() { "Atf4.pdf" },
    content  = function(file) { export(plotInput(), file = file)  })
  
  output$imagesave2 <- downloadHandler(
    filename = function() { "Atf4.png" },
    content  = function(file) { export(plotInput(), file = file, vheight = 650, vwidth = 1000, zoom = 4)  })
  
}
ui <- fluidPage(   fluidRow(div(column(width=12, plotlyOutput("coverage_plot", height = 400))), style = "padding:25px;"),
                   fluidRow( downloadButton("imagesave1", "Save the image as pdf", style = "background-color:#084594; color:white"),
                             downloadButton("imagesave2", "Save the image as png", style = "background-color:#005a32; color:white"))
                   
)

shinyApp(ui = ui, server = server)






