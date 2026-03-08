### Baku Konferansi ici Kodlar
library(MASS)
library(mclust)
library(e1071)
library(ggplot2)
library(gridExtra)
library(clusterGeneration)
library(ppclust)
library(cluster) 
library(factoextra)
library(fcvalid)
library(dplyr)
library(pastecs)
library(patchwork)

###OUTLINE
#1) Kwon_mahalanobis 34.line
#2.1) OUTLIER EXPERIMENT 90.
#2.2) OUTLIER EXPERIMENT VISUALIZATION 127.
#2.3) OUTLIER EXPERIMENT COMPARISON BETWEEN INDICES 189.
#3) hangi sonuclari kullanalim?
#3.1) Monte Carlo 30 Trial Stress Test 266.
#3.2) Monte Carlo 50 Trial ST COMPARISON OF INDICES 314.
###CAUTION 3.2) takes a long time to run!
#3.3) STRESS TEST 427.
#3.4) STRESS TEST VISUALIZATION 468.
#4.1) Kwon Manual Calculation 533.
#4.2) Kwon_Mahalanobis Calculation 573.
#4.3) Kwon_Mahalanobis Calculation and Others Performance Comparison 624.
#5) Mahalanobis vs Euclidean the Visualization 665.

##Kwon_mahalanobis

kwon_mahalanobis <- function(data, k, U, centroids, m) {
  n <- nrow(data)
  p <- ncol(data)
  ## the numerator (ensures within-cluster compactness + penalty term)
  local_cov_matrices <- list()
  for (j in 1:k) {
    cluster_weights <- U[, j]^m
    S_j <- cov.wt(data, wt = cluster_weights)$cov
    S_j_reg <- S_j + diag(p) * 1e-6
    local_cov_matrices[[j]] <- S_j_reg
  }
  # the penalty calculation 
  S_global <- cov(data)
  centroid_mean <- colMeans(centroids)
  penalty_term <- 0
  for (j in 1:k) {
    penalty_term <- penalty_term + mahalanobis(centroids[j, ], center = centroid_mean, cov = S_global)
  }
  # within-cluster compactness
  compactness <- 0
  for (j in 1:k) {
    S_j_reg <- local_cov_matrices[[j]]
    for (i in 1:n) {
      d2_ij <- mahalanobis(data[i, ], center = centroids[j, ], cov = S_j_reg)
      compactness <- compactness + (U[i, j]^m) * d2_ij
    }
  }
  numerator <- compactness + (penalty_term / k)
  # denaminator
  min_centroid_dist_sq <- Inf
  for (j in 1:(k - 1)) {
    S_j_reg <- local_cov_matrices[[j]]
    for (l in (j + 1):k) {
      d2_jl <- mahalanobis(centroids[j, ], center = centroids[l, ], cov = S_global)
      if (d2_jl < min_centroid_dist_sq) {
        min_centroid_dist_sq <- d2_jl
      }
    }
  }
  if (min_centroid_dist_sq <= 0) {
    min_centroid_dist_sq <- .Machine$double.eps
  }
  return(numerator / min_centroid_dist_sq)
}




###a) hangi sonuclari kullanalim 
#monte carloyla 30 ve 50 kere kod tekrar ediyor ve uzun zaman aliyor




###3) aykiri deger deneyi icin kod 

###3.1) OUTLIER EXPERIMENT

outlier_rates <- c(0, 0.01, 0.05, 0.1, 0.15, 0.25) 
final_outlier_results <- data.frame()
set.seed(123)

for(rate in outlier_rates) {
  # Data generation (K=3)
  data_gen <- genRandomClust(numClust = 3, sepVal = 0.25, clustSizes = rep(200, 3), numNonNoisy = 2)
  X_clean <- data_gen$datList[[1]]
  
  # Injection of outliers
  num_outliers <- round(nrow(X_clean) * rate)
  if(num_outliers > 0) {
    # random points 
    outliers <- matrix(runif(num_outliers * 2, min = -50, max = 50), ncol = 2)
    X_noisy <- rbind(X_clean, outliers)
  } else {
    X_noisy <- X_clean
  }
  
  # Test for K values
  for(k_try in 2:5) {
    fcm <- cmeans(X_noisy, centers = k_try, m = 2)
    score <- kwon_mahalanobis(X_noisy, k_try, fcm$membership, fcm$centers, 2)
    final_outlier_results <- rbind(final_outlier_results, data.frame(Outlier_Rate = rate, K = k_try, Score = score))
  }
}

# Visualization
ggplot(final_outlier_results, aes(x=K, y=Score, color=as.factor(Outlier_Rate), group=Outlier_Rate)) +
  geom_line() + geom_point() + scale_y_log10() +
  labs(title="K Values For The Kwon_mahalanobis Index", color="Outlier Percentage") +
  theme_minimal()



### 3.2) OUTLIER EXPERIMENT VISUALIZATION (to see the geometric shape of clusters)
outlier_rates <- c(0, 0.01, 0.05, 0.10, 0.15, 0.25) 
viz_plots_outlier <- list()
counter <- 1

# LOOP
for(rate in outlier_rates) {
  set.seed(123) 
  # Cluster generation (N=600)
  # sepVal=0.25 for every cluster 
  data_gen <- genRandomClust(numClust = 3, sepVal = 0.25, clustszind = 3,
                             clustSizes = rep(200, 3), numNonNoisy = 2)
  X_clean <- data_gen$datList[[1]]
  labels_clean <- data_gen$memList[[1]] 
  
  # Outlier Injection
  # Aggregate to the initial dataset (control case n=600)
  num_outliers <- round(nrow(X_clean) * rate)
  if(num_outliers > 0) {
    # Noise spread around min = -40, max = 40
    outliers <- matrix(runif(num_outliers * 2, min = -40, max = 40), ncol = 2)
    # Merge
    X_viz <- rbind(X_clean, outliers)
    labels_viz <- c(labels_clean, rep("Noise", num_outliers))
  } else {
    X_viz <- X_clean
    labels_viz <- labels_clean
  }
  
  # Details of Drawing
  df_viz <- as.data.frame(X_viz)
  colnames(df_viz) <- c("V1", "V2")
  df_viz$Cluster <- as.factor(labels_viz)
  
  # Graph details
  p <- ggplot(df_viz, aes(x=V1, y=V2, color=Cluster)) +
    geom_point(alpha=0.7, size=1.2) + 
    
    # colors
    scale_color_manual(values = c("1" = "#F8766D",  
                                  "2" = "#00BA38",  
                                  "3" = "#619CFF",  
                                  "Noise" = "grey40")) + 
    
    labs(title = paste0("Outlier Percentage = %", rate * 100), 
         subtitle = paste("Total Data:", nrow(df_viz), "| Noise:", num_outliers)) +
    theme_minimal() +
    theme(legend.position = "none", 
          axis.title = element_blank(), 
          axis.text = element_blank(),  
          panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
          plot.title = element_text(size=11, face="bold"))
  
  viz_plots_outlier[[counter]] <- p
  counter <- counter + 1
}

# Outcome
grid.arrange(grobs = viz_plots_outlier, ncol = 3, 
             top = "Cluster Geometry In Increasing Outlier Percentage")


#3.3) OUTLIER EXPERIMENT COMPARISON BETWEEN INDICES

set.seed(123)
outlier_rates <- c(0, 0.01, 0.05, 0.10, 0.15, 0.25)
validity_results <- data.frame()

# current_rate instead of current_sep
for(current_rate in outlier_rates) {
  cat("Initializing: Outlier Rate = %", current_rate * 100, "...\n") 
  
  # Sep=0.25 ve N=600 (200x3) 
  data_gen <- genRandomClust(numClust = 3, sepVal = 0.25, clustszind = 3,
                             clustSizes = rep(200, 3), numNonNoisy = 2)
  X_clean <- data_gen$datList[[1]]
  
  # OUTLIER INJECTION
  num_outliers <- round(nrow(X_clean) * current_rate)
  if(num_outliers > 0) {
    # Random data points min = -40, max = 40
    outliers <- matrix(runif(num_outliers * 2, min = -40, max = 40), ncol = 2)
    # Aggregate two datasets
    X_noisy <- rbind(X_clean, outliers)
  } else {
    X_noisy <- X_clean }
  
  if(TRUE) { 
    dist_matrix <- dist(X_noisy) 
  } else {
    dist_matrix <- NULL }
  
  for(k_try in 2:5) {
    # The model (X_noisy)
    res_fcm <- ppclust::fcm(X_noisy, centers = k_try, m = 2)
    val_pc <- as.numeric(pc(res_fcm))
    val_pe <- as.numeric(pe(res_fcm))
    val_fs <- as.numeric(fs(res_fcm))
    val_xb <- as.numeric(xb(res_fcm))
    val_kwon_euc <- as.numeric(kwon(res_fcm))
    
    # silhouette
    if(!is.null(dist_matrix)) {
      sil_obj <- cluster::silhouette(res_fcm$cluster, dist_matrix)
      val_sil <- summary(sil_obj)$avg.width
    } else {
      val_sil <- NA
    }
    # Mahalanobis_kwon 
    val_kwon_mah <- kwon_mahalanobis(X_noisy, k_try, res_fcm$u, res_fcm$v, 2)
    
    # Outlier rates
    validity_results <- rbind(validity_results, data.frame(
      Outlier_Rate = current_rate,
      K = k_try,
      PC = val_pc,
      PE = val_pe,
      Fukuyama_Sugeno = val_fs,
      Xie_Beni = val_xb,
      Kwon_Euclidean = val_kwon_euc,
      Silhouette = val_sil,
      Kwon_Mahalanobis = val_kwon_mah
    ))
  }
}

# The results
unique_rates <- unique(validity_results$Outlier_Rate)
for(r_val in unique_rates) {
  cat("\n=========================================\n")
  cat(" Results: Outlier Percentage = %", r_val * 100, "\n")
  cat("=========================================\n")
  subset_table <- validity_results[validity_results$Outlier_Rate == r_val, ]
  print(subset_table, row.names = FALSE)
  cat("\n") 
}



###3.1) Monte Carlo 30 Trial Stress Test  

n_clean <- 600
k_true <- 3
sep_val <- 0.25 
ranges <- c(50, 75, 100) # Noise Range
rates <- c(0.15, 0.30, 0.45) # Noise Percentage 

n_sim <- 30 # Monte Carlo Trial Number 
final_stress_test <- data.frame()
set.seed(42)

data_clean <- genRandomClust(numClust = k_true, sepVal = sep_val, clustszind = 3,
                             clustSizes = rep(n_clean/k_true, k_true), 
                             numNonNoisy = 2)$datList[[1]]
for(r_val in ranges) {
  for(rate in rates) {
    for(sim_i in 1:n_sim) {
      #  Noise Injection
      n_out <- round(n_clean * rate)
      outliers <- matrix(runif(n_out * 2, min = -r_val, max = r_val), ncol = 2)
      X_noisy <- rbind(data_clean, outliers)
      # k=3
      scores <- sapply(2:5, function(k) {
        fcm <- cmeans(X_noisy, centers = k, m = 2)
        kwon_mahalanobis(X_noisy, k, fcm$membership, fcm$centers, 2)
      })
      # Optimal value of kwon_mahal
      k_predicted <- (2:5)[which.min(scores)]
      
      final_stress_test <- rbind(final_stress_test, 
                                 data.frame(Noise_Range = r_val, Noise_Rate_yuzdelik = rate*100, 
                                            Predicted_K = k_predicted,
                                            Success = (k_predicted == 3)))}}}

print("--- RESULT TABLE ---")
print(final_stress_test)


sonuc_tablosu <- final_stress_test %>%
  group_by(Noise_Range, Noise_Rate_yuzdelik) %>%
  summarise(
    Basari_Yuzdesi = mean(Success) * 100 
  )
print("--- RESULT TABLE (SUCCESS RATES) ---")
print(sonuc_tablosu)


###3.2) Monte Carlo 50 Trial ST COMPARISON OF INDICES

n_clean <- 600
k_true <- 3
sep_val <- 0.25 
ranges <- c(50, 75, 100) 
rates <- c(0.15, 0.30, 0.45) 
n_sim <- 50

results_df <- data.frame() 

set.seed(123) 

total_iter <- length(ranges) * length(rates) * n_sim
counter <- 0

print("Simulasyon basliyor...")

for(r_val in ranges) {
  for(rate in rates) {
    
    # Bu senaryo icin basari sayaclari
    success_counts <- c(PC=0, PE=0, FS=0, XB=0, SIL=0, KWON_EUC=0, KWON_MAH=0)
    
    for(sim in 1:n_sim) {
      counter <- counter + 1
      if(counter %% 10 == 0) cat(paste0("\rTamamlanan: ", counter, "/", total_iter, " ... "))
      
      # 1. Veri Uretimi
      data_clean <- genRandomClust(numClust = k_true, sepVal = sep_val, clustszind = 3, 
                                   clustSizes = rep(n_clean/k_true, k_true), 
                                   numNonNoisy = 2)$datList[[1]]
      
      # 2. Gurultu Ekleme
      n_out <- round(n_clean * rate)
      outliers <- matrix(runif(n_out * 2, min = -r_val, max = r_val), ncol = 2)
      X_sim <- rbind(data_clean, outliers) # Senin kodundaki degisken ismi X_sim
      
      # 3. Mesafe Matrisi (Senin kodun icin gerekli)
      dist_matrix <- dist(X_sim)
      
      # Skorlari tutacak gecici tablo
      indices_res <- data.frame(k=2:5, PC=NA, PE=NA, FS=NA, XB=NA, SIL=NA, KE=NA, KM=NA)
      
      for(k_idx in 1:4) { 
        k_try <- k_idx + 1 # k=2,3,4,5
        
        # --- SENIN GONDERDIGIN KOD BLOGU ---
        
        # Modelleme
        res_fcm <- ppclust::fcm(X_sim, centers = k_try, m = 2)
        
        val_pc <- as.numeric(pc(res_fcm))
        val_pe <- as.numeric(pe(res_fcm))
        val_fs <- as.numeric(fs(res_fcm))
        val_xb <- as.numeric(xb(res_fcm))
        val_kwon_euc <- as.numeric(kwon(res_fcm)) # ppclust icindeki standart Kwon
        
        # Siluet cluster paketi
        if(!is.null(dist_matrix)) {
          sil_obj <- cluster::silhouette(res_fcm$cluster, dist_matrix)
          val_sil <- summary(sil_obj)$avg.width
        } else {
          val_sil <- NA
        }
        
        # --- OZEL FONKSIYONUNUZ (Environment'da yuklu oldugu varsayiliyor) ---
        # Burada sadece cagri yapiyoruz, tanimlama yok.
        val_kwon_mah <- kwon_mahalanobis(X_sim, k_try, res_fcm$u, res_fcm$v, 2)
        
        # --- SONUCLARI TABLOYA YAZ ---
        indices_res$PC[k_idx] <- val_pc
        indices_res$PE[k_idx] <- val_pe
        indices_res$FS[k_idx] <- val_fs
        indices_res$XB[k_idx] <- val_xb
        indices_res$SIL[k_idx] <- val_sil
        indices_res$KE[k_idx] <- val_kwon_euc
        indices_res$KM[k_idx] <- val_kwon_mah
      }
      
      # --- KARAR MEKANIZMASI ---
      # PC ve SIL -> Maksimum, Digerleri -> Minimum
      pred_k <- c(
        PC = indices_res$k[which.max(indices_res$PC)],
        PE = indices_res$k[which.min(indices_res$PE)],
        FS = indices_res$k[which.min(indices_res$FS)],
        XB = indices_res$k[which.min(indices_res$XB)],
        SIL = indices_res$k[which.max(indices_res$SIL)],
        KWON_EUC = indices_res$k[which.min(indices_res$KE)],
        KWON_MAH = indices_res$k[which.min(indices_res$KM)]
      )
      
      # Basariyi kaydet (Dogru cevap 3 ise 1, degilse 0 ekle)
      success_counts <- success_counts + (pred_k == k_true)
    }
    
    # Yuzdelik basari hesabi
    success_rates <- (success_counts / n_sim) * 100
    temp_row <- c(Range=r_val, Rate=rate, success_rates)
    results_df <- rbind(results_df, temp_row)
  }
}

# Sutun isimlerini duzenle
colnames(results_df) <- c("Range", "Rate", "Acc_PC", "Acc_PE", "Acc_FS", "Acc_XB", "Acc_SIL", "Acc_KWON_EUC", "Acc_KWON_MAH")

print(" ")
print("--- STRES TESTİ SONUÇLARI ---")
print(results_df)




###3.3) STRESS TEST 
n_clean <- 600
k_true <- 3
sep_val <- 0.25 # Seperation value
ranges <- c(50, 75, 100) # Noise Range
rates <- c(0.15, 0.30, 0.45) # Noise Percentage

final_stress_test <- data.frame()
set.seed(42)

data_clean <- genRandomClust(numClust = k_true, sepVal = sep_val, clustszind = 3,
                             clustSizes = rep(n_clean/k_true, k_true), 
                             numNonNoisy = 2)$datList[[1]]

for(r_val in ranges) {
  for(rate in rates) {
    #  Noise
    n_out <- round(n_clean * rate)
    outliers <- matrix(runif(n_out * 2, min = -r_val, max = r_val), ncol = 2)
    X_noisy <- rbind(data_clean, outliers)
    
    # k=3
    scores <- sapply(2:5, function(k) {
      fcm <- cmeans(X_noisy, centers = k, m = 2)
      kwon_mahalanobis(X_noisy, k, fcm$membership, fcm$centers, 2)
    })
    
    # Optimum values 
    k_predicted <- (2:5)[which.min(scores)]
    
    final_stress_test <- rbind(final_stress_test, 
                               data.frame(Noise_Range = r_val, Noise_Rate_yuzdelik = rate*100, 
                                          Predicted_K = k_predicted,
                                          Success = (k_predicted == 3)))
  }
}
# Results
print(final_stress_test)



###3.4) STRESS TEST VISUALIZATION 

# Empty list
plot_list <- list()

for(i in 1:nrow(final_stress_test)) {
  
  # The Parameters
  current_range <- final_stress_test$Noise_Range[i]
  current_rate <- final_stress_test$Noise_Rate_yuzdelik[i] / 100 # Yüzdeyi orana çevir
  pred_k <- final_stress_test$Predicted_K[i]
  is_success <- final_stress_test$Success[i]
  
  # Noise Injection
  # Seed changes in every step, the points may not be the same 
  # Range ve Rate stays the same 
  n_out <- round(n_clean * current_rate)
  outliers <- matrix(runif(n_out * 2, min = -current_range, max = current_range), ncol = 2)
  
  # Merge datasets
  X_display <- rbind(data_clean, outliers)
  
  #Data frame
  df_plot <- as.data.frame(X_display)
  colnames(df_plot) <- c("V1", "V2")
  
  # Colors
  plot_color <- if(is_success) "steelblue" else "firebrick"
  # If wrong give the true k value
  k_label <- if(is_success) paste0("K=", pred_k) else paste0("ERROR (K=", pred_k, ")")
  
  main_title <- paste0("Range: ", current_range, " | Percentage: %", current_rate*100)
  sub_title <- k_label
  
  # The Graph
  p <- ggplot(df_plot, aes(x=V1, y=V2)) +
    geom_point(alpha=0.6, size=0.5, color=plot_color) +
    labs(title = main_title, subtitle = sub_title) +
    theme_minimal() +
    theme(
      plot.title = element_text(size=9, face="bold"),
      plot.subtitle = element_text(size=8, color=if(is_success) "black" else "red"),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    )
  
  # Add to list
  plot_list[[i]] <- p
}

# 3x3 
grid.arrange(grobs = plot_list, ncol = 3, 
             top = "Kwon-Mahalanobis Stress Test: Visualization Of 9 Datasets")





###4) Seeds Dataset Analysis & Others 
numeric_2 <- subset(seeds_dataset.2,select = -c(V8), drop=TRUE )  
names(numeric_2) <- c('area', 'parameter', 'compactness', 'kernel lenght', 'kernel width', 'assymetryCoef', 'grooveLength')
View(numeric_2)

###4.1) Kwon Manual Calculation 
set.seed(123)
kwon_euclidean <- function(numeric_2, k, U, centroids, m) {
  n <- nrow(numeric_2)
  p <- ncol(numeric_2)
  ## the numerator (ensures within-cluster compactness + penalty term)
  # the penalty calculation 
  centroid_mean <- colMeans(centroids)
  penalty_term <- 0
  for (j in 1:k) {
    penalty_term <- penalty_term + sum((centroids[j, ] - centroid_mean)^2)
  }
  penalty_term <- penalty_term / k
  # within-cluster compactness
  numerator <- 0
  for (i in 1:n) {
    for (j in 1:k) {
      dist_ij <- sqrt(sum((numeric_2[i, ] - centroids[j, ])^2))
      numerator <- numerator + (U[i, j]^m) * (dist_ij^2)
    }
  }
  numerator <- numerator + penalty_term 
  # denaminator
  min_centroid_dist_sq <- Inf
  for (j in 1:(k - 1)) {
    for (l in (j + 1):k) {
      dist_jl_sq <- sum((centroids[j, ] - centroids[l, ])^2)
      if (dist_jl_sq < min_centroid_dist_sq) {
        min_centroid_dist_sq <- dist_jl_sq
      }
    }
  }
  
  if (min_centroid_dist_sq <= 0) {
    min_centroid_dist_sq <- .Machine$double.eps
  }
  return(numerator / min_centroid_dist_sq)
}


###4.2) Kwon_Mahalanobis Calculation 
kwon_mahalanobis <- function(numeric_2, k, U, centroids, m) {
  n <- nrow(numeric_2)
  p <- ncol(numeric_2)
  
  ## the numerator (ensures within-cluster compactness + penalty term)
  local_cov_matrices <- list()
  for (j in 1:k) {
    cluster_weights <- U[, j]^m
    S_j <- cov.wt(numeric_2, wt = cluster_weights)$cov
    S_j_reg <- S_j + diag(p) * 1e-6
    local_cov_matrices[[j]] <- S_j_reg
  }
  
  # the penalty calculation 
  S_global <- cov(numeric_2)
  centroid_mean <- colMeans(centroids)
  penalty_term <- 0
  for (j in 1:k) {
    penalty_term <- penalty_term + mahalanobis(centroids[j, ], center = centroid_mean, cov = S_global)
  }
  # within-cluster compactness
  compactness <- 0
  for (j in 1:k) {
    S_j_reg <- local_cov_matrices[[j]]
    for (i in 1:n) {
      d2_ij <- mahalanobis(numeric_2[i, ], center = centroids[j, ], cov = S_j_reg)
      compactness <- compactness + (U[i, j]^m) * d2_ij
    }
  }
  numerator <- compactness + (penalty_term / k)
  
  # denaminator
  min_centroid_dist_sq <- Inf
  for (j in 1:(k - 1)) {
    S_j_reg <- local_cov_matrices[[j]]
    for (l in (j + 1):k) {
      d2_jl <- mahalanobis(centroids[j, ], center = centroids[l, ], cov = S_global)
      if (d2_jl < min_centroid_dist_sq) {
        min_centroid_dist_sq <- d2_jl
      }
    }
  }
  
  if (min_centroid_dist_sq <= 0) {
    min_centroid_dist_sq <- .Machine$double.eps
  }
  return(numerator / min_centroid_dist_sq)
}


###4.3) Kwon_Mahalanobis Calculation and Others Performance Comparison 
#comparison with ppclust library 
kwon_values <- numeric() #very important!
kwon_builtin_values <- numeric() #very important!
kwon_values_fcm_mahal <- numeric()
kwon_values_gk_mahal <- numeric()
m <- 2
set.seed(123)
for (k in 2:7) {
  fcm_result <- fcm(numeric_2, centers = k, m = m)
  gk_result <- gk(numeric_2, centers = k, m = m)
  U <- fcm_result$u
  V <- fcm_result$v
  # Manual Kwon calculation
  kwon_val_manual <- kwon_euclidean(numeric_2, k, U, V, m)
  kwon_values <- c(kwon_values, kwon_val_manual)
  # Built-in Kwon index from ppclust
  kwon_val_builtin <- kwon(fcm_result)
  kwon_builtin_values <- c(kwon_builtin_values, kwon_val_builtin)
  # Mahalanobis Kwon calculation with fcm
  kwon_val_kwon_fcm <- kwon_mahalanobis(numeric_2, k, U, V, m)
  kwon_values_fcm_mahal <- c(kwon_values_fcm_mahal, kwon_val_kwon_fcm)
  # Mahalanobis Kwon calculation with gk
  kwon_val_kwon <- kwon_mahalanobis(numeric_2, k, gk_result$u, gk_result$v, m)
  kwon_values_gk_mahal <- c(kwon_values_gk_mahal, kwon_val_kwon)
}


results_df_seeds <- data.frame(
  Clusters = 2:7,
  Kwon_Euclidean = round(kwon_values, 4),
  Kwon_Builtin = round(kwon_builtin_values, 4),
  Kwon_Mahalanobis_fcm = round(kwon_values_fcm_mahal, 4),
  Kwon_Mahalanobis_gk = round(kwon_values_gk_mahal, 4)
) ## resetting the object kwon values is very important

print(results_df_seeds, row.names = FALSE)




###5) Mahalanobis vs Euclidean the Visualization
# 
set.seed(42) 
n_points <- 100 
mu <- c(10, 10) 

# the X vector
x_vec <- rnorm(n_points, mean = mu[1], sd = 2) 

# Dataset creation (elipsoid with positive correlation)
data_scatter <- data.frame(
  x = x_vec,
  y = 0.7 * x_vec + rnorm(n_points, mean = 3, sd = 0.8) 
)

# Center Point
center_point <- data_scatter %>%
  summarise(mean_x = mean(x), mean_y = mean(y))
center_point$name <- "Center Point"
center_point$color <- "orange"

# Target Points (green & pink)
euc_dist_target <- 3.0 
# Euclidean: d = sqrt(dx^2 + dy^2)

# Target 1 (green): Has less value in terms of Mahalanobis Large dx, small dy so it moves along the long axis of the ellipse
dx_green <- 2.5 # large displacement in X (along the direction data naturally varies)
dy_green <- sqrt(euc_dist_target^2 - dx_green^2) # small displacement in Y (solved from Euclidean formula)
target_green <- data.frame(
  mean_x = center_point$mean_x + dx_green,
  mean_y = center_point$mean_y + dy_green,
  name = "Target Green", color = "green"
)

# Target 2 (pink): Euclidean stays the same, but mahal dist is high
dx_pink <- 1.0 # small displacement in X
dy_pink <- sqrt(euc_dist_target^2 - dx_pink^2) #large displacement in Y (≈2.83) → pushes point outside ellipse 
target_pink <- data.frame(
  mean_x = center_point$mean_x + dx_pink,
  mean_y = center_point$mean_y + dy_pink,
  name = "Target Pink", color = "pink"
)


special_points <- bind_rows(center_point, target_green, target_pink)

# Distance Calculation
S <- cov(data_scatter[, c("x", "y")])
S_inv <- solve(S) # Inverse of Cov matrix

mu_vec <- as.matrix(center_point[, c("mean_x", "mean_y")])

calculate_distances <- function(target_point, mu_vec, S_inv) {
  x_vec <- as.matrix(target_point[, c("mean_x", "mean_y")])
  
  #  The Euclidean Calculation
  diff_euc <- x_vec - mu_vec
  euc_dist <- sqrt(sum(diff_euc^2))
  
  # The Mahalanobis Calculation
  diff_mah_row <- x_vec - mu_vec
  intermediate_result <- diff_mah_row %*% S_inv
  mah_sq <- intermediate_result %*% t(diff_mah_row)
  mah_dist <- sqrt(mah_sq[1, 1])
  
  return(data.frame(
    Target = target_point$name,
    Euclidean = round(euc_dist, 2),
    Mahalanobis = round(mah_dist, 2)
  ))
}

distances_green <- calculate_distances(target_green, mu_vec, S_inv)
distances_pink <- calculate_distances(target_pink, mu_vec, S_inv)
results_table <- bind_rows(distances_green, distances_pink)

# Visualization

# Scatter Plot
p_scatter <- ggplot(data_scatter, aes(x = x, y = y)) +
  geom_point(color = "#4c84c4", size = 2, alpha = 0.7) +
  
  # The ellipse
  stat_ellipse(type = "norm", level = 0.95, color = "gray50", linetype = "dashed") +
  
  # Special features and the arrows
  geom_point(data = special_points, aes(x = mean_x, y = mean_y),
             color = special_points$color, size = 4) +
  
  geom_segment(aes(x = center_point$mean_x, y = center_point$mean_y,
                   xend = target_green$mean_x, yend = target_green$mean_y),
               arrow = arrow(length = unit(0.3, "cm")), color = "green", linewidth = 0.8) +
  
  geom_segment(aes(x = center_point$mean_x, y = center_point$mean_y,
                   xend = target_pink$mean_x, yend = target_pink$mean_y),
               arrow = arrow(length = unit(0.3, "cm")), color = "pink", linewidth = 0.8) +
  
  labs(x = "X", y = "Y", title = "The Comparison of Euclidean and Mahalanobis") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# The outcome table 
table_data <- results_table
# The euclidean scores
table_data$Euc_Color <- "red"
table_data$Mah_Color <- ifelse(table_data$Mahalanobis == min(table_data$Mahalanobis), "darkgreen", "red")

p_table <- ggplot() +
  # Titles
  annotate("text", x = 1, y = 3.5, label = "Euclidean", size = 4.5, fontface = "bold") +
  annotate("text", x = 2, y = 3.5, label = "Mahalanobis", size = 4.5, fontface = "bold") +
  
  # Target green
  annotate("point", x = 0.5, y = 2.5, color = "green", size = 5) +
  annotate("text", x = 1, y = 2.5, label = table_data$Euclidean[1], size = 5, color = table_data$Euc_Color[1]) +
  annotate("text", x = 2, y = 2.5, label = table_data$Mahalanobis[1], size = 5, color = table_data$Mah_Color[1]) +
  
  # Target pink
  annotate("point", x = 0.5, y = 1.5, color = "pink", size = 5) +
  annotate("text", x = 1, y = 1.5, label = table_data$Euclidean[2], size = 5, color = table_data$Euc_Color[2]) +
  annotate("text", x = 2, y = 1.5, label = table_data$Mahalanobis[2], size = 5, color = table_data$Mah_Color[2]) +
  
  
  # The table
  geom_segment(aes(x = 1.5, y = 0, xend = 1.5, yend = 4), linetype = "dotted") +
  
  xlim(0, 3) + ylim(0, 4) +
  theme_void()

# Result
final_plot <- p_scatter / p_table + 
  plot_layout(ncol = 1, heights = c(3, 1.5))

print(final_plot)