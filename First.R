data <- read.csv("survey lung cancer.csv", sep=";", header=TRUE, stringsAsFactors = FALSE)
colnames(data) <- trimws(colnames(data))
data$LUNG_CANCER <- ifelse(data$LUNG_CANCER == "YES", 1, 0)
Y_list <- c("LUNG_CANCER", "AGE", "SMOKING", "YELLOW_FINGERS", "CHEST_PAIN", "SHORTNESS.OF.BREATH")
dir.create("lab1", showWarnings = FALSE)
results_table <- data.frame()
for (Yi in Y_list) {
  current_val <- as.numeric(data[[Yi]])
  current_val <- current_val[!is.na(current_val)]
  if (all(current_val %in% c(1, 2))) {
    current_val <- current_val - 1
  }
  n <- length(current_val)
  if (Yi == "AGE") {
    sorted_val <- sort(current_val)
    Q1 <- sorted_val[floor(n * 0.25)]
    Q3 <- sorted_val[ceiling(n * 0.75)]
    IQR_val <- Q3 - Q1
    lower_bound <- Q1 - 1.5 * IQR_val
    upper_bound <- Q3 + 1.5 * IQR_val
    current_val <- current_val[current_val >= lower_bound & current_val <= upper_bound]
    n <- length(current_val)
    cat("AGE: удалено", 200 - n, "выбросов, осталось", n, "наблюдений\n")
  }
  m_mean   <- sum(current_val) / n
  m_var    <- sum((current_val - m_mean)^2) / (n - 1)
  m_sd     <- sqrt(m_var)
  sorted_val <- sort(current_val)
  m_median   <- (sorted_val[n/2] + sorted_val[n/2 + 1]) / 2
  tbl <- table(current_val)
  m_mode <- as.numeric(names(tbl)[tbl == max(tbl)])[1]
  is_01 <- all(current_val %in% c(0, 1))
  png(paste0("lab1/", Yi, ".png"), width = 800, height = 600)
  if (is_01) {
    brks <- c(-0.5, 0.5, 1.5)
    hist(current_val, breaks = brks,
         main = paste("Гистограмма:", Yi),
         col = "lightgreen", border = "black",
         xlab = Yi, ylab = "Частота", xaxt = "n")
    axis(1, at = c(0, 1), labels = c("0", "1"))
    h <- hist(current_val, breaks = brks, plot = FALSE)
    lines(h$mids, h$counts, type = "b", col = "red", lwd = 2)
  } else {
    n_breaks <- 15
    brks <- seq(min(current_val), max(current_val), length.out = n_breaks + 1)
    h <- hist(current_val, breaks = brks, plot = FALSE)
    xlim <- range(brks) + c(-1, 1) * diff(range(brks)) * 0.02
    ylim <- c(0, max(h$counts) * 1.15)
    hist(current_val, breaks = brks,
         main = paste("Гистограмма:", Yi),
         col = "lightgreen", border = "black",
         xlab = Yi, ylab = "Частота",
         xlim = xlim, ylim = ylim, xaxt = "n")
    axis(1, at = round(brks, 0))
    lines(h$mids, h$counts, type = "b", col = "red", lwd = 2)
  }
  dev.off()
  new_row <- data.frame(
    Признак = Yi,
    Среднее = round(m_mean, 3),
    Медиана = round(m_median, 3),
    Мода = round(m_mode, 3),
    СКО = round(m_sd, 3),
    Дисперсия = round(m_var, 3)
  )
  results_table <- rbind(results_table, new_row)
}
cat("\n====== МЕТРИКИ ======\n")
print(results_table, row.names = FALSE)

write.csv(results_table, "lab1/metrics_table.csv", row.names = FALSE)
