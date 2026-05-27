# Лабораторная работа №2: Интервальное оценивание и проверка гипотез
# Выполнил: [Имя]

# 1. Загрузка и предобработка данных
data <- read.csv("survey lung cancer.csv", sep=";", header=TRUE, stringsAsFactors = FALSE)
colnames(data) <- trimws(colnames(data))

# Рекодирование целевого признака X (LUNG_CANCER)
data$LUNG_CANCER <- ifelse(data$LUNG_CANCER == "YES", 1, 0)

# Список влияющих признаков Yi
Yi_list <- c("AGE", "SMOKING", "YELLOW_FINGERS", "CHEST_PAIN", "SHORTNESS.OF.BREATH")

# Очистка и преобразование (аналогично Лаб 1 для консистентности)
processed_data <- data[, c("LUNG_CANCER", Yi_list)]
for (col in Yi_list) {
  processed_data[[col]] <- as.numeric(processed_data[[col]])
  # Рекодирование {1,2} в {0,1} для бинарных
  if (all(processed_data[[col]] %in% c(1, 2), na.rm = TRUE)) {
    processed_data[[col]] <- processed_data[[col]] - 1
  }
}

# Обработка выбросов для AGE (1.5 * IQR)
age_vals <- processed_data$AGE
n_initial <- length(age_vals)
q1_idx <- floor(n_initial * 0.25)
q3_idx <- ceiling(n_initial * 0.75)
sorted_age <- sort(age_vals)
Q1 <- sorted_age[q1_idx]
Q3 <- sorted_age[q3_idx]
IQR_val <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val

# Фильтруем весь датасет по выбросам AGE
processed_data <- processed_data[processed_data$AGE >= lower_bound & processed_data$AGE <= upper_bound, ]
n <- nrow(processed_data)
cat("Удалено выбросов AGE:", n_initial - n, "\n")

# ==============================================================================
# ЧАСТЬ 1: ПРОВЕРКА ГИПОТЕЗ (ДЛЯ X = LUNG_CANCER)
# ==============================================================================

# Гипотеза H0: Среднее значение X (доля больных) равно 0.5 (равновероятность)
# H1: Среднее значение X не равно 0.5
mu_0 <- 0.5
alpha <- 0.05
X <- processed_data$LUNG_CANCER

# Ручной расчет метрик
x_mean <- sum(X) / n
x_var  <- sum((X - x_mean)^2) / (n - 1)
x_sd   <- sqrt(x_var)
x_se   <- x_sd / sqrt(n)

# t-статистика
t_stat <- (x_mean - mu_0) / x_se

# Границы доверительного интервала (alpha = 0.05)
# Используем qt для получения критического значения (согласно методичке допускается для проверки)
t_crit <- qt(1 - alpha/2, df = n - 1)
ci_lower <- x_mean - t_crit * x_se
ci_upper <- x_mean + t_crit * x_se

# P-value
p_value <- 2 * (1 - pt(abs(t_stat), df = n - 1))

cat("\n--- ПРОВЕРКА ГИПОТЕЗЫ (X = LUNG_CANCER) ---\n")
cat("H0: mu =", mu_0, "\n")
cat("Выборочное среднее:", round(x_mean, 4), "\n")
cat("t-статистика:", round(t_stat, 4), "\n")
cat("P-value:", round(p_value, 6), "\n")
cat("Доверительный интервал [", alpha, "]: [", round(ci_lower, 4), ",", round(ci_upper, 4), "]\n")

if (p_value < alpha) {
  cat("Вывод: Гипотеза H0 отвергается (среднее значимо отличается от", mu_0, ")\n")
} else {
  cat("Вывод: Гипотеза H0 не отвергается\n")
}

# Дополнительно: проверка для других alpha
alphas <- c(0.01, 0.1)
for (a in alphas) {
  tc <- qt(1 - a/2, df = n - 1)
  cat("При alpha =", a, "критическое t =", round(tc, 4), 
      "Результат:", ifelse(abs(t_stat) > tc, "Отвергаем H0", "Принимаем H0"), "\n")
}

# ==============================================================================
# ЧАСТЬ 2: ЛИНЕЙНАЯ РЕГРЕССИЯ С НУЛЯ
# ==============================================================================

cat("\n--- ЛИНЕЙНАЯ РЕГРЕССИЯ (Градиентный спуск) ---\n")

# 5. Нормализация признаков (Z-score)
normalize <- function(vec) {
  m <- sum(vec) / length(vec)
  s <- sqrt(sum((vec - m)^2) / (length(vec) - 1))
  if (s == 0) return(vec - m)
  return((vec - m) / s)
}

# Матрица признаков Y (Yi_list)
Y_raw <- as.matrix(processed_data[, Yi_list])
Y_norm <- apply(Y_raw, 2, normalize)

# Добавляем колонку-смещение (intercept)
Y_norm <- cbind(1, Y_norm)
colnames(Y_norm)[1] <- "Intercept"

# Целевая переменная X
X_target <- processed_data$LUNG_CANCER

# Разделение на обучающую и тестовую выборки (80/20)
set.seed(42) # Для воспроизводимости
train_idx <- sample(1:n, size = floor(0.8 * n))
Y_train <- Y_norm[train_idx, ]
X_train <- X_target[train_idx]
Y_test  <- Y_norm[-train_idx, ]
X_test  <- X_target[-train_idx]

# 6. Реализация градиентного спуска
# Функция потерь: MSE = (1/2m) * sum((pred - actual)^2)
compute_cost <- function(Y, X, theta) {
  m <- length(X)
  preds <- Y %*% theta
  cost <- sum((preds - X)^2) / (2 * m)
  return(cost)
}

# Параметры градиентного спуска
alpha_lr <- 0.01  # Скорость обучения
epsilon  <- 1e-6 # Порог остановки
max_iter <- 5000
theta <- rep(0, ncol(Y_train)) # Инициализация весов
cost_history <- numeric(max_iter)
m_train <- length(X_train)

iter <- 0
delta_cost <- Inf

while (iter < max_iter && delta_cost > epsilon) {
  iter <- iter + 1
  preds <- Y_train %*% theta
  errors <- preds - X_train
  
  # Расчет градиента: (1/m) * Y^T * (Y*theta - X)
  gradient <- (t(Y_train) %*% errors) / m_train
  
  # Обновление весов
  theta_new <- theta - alpha_lr * gradient
  
  new_cost <- compute_cost(Y_train, X_train, theta_new)
  cost_history[iter] <- new_cost
  
  if (iter > 1) {
    delta_cost <- abs(cost_history[iter-1] - cost_history[iter])
  }
  
  theta <- theta_new
}

cost_history <- cost_history[1:iter]

# Визуализация функции потерь
dir.create("lab2", showWarnings = FALSE)
png("lab2/loss_history.png", width = 800, height = 600)
plot(1:iter, cost_history, type = "l", col = "blue", lwd = 2,
     main = "История функции потерь (MSE)",
     xlab = "Итерации", ylab = "Потери")
dev.off()

cat("Остановка после", iter, "итераций.\n")
cat("Веса (theta):\n")
print(round(theta, 4))

# 7. Функция предсказания
predict_manual <- function(Y, theta) {
  return(Y %*% theta)
}

# 8. Метрики качества на тестовой выборке
preds_test <- predict_manual(Y_test, theta)
m_test <- length(X_test)

mse_test  <- sum((preds_test - X_test)^2) / m_test
rmse_test <- sqrt(mse_test)

# R^2 = 1 - SSR/SST
rss <- sum((X_test - preds_test)^2)
mean_x_test <- sum(X_test) / m_test
tss <- sum((X_test - mean_x_test)^2)
r_squared <- 1 - (rss / tss)

cat("\nМетрики на тестовой выборке:\n")
cat("MSE:", round(mse_test, 6), "\n")
cat("RMSE:", round(rmse_test, 6), "\n")
cat("R^2:", round(r_squared, 6), "\n")

# 9. Сравнение с встроенной функцией lm()
model_lm <- lm(LUNG_CANCER ~ ., data = as.data.frame(cbind(LUNG_CANCER = X_train, Y_train[, -1])))
cat("\nКоэффициенты встроенной модели lm():\n")
print(round(coef(model_lm), 4))

# Предсказание lm на тесте
preds_lm <- predict(model_lm, as.data.frame(Y_test[, -1]))
mse_lm <- sum((preds_lm - X_test)^2) / m_test
cat("MSE (lm):", round(mse_lm, 6), "\n")

write.csv(data.frame(Metric=c("MSE", "RMSE", "R2"), Value=c(mse_test, rmse_test, r_squared)), 
          "lab2/metrics_regression.csv", row.names = FALSE)
