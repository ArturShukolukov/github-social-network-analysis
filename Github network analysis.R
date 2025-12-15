# ------------------------------------------------------------------
# Анализ социальной сети GitHub (igraph + poweRlaw)
# Данные: musae_git_edges.csv, musae_git_target.csv
# Сеть: неориентированная, невзвешенная. Узлы — пользователи, рёбра — взаимная подписка.
# ------------------------------------------------------------------


library(igraph)
library(poweRlaw)


# ---------------------- 1) Загрузка данных ----------------------

edges <- read.csv("musae_git_edges.csv", stringsAsFactors = FALSE)


g <- graph_from_data_frame(edges, directed = FALSE)

cat("Вершин:", vcount(g), "\n")
cat("Рёбер:", ecount(g), "\n")

# Проверка на связность
cat(is_connected(g))

# ---------------------- 2) Базовые метрики ----------------------

deg <- degree(g)

diam <- diameter(g, directed = FALSE)
apl  <- mean_distance(g, directed = FALSE)
avg_deg <- mean(deg)

clust_global <- transitivity(g, type = "global")
assort_deg <- assortativity_degree(g, directed = FALSE)
dens <- edge_density(g)

cat("Плотность:", dens, "\n")
cat("Средняя степень:", avg_deg, "\n")
cat("Диаметр:", diam, "\n")
cat("Средняя длина пути:", apl, "\n")
cat("Коэффициент кластеризации (global):", clust_global, "\n")
cat("Ассортативность по степени:", assort_deg, "\n\n")

# ---------------------- 2) Распределение степеней (CCDF) ----------------------

sorted_deg <- sort(deg)
ccdf_y <- 1 - ecdf(sorted_deg)(sorted_deg)
ccdf_x <- sorted_deg

ccdf_df <- data.frame(k = ccdf_x, ccdf = ccdf_y)
ccdf_df <- ccdf_df[!duplicated(ccdf_df$k), ]
ccdf_df <- ccdf_df[ccdf_df$ccdf > 0, ]

plot(ccdf_df$k, ccdf_df$ccdf, log = "xy",
     xlab = "k (степень вершины)",
     ylab = "CCDF",
     main = "CCDF распределения степеней (GitHub)")

# ---------------------- 3) Оценка хвоста степенного распределения (poweRlaw) ----------------------

m_pl <- displ$new(deg)
est_xmin <- estimate_xmin(m_pl)
m_pl$setXmin(est_xmin$xmin)

est_pars <- estimate_pars(m_pl)
gamma_hat <- est_pars$pars

cat("Оценка gamma (ММП) на хвосте:", gamma_hat, "\n")
cat("Оценённый xmin:", est_xmin$xmin, "\n\n")

# ---------------------- 4) Кластеризация как функция степени C(k) ----------------------

clust_local <- transitivity(g, type = "local")
Ck <- tapply(clust_local, deg, mean, na.rm = TRUE)

plot(as.numeric(names(Ck)), Ck, log = "xy",
     xlab = "k (степень)",
     ylab = "C(k) (средняя локальная кластеризация)",
     main = "Зависимость C(k)")

# ---------------------- 5) Сообщества (Louvain) и атрибуты Web/ML ----------------------
# Чтобы получить воспроизводимый один пример разбиения, фиксируем seed.

set.seed(123)
comm <- cluster_louvain(g)

cat("\nМодульность Louvain:", modularity(comm), "\n")
cat("Число сообществ:", length(sizes(comm)), "\n")
cat("Размер крупнейшего сообщества:", max(sizes(comm)), "\n\n")

# Загрузка атрибутов (0 = web, 1 = ML).
targets <- read.csv("musae_git_target.csv", stringsAsFactors = FALSE)

# Сопоставление ml_target к вершинам по id.
V(g)$ml_target <- targets$ml_target[match(V(g)$name, targets$id)]

# Проверим, сколько вершин успешно сопоставилось.
cat("Сопоставлено ml_target (не NA):", sum(!is.na(V(g)$ml_target)), "из", vcount(g), "\n\n")

# Один агрегированный график: покажем топ-5 крупнейших сообществ по доле типов.
mem <- membership(comm)
tab <- table(mem, V(g)$ml_target)

# Найдём 5 крупнейших сообществ
comm_sizes <- rowSums(tab)
top_idx <- order(comm_sizes, decreasing = TRUE)[1:min(5, length(comm_sizes))]

tab_top <- tab[top_idx, , drop = FALSE]
barplot(t(tab_top),
        beside = TRUE,
        legend.text = c("Web (0)", "ML (1)"),
        main = "Состав типов разработчиков в топ-5 сообществах (Louvain)",
        xlab = "Сообщество (по размеру)",
        ylab = "Число участников")


# ---------------------- 6) Визуализация хабов (необязательная иллюстрация) ----------------------
# Наглядный пример структуры вокруг узлов с высокой степенью.

top_n <- 30
top_nodes <- order(deg, decreasing = TRUE)[1:top_n]
sub_hubs <- induced_subgraph(g, vids = top_nodes)

plot(sub_hubs,
     vertex.label = V(sub_hubs)$name,
     vertex.size = 12,
     vertex.label.cex = 0.6,
     main = paste(top_n, "хабов по степени"))

# ---------------------- 7) Надёжность (robustness) ----------------------
# Оцениваем долю вершин в крупнейшей компоненте:
# 1) от оставшихся вершин (после удаления k узлов)
# 2) от исходного числа вершин N0 (удобно для сравнения с отчётами/статьями)

N0 <- vcount(g)

largest_comp_stats <- function(h, N0) {
  comps <- components(h)
  gcc <- max(comps$csize)
  c(gcc_of_remaining = gcc / vcount(h),
    gcc_of_original  = gcc / N0)
}

set.seed(123)
removals <- c(0, 10, 50, 100, 500, 1000, 2500, 5000, 10000)
removals <- removals[removals <= vcount(g)]

# Случайные удаления
gcc_random_rem  <- numeric(length(removals))
gcc_random_orig <- numeric(length(removals))

for (i in seq_along(removals)) {
  k <- removals[i]
  g_tmp <- if (k == 0) g else delete_vertices(g, sample(V(g), k))
  s <- largest_comp_stats(g_tmp, N0)
  gcc_random_rem[i]  <- s["gcc_of_remaining"]
  gcc_random_orig[i] <- s["gcc_of_original"]
}

# Удаление хабов (фиксированный список по степени исходного графа)
deg <- degree(g)
order_by_deg <- order(deg, decreasing = TRUE)

gcc_target_rem  <- numeric(length(removals))
gcc_target_orig <- numeric(length(removals))

for (i in seq_along(removals)) {
  k <- removals[i]
  g_tmp <- if (k == 0) g else delete_vertices(g, order_by_deg[1:k])
  s <- largest_comp_stats(g_tmp, N0)
  gcc_target_rem[i]  <- s["gcc_of_remaining"]
  gcc_target_orig[i] <- s["gcc_of_original"]
}

# Таблица результатов
robust_df <- data.frame(
  removed = removals,
  gcc_random_remaining = round(gcc_random_rem, 4),
  gcc_random_original  = round(gcc_random_orig, 4),
  gcc_hubs_remaining   = round(gcc_target_rem, 4),
  gcc_hubs_original    = round(gcc_target_orig, 4)
)
print(robust_df)

ylim_all <- range(c(gcc_random_orig, gcc_target_orig))
plot(removals, gcc_random_orig, type = "b", ylim = ylim_all,
     xlab = "Число удалённых узлов",
     ylab = "Доля вершин в GCC (от исходного N)",
     main = "Robustness: случайные удаления vs атака на хабы (GCC от исходных)")

lines(removals, gcc_target_orig, type = "b", lty = 2)
legend("topright",
       legend = c("Случайное удаление", "Удаление хабов"),
       lty = c(1, 2))


# ---------------------- 8) Пересвязывание рёбер при сохранении степеней ----------------------
# тест если перемешать связи, сохранив степени вершин, меняется ли размер GCC

niter_values <- c(0, 1000, 10000, 100000, 200000)
rewire_res <- data.frame(niter = niter_values, gcc_fraction = NA_real_)

for (i in seq_along(niter_values)) {
  niter <- niter_values[i]
  g_rw <- rewire(g, with = keeping_degseq(niter = niter))
  rewire_res$gcc_fraction[i] <- largest_comp_fraction(g_rw)
}

plot(rewire_res$niter, rewire_res$gcc_fraction, type = "b",
     xlab = "n_iter (итерации пересвязывания)",
     ylab = "Доля в крупнейшей компоненте",
     main = "GCC при пересвязывании (keeping_degseq)")

cat("\nRewiring results:\n")
print(rewire_res)
