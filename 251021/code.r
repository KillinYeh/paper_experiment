library(ggplot2)
library(dplyr)
library(scales)
library(ggrepel)
plot_utilization_distribution <- function(x,
                                              title_prefix,
                                              show_rug = TRUE,
                                              label_digits = 3) {
    # --- 清理：只移除非有限值，不做任何縮放或夾值 ---
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) stop("沒有可用的數值。")
    
    # 邊界檢查（只提醒，不修改）
    n_lt0 <- sum(x < 0, na.rm = TRUE)
    n_gt1 <- sum(x > 1, na.rm = TRUE)
    if (n_lt0 + n_gt1 > 0) {
        message(sprintf("[提醒] 有 %d 筆 <0、%d 筆 >1 的值（未被修改）。", n_lt0, n_gt1))
    }
    
    # --- deciles（每10%分位） ---
    probs_dec <- seq(0.1, 0.9, 0.1)
    qs_dec    <- as.numeric(quantile(x, probs_dec, na.rm = TRUE, names = FALSE))
    df_dec <- data.frame(
        p   = probs_dec,
        q   = qs_dec,
        lab = paste0(percent(probs_dec), "\n", round(qs_dec, label_digits))
    )
    
    # --- 密度（在資料範圍內估計） ---
    xr <- range(x, na.rm = TRUE)
    pad <- diff(xr) * 0.02
    den <- tryCatch(density(x, from = xr[1], to = xr[2], na.rm = TRUE),
                    error = function(e) NULL)
    if (is.null(den) || !any(is.finite(den$y))) {
        # 後備：用直方圖近似密度
        brks <- pretty(xr, n = 30)
        h <- hist(x, breaks = brks, plot = FALSE)
        df_den <- data.frame(x = h$mids, y = h$density)
        f <- approxfun(df_den$x, df_den$y, rule = 2)
    } else {
        df_den <- data.frame(x = den$x, y = den$y)
        f <- approxfun(df_den$x, df_den$y, rule = 2)
    }
    
    df_dec$y <- f(df_dec$q)
    ymax <- max(df_den$y, na.rm = TRUE); if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    # --- 圖1：密度 + deciles（標籤避讓） ---
    p_density <- ggplot(df_den, aes(x, y)) +
        geom_area(alpha = 0.2) +
        geom_line(linewidth = 1) +
        geom_segment(data = df_dec,
                     aes(x = q, xend = q, y = 0, yend = y),
                     inherit.aes = FALSE, linewidth = 0.7) +
        geom_label_repel(data = df_dec, aes(x = q, y = y, label = lab),
                         inherit.aes = FALSE, size = 3,
                         direction = "y", box.padding = 0.2, point.padding = 0.2,
                         min.segment.length = 0, seed = 123) +
        { if (show_rug) geom_rug(data = data.frame(x = x), aes(x = x, y = NULL),
                                 sides = "b", alpha = 0.25) } +
        coord_cartesian(xlim = xr + c(-pad, pad), expand = FALSE, clip = "off") +
        labs(title = paste0(title_prefix, " 分布（含每10%百分點）"),
             x = paste0(title_prefix, " (實際範圍)"), y = "密度") +
        theme_minimal()
    
    # --- 圖2：ECDF + deciles（標籤避讓） ---
    # --- 圖2：ECDF + deciles（把標籤放到右側邊界外，不會擋線） ---
    x_sorted <- sort(x)
    Fhat <- ecdf(x)
    df_ecdf <- data.frame(x = x_sorted, F = Fhat(x_sorted))
    
    xr  <- range(x, na.rm = TRUE)
    pad <- max(1e-8, diff(xr) * 0.02)          # 內側少許邊距
    x_lab <- xr[2] + pad * 0.8                 # 右側放標籤的位置（在繪圖區外）
    
    df_labels <- data.frame(
        x0 = df_dec$q,            # decile 的 x
        y  = df_dec$p,            # decile 的 F(x)
        x1 = x_lab,               # 標籤 x 位置（右側）
        lab = df_dec$lab
    )
    
    p_ecdf <- ggplot(df_ecdf, aes(x, F)) +
        geom_step() +
        geom_point(data = data.frame(x = df_dec$q, F = df_dec$p), size = 2) +
        # 從 decile 點畫水平線到右側標籤
        geom_segment(data = df_labels,
                     aes(x = x0, xend = x1, y = y, yend = y),
                     linewidth = 0.4, alpha = 0.7, inherit.aes = FALSE) +
        # 右側標籤（完全不與線重疊）
        geom_label(data = df_labels,
                   aes(x = x1, y = y, label = lab),
                   inherit.aes = FALSE, size = 3,
                   label.size = 0.2, label.padding = unit(0.12, "lines")) +
        # 擴大右側視窗讓標籤有空間；clip="off" 讓標籤能畫在外側
        coord_cartesian(xlim = xr + c(-pad, pad * 1.8), ylim = c(0, 1), clip = "off") +
        labs(title = paste0(title_prefix, " 的 ECDF（每10%分位）"),
             x = paste0(title_prefix, " (實際範圍)"), y = "F(x)") +
        theme_minimal() +
        theme(plot.margin = margin(5.5, 40, 5.5, 5.5, "pt"))  # 右邊多留一點空白
    
    print(p_ecdf)
    # --- 圖2：ECDF + deciles（把標籤放到右側邊界外，不會擋線） ---
    x_sorted <- sort(x)
    Fhat <- ecdf(x)
    df_ecdf <- data.frame(x = x_sorted, F = Fhat(x_sorted))
    
    xr  <- range(x, na.rm = TRUE)
    pad <- max(1e-8, diff(xr) * 0.02)          # 內側少許邊距
    x_lab <- xr[2] + pad * 0.8                 # 右側放標籤的位置（在繪圖區外）
    
    df_labels <- data.frame(
        x0 = df_dec$q,            # decile 的 x
        y  = df_dec$p,            # decile 的 F(x)
        x1 = x_lab,               # 標籤 x 位置（右側）
        lab = df_dec$lab
    )
    
    p_ecdf <- ggplot(df_ecdf, aes(x, F)) +
        geom_step() +
        geom_point(data = data.frame(x = df_dec$q, F = df_dec$p), size = 2) +
        # 從 decile 點畫水平線到右側標籤
        geom_segment(data = df_labels,
                     aes(x = x0, xend = x1, y = y, yend = y),
                     linewidth = 0.4, alpha = 0.7, inherit.aes = FALSE) +
        # 右側標籤（完全不與線重疊）
        geom_label(data = df_labels,
                   aes(x = x1, y = y, label = lab),
                   inherit.aes = FALSE, size = 3,
                   label.size = 0.2, label.padding = unit(0.12, "lines")) +
        # 擴大右側視窗讓標籤有空間；clip="off" 讓標籤能畫在外側
        coord_cartesian(xlim = xr + c(-pad, pad * 1.8), ylim = c(0, 1), clip = "off") +
        labs(title = paste0(title_prefix, "  CDF"),
             x = paste0(title_prefix), y = "F(x)") +
        theme_minimal() +
        theme(plot.margin = margin(5.5, 40, 5.5, 5.5, "pt"))  # 右邊多留一點空白
    
    print(p_ecdf)
    
    
    invisible(list(
        p_density = p_density,
        p_ecdf    = p_ecdf,
        deciles   = df_dec[, c("p", "q")]
    ))
}
--------------------experiment------------------------

plot_utilization_distribution(utilization_email$utilization,"utilization_email")
plot_utilization_distribution(utilization_fetal$utilization,"utilization_fetal")
plot_utilization_distribution(utilization_gene$utilization,"utilization_gene")
plot_utilization_distribution(utilization_mush$utilization,"utilization_mush")
plot_utilization_distribution(utilization_loan$utilization,"utilization_loan")
plot_utilization_distribution(utilization_gesture$utilization,"utilization_gesture")

