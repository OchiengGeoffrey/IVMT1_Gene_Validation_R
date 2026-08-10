# Identify the 41 pair_ids with HSP support in both flanks

gap_start <- 649L
gap_end <- 732L

n_flank_min <- 550L
n_flank_max <- 648L

c_flank_min <- 733L
c_flank_max <- 830L

both_flank_pairs <- character(0)

pair_summary <- lapply(split(tac, tac$pair_id), function(sub_df) {

  has_n <- FALSE
  has_c <- FALSE

  for (i in seq_len(nrow(sub_df))) {

    starts <- as.integer(
      base::strsplit(sub_df$qstart_union[i], ";", fixed = TRUE)[[1]]
    )

    ends <- as.integer(
      base::strsplit(sub_df$qend_union[i], ";", fixed = TRUE)[[1]]
    )

    if (any(starts <= n_flank_max & ends >= n_flank_min))
      has_n <- TRUE

    if (any(starts <= c_flank_max & ends >= c_flank_min))
      has_c <- TRUE
  }

  data.frame(
    pair_id = unique(sub_df$pair_id),
    has_N_flank = has_n,
    has_C_flank = has_c,
    stringsAsFactors = FALSE
  )
})

pair_summary <- do.call(rbind, pair_summary)

both_flank <- pair_summary[
  pair_summary$has_N_flank & pair_summary$has_C_flank,
  ,
  drop = FALSE
]

nrow(both_flank)