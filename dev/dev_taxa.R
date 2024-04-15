# [FUNCTIONS] --------------------------------------------------------------
# - fun_taxa_hclust ---------------------------------------------------------
fun_taxa_hclust <- function(
    mtx_similarity
    , int_levels = 7
    , chr_levels = NULL
    , dbl_height = NULL
    , chr_method = c(
      'complete',
      'centroid',
      'single',
      'average',
      'median',
      'ward.D',
      'ward.D2',
      'mcquitty'
    )
){

  # arguments validation
  stopifnot(
    "'mtx_similarity' must be a square matrix of similarity scores between 0 and 1." =
      all(
        is.matrix(mtx_similarity),
        is.numeric(mtx_similarity),
        all(diag(mtx_similarity) == 1),
        all(mtx_similarity >= 0),
        all(mtx_similarity <= 1)
      )
  )

  stopifnot(
    "'int_levels' must be an integer greater between 1 and the number of rows in 'mtx_similarity', indicating how many taxa to extract." =
      all(
        is.numeric(int_levels),
        int_levels <= nrow(mtx_similarity),
        int_levels >= 1
      )
  )

  stopifnot(
    "'dbl_height' must be either NULL or a numeric vector that can be feature scaled to indicate the percentage of maximum height on which to cut the hierarchical tree." =
      any(
        is.null(dbl_height),
        is.numeric(dbl_height)
      )
  )

  stopifnot(
    "'chr_levels' must be either NULL or a character vector with the same length as the number of taxa to extract." =
      any(
        is.null(chr_levels),
        all(
          is.character(chr_levels),
          any(
            length(chr_levels) == length(dbl_height),
            length(chr_levels) == as.integer(int_levels[[1]])
          )
        )
      )
  )

  # data wrangling
  int_levels[[1]] -> int_levels

  as.integer(int_levels) -> int_levels

  chr_method[[1]] -> chr_method

  unique(dbl_height) -> dbl_height

  # hierarchical clustering
  as.dist(
    1 - mtx_similarity
  ) -> dist_distance

  rm(mtx_similarity)

  hclust(
    d = dist_distance,
    method = chr_method
  ) -> hclust_taxa

  rm(dist_distance)

  # heights to cut the tree
  if(!length(dbl_height)){

    seq(
      1, 0, length.out = int_levels
      # 0, 1, length.out = int_levels
    ) -> dbl_height

  }

  if(length(dbl_height) == 1){

    1 -> dbl_height

  }

  # feature scale heights
  if(length(dbl_height) > 1){

    as.numeric(scale(
      dbl_height,
      center = min(dbl_height),
      scale = diff(range(dbl_height))
    )) -> dbl_height

  }

  # breaks as percentage height range
  dbl_height * (
    max(hclust_taxa$height) -
      min(hclust_taxa$height)
  ) -> dbl_height

  # sort breaks
  sort(
    dbl_height,
    decreasing = T
    # decreasing = F
  ) -> dbl_height

  rm(int_levels)

  cutree(
    tree = hclust_taxa,
    h = dbl_height
  ) -> mtx_taxa

  as.matrix(
    mtx_taxa
  ) -> mtx_taxa

  rm(hclust_taxa)

  # taxa labels
  if(!length(chr_levels)){

    paste0(
      'taxon'
      , 1:length(dbl_height)
    ) -> chr_levels

  }

  rm(dbl_height)

  chr_levels ->
    colnames(mtx_taxa)

  rm(chr_levels)

  # taxonomic data frame
  suppressWarnings(
    mtx_taxa %>%
      as_tibble(
        rownames = 'set'
      ) %>%
      pivot_longer(
        cols = -1,
        names_to = 'taxon',
        values_to = 'taxon_id'
      ) %>%
      pivot_wider(
        names_from = 'taxon',
        values_from = 'set'
      ) %>%
      pivot_longer(
        cols = -1,
        names_to = 'taxon',
        values_to = 'set'
      ) %>%
      mutate(
        taxon = factor(
          taxon
          , levels =
            colnames(mtx_taxa)
        )
        , n = map_dbl(
          set,
          length
        )
      ) %>%
      filter(
        n > 0
      ) %>%
      relocate(
        taxon
      ) %>%
      arrange(
        taxon,
        desc(n)
      ) -> df_taxa
  )

  rm(mtx_taxa)

  # data frame subclass
  new_data_frame(
    df_taxa
    , class = c(
      class(df_taxa)
      , 'df_taxa'
    )
  ) -> df_taxa

  # output
  return(df_taxa)

}

# - fun_taxa_desc ---------------------------------------------------------
fun_taxa_desc <- function(df_taxa){

  # arguments validation
  stopifnot(
    "'df_taxa' must be a data frame with the 'df_taxa' subclass." =
      all(
        is.data.frame(df_taxa),
        any(class(df_taxa) == 'df_taxa')
      )
  )

  # output
  return(
    df_taxa %>%
      group_by(
        taxon
      ) %>%
      reframe(
        sets = max(taxon_id),
        n_min = min(n),
        n_mean = mean(n),
        n_max = max(n)
      )
  )

}

# - fun_taxa_unnest ---------------------------------------------------------
fun_taxa_unnest <- function(df_taxa){

  # arguments validation
  stopifnot(
    "'df_taxa' must be a data frame with the 'df_taxa' subclass." =
      all(
        is.data.frame(df_taxa),
        any(class(df_taxa) == 'df_taxa')
      )
  )

  # output
  return(
    df_taxa %>%
      unnest(set) %>%
      pivot_wider(
        id_cols = 'set',
        names_from = 'taxon',
        values_from = 'taxon_id'
      )
  )

}

# - fun_taxa_list ---------------------------------------------------------
fun_taxa_list <- function(df_taxa, lgc_unnest = F){

  # arguments validation
  stopifnot(
    "'df_taxa' must be a data frame with the 'df_taxa' subclass." =
      all(
        is.data.frame(df_taxa),
        any(class(df_taxa) == 'df_taxa')
      )
  )

  stopifnot(
    "'lgc_unnest' must be either TRUE or FALSE." =
      all(
        is.logical(lgc_unnest),
        !is.na(lgc_unnest)
      )
  )

  # split taxonomic data frame
  df_taxa %>%
    split(.$taxon) ->
    list_taxa

  rm(df_taxa)

  # split taxonomic list
  if(lgc_unnest){

    list_taxa %>%
      map(
        ~ .x %>%
          group_by(
            taxon_id
          ) %>%
          unnest(set) %>%
          ungroup() %>%
          split(.$taxon_id) %>%
          set_names(
            map(
              .x = .
              , ~ paste0(
                first(.x$taxon),
                first(.x$taxon_id)
              )
            )
          )
      ) -> list_taxa

    rm(lgc_unnest)

  }

  # output
  return(list_taxa)

}

# # - [dev] fun_taxa_filter -------------------------------------------------------
# list_taxa %>%
#   map(
#     ~ .x %>%
#       filter(
#         map_lgl(
#           set,
#           ~ chr_filter %in% .x
#         )
#       ) %>%
#       unnest(set)
#   )
