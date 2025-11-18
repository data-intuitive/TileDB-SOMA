#' Write an AnnData object to a SOMA
#'
#' Convert an AnnData object to its SOMA counterpart and save it.
#'
#' @param x An `AnnData` object from the `anndataR` package.
#' @param uri URI for the resulting SOMA object.
#' @inheritParams write_soma_objects
#' 
#' @return The URI to the resulting [`SOMAExperiment`] generated from
#' the data contained in `x`.
#'
#' @section Writing AnnData objects:
#' AnnData objects are written out as [`SOMAExperiment`] objects:
#' \itemize{
#'  \item `obs` is written out as a [`SOMADataFrame`] called "obs" at
#'   the experiment level.
#'  \item `var` is written out as a [`SOMADataFrame`] called "var" within
#'   the measurement.
#'  \item `X` matrix is written out as a [`SOMASparseNDArray`] called
#'   "X" within the measurement's X group.
#'  \item `layers` are written out as [`SOMASparseNDArray`] objects within 
#'   the measurement's X group.
#'  \item `obsm` matrices are written out as [`SOMASparseNDArray`] objects
#'   within the measurement's obsm group.
#'  \item `varm` matrices are written out as [`SOMASparseNDArray`] objects
#'   within the measurement's varm group.
#'  \item `obsp` matrices are written out as [`SOMASparseNDArray`] objects
#'   within the measurement's obsp group.
#'  \item `varp` matrices are written out as [`SOMASparseNDArray`] objects
#'   within the measurement's varp group.
#' }
#' Expression matrices are transposed (cells as rows) prior to writing.
#'
#' @method write_soma AbstractAnnData
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE) && requireNamespace("anndataR", quietly = TRUE)
#' \donttest{
#' uri <- withr::local_tempfile(pattern = "anndata")
#' 
#' # Create a simple AnnData object
#' adata <- anndataR::AnnData(
#'   X = matrix(1:15, nrow = 3),
#'   obs = data.frame(cell_type = c("A", "B", "C")),
#'   var = data.frame(gene_name = paste0("gene_", 1:5))
#' )
#' 
#' uri <- write_soma(adata, uri)
#' 
#' (exp <- SOMAExperimentOpen(uri))
#' exp$obs
#' (ms <- exp$ms$get("RNA"))
#' ms$var
#' ms$X$names()
#' 
#' exp$close()
#' }
#'
write_soma.AbstractAnnData <- function(
  x,
  uri,
  ms_name = "RNA",
  ...,
  ingest_mode = "write",
  platform_config = NULL,
  tiledbsoma_ctx = NULL
) {
  check_package("anndataR")
  stopifnot(
    "'uri' must be a single character value" = is.null(uri) ||
      is_scalar_character(uri),
    "'ms_name' must be a single character value" = is_scalar_character(ms_name) &&
      nzchar(ms_name),
    "'x' must be an AbstractAnnData object" = inherits(x, "AbstractAnnData")
  )
  ingest_mode <- match.arg(arg = ingest_mode, choices = c("write", "resume"))
  
  # Create the experiment
  experiment <- SOMAExperimentCreate(
    uri = uri,
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  on.exit(experiment$close(), add = TRUE, after = FALSE)
  
  # Write cell-level metadata (obs)
  spdl::info("Adding obs")
  obs_df <- x$obs
  if (is.null(rownames(obs_df))) {
    rownames(obs_df) <- paste0("obs", seq_len(nrow(obs_df)))
  }
  obs_df <- .df_index(obs_df, axis = "obs")
  obs_df[[attr(obs_df, "index")]] <- rownames(x$obs)
  
  write_soma(
    x = obs_df,
    uri = "obs",
    soma_parent = experiment,
    key = "obs",
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  
  # Create measurements collection
  spdl::info("Creating measurements collection")
  expms <- SOMACollectionCreate(
    file_path(experiment$uri, "ms"),
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  withCallingHandlers(
    expr = .register_soma_object(expms, soma_parent = experiment, key = "ms"),
    existingKeyWarning = .maybe_muffle
  )
  
  # Create measurement
  ms_uri <- .check_soma_uri(uri = ms_name, soma_parent = expms)
  ms <- SOMAMeasurementCreate(
    uri = ms_uri,
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  on.exit(ms$close(), add = TRUE, after = FALSE)
  
  # Write feature-level metadata (var)
  spdl::info("Adding var")
  var_df <- x$var
  if (is.null(rownames(var_df))) {
    rownames(var_df) <- paste0("var", seq_len(nrow(var_df)))
  }
  var_df <- .df_index(var_df, axis = "var")
  var_df[[attr(var_df, "index")]] <- rownames(x$var)
  
  write_soma(
    x = var_df,
    uri = "var",
    soma_parent = ms,
    key = "var",
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  
  # Create X collection
  X <- if (!"X" %in% ms$names()) {
    SOMACollectionCreate(
      file_path(ms$uri, "X"),
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    )
  } else {
    SOMACollectionOpen(file_path(ms$uri, "X"), mode = "WRITE")
  }
  withCallingHandlers(
    .register_soma_object(X, soma_parent = ms, key = "X"),
    existingKeyWarning = .maybe_muffle
  )
  on.exit(X$close(), add = TRUE, after = FALSE)
  
  # Write main X matrix
  if (!is.null(x$X)) {
    spdl::info("Adding X matrix")
    write_soma(
      x = x$X,
      uri = "X",
      soma_parent = X,
      sparse = TRUE,
      transpose = TRUE,  # AnnData uses genes x cells, SOMA uses cells x genes
      key = "X",
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    )
  }
  
  # Write layers
  if (length(x$layers_keys()) > 0) {
    spdl::info("Adding layers")
    for (layer_name in x$layers_keys()) {
      spdl::info("Adding layer {}", layer_name)
      layer_data <- x$layers[[layer_name]]
      write_soma(
        x = layer_data,
        uri = layer_name,
        soma_parent = X,
        sparse = TRUE,
        transpose = TRUE,  # AnnData uses genes x cells, SOMA uses cells x genes
        key = layer_name,
        ingest_mode = ingest_mode,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      )
    }
  }
  
  # Write obsm matrices
  if (length(x$obsm_keys()) > 0) {
    spdl::info("Adding obsm matrices")
    obsm <- if (!"obsm" %in% ms$names()) {
      SOMACollectionCreate(
        file_path(ms$uri, "obsm"),
        ingest_mode = ingest_mode,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      )
    } else {
      SOMACollectionOpen(file_path(ms$uri, "obsm"), mode = "WRITE")
    }
    withCallingHandlers(
      .register_soma_object(obsm, soma_parent = ms, key = "obsm"),
      existingKeyWarning = .maybe_muffle
    )
    on.exit(obsm$close(), add = TRUE, after = FALSE)
    
    for (obsm_name in x$obsm_keys()) {
      spdl::info("Adding obsm matrix {}", obsm_name)
      obsm_data <- x$obsm[[obsm_name]]
      write_soma(
        x = obsm_data,
        uri = obsm_name,
        soma_parent = obsm,
        sparse = TRUE,
        transpose = FALSE,  # obsm should keep original orientation
        key = obsm_name,
        ingest_mode = ingest_mode,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      )
    }
  }
  
  # Write varm matrices
  if (length(x$varm_keys()) > 0) {
    spdl::info("Adding varm matrices")
    varm <- if (!"varm" %in% ms$names()) {
      SOMACollectionCreate(
        file_path(ms$uri, "varm"),
        ingest_mode = ingest_mode,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      )
    } else {
      SOMACollectionOpen(file_path(ms$uri, "varm"), mode = "WRITE")
    }
    withCallingHandlers(
      .register_soma_object(varm, soma_parent = ms, key = "varm"),
      existingKeyWarning = .maybe_muffle
    )
    on.exit(varm$close(), add = TRUE, after = FALSE)
    
    for (varm_name in x$varm_keys()) {
      spdl::info("Adding varm matrix {}", varm_name)
      varm_data <- x$varm[[varm_name]]
      write_soma(
        x = varm_data,
        uri = varm_name,
        soma_parent = varm,
        sparse = TRUE,
        transpose = FALSE,  # varm should keep original orientation
        key = varm_name,
        ingest_mode = ingest_mode,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      )
    }
  }
  
  # Write obsp matrices
  if (length(x$obsp_keys()) > 0) {
    spdl::info("Adding obsp matrices")
    obsp <- if (!"obsp" %in% ms$names()) {
      SOMACollectionCreate(
        file_path(ms$uri, "obsp"),
        ingest_mode = ingest_mode,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      )
    } else {
      SOMACollectionOpen(file_path(ms$uri, "obsp"), mode = "WRITE")
    }
    withCallingHandlers(
      .register_soma_object(obsp, soma_parent = ms, key = "obsp"),
      existingKeyWarning = .maybe_muffle
    )
    on.exit(obsp$close(), add = TRUE, after = FALSE)
    
    for (obsp_name in x$obsp_keys()) {
      spdl::info("Adding obsp matrix {}", obsp_name)
      obsp_data <- x$obsp[[obsp_name]]
      write_soma(
        x = obsp_data,
        uri = obsp_name,
        soma_parent = obsp,
        sparse = TRUE,
        transpose = FALSE,  # obsp should keep original orientation
        key = obsp_name,
        ingest_mode = ingest_mode,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      )
    }
  }
  
  # Write varp matrices
  if (length(x$varp_keys()) > 0) {
    spdl::info("Adding varp matrices")
    varp <- if (!"varp" %in% ms$names()) {
      SOMACollectionCreate(
        file_path(ms$uri, "varp"),
        ingest_mode = ingest_mode,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      )
    } else {
      SOMACollectionOpen(file_path(ms$uri, "varp"), mode = "WRITE")
    }
    withCallingHandlers(
      .register_soma_object(varp, soma_parent = ms, key = "varp"),
      existingKeyWarning = .maybe_muffle
    )
    on.exit(varp$close(), add = TRUE, after = FALSE)
    
    for (varp_name in x$varp_keys()) {
      spdl::info("Adding varp matrix {}", varp_name)
      varp_data <- x$varp[[varp_name]]
      write_soma(
        x = varp_data,
        uri = varp_name,
        soma_parent = varp,
        sparse = TRUE,
        transpose = FALSE,  # varp should keep original orientation
        key = varp_name,
        ingest_mode = ingest_mode,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      )
    }
  }
  
  withCallingHandlers(
    .register_soma_object(ms, soma_parent = expms, key = ms_name),
    existingKeyWarning = .maybe_muffle
  )
  
  return(experiment$uri)
}

