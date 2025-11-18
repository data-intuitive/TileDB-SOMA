test_that("to_anndata with obs_index parameter", {
  skip_if(!extended_tests())
  skip_if_not_installed("anndataR")
  
  # Create test data with string cell names
  # AnnData expects cells x genes format
  X <- matrix(rpois(50, 3), nrow = 5, ncol = 10, 
              dimnames = list(paste0("cell", 1:5), paste0("gene", 1:10)))
  
  obs <- data.frame(
    cell_id = paste0("cell_", 1:5),
    cell_type = c("A", "B", "C", "A", "B"),
    n_genes = c(100, 150, 120, 80, 110),
    row.names = rownames(X)
  )
  
  var <- data.frame(
    gene_name = colnames(X),
    highly_variable = c(rep(TRUE, 5), rep(FALSE, 5)),
    row.names = colnames(X)
  )
  
  # Create and write AnnData object
  adata <- anndataR::AnnData(
    X = Matrix::Matrix(X, sparse = TRUE),
    obs = obs,
    var = var
  )
  
  uri <- tempfile(pattern = "to-anndata-obs-index-test")
  experiment_uri <- write_soma(adata, uri)
  exp <- SOMAExperimentOpen(experiment_uri)
  
  # Create a query object
  query <- SOMAExperimentAxisQuery$new(exp, "RNA")
  
  # Test with obs_index parameter
  expect_no_condition(adata_result <- query$to_anndata(obs_index = "cell_id"))
  expect_s3_class(adata_result, "AbstractAnnData")
  
  # Verify that the obs index was used
  obs_result <- as.data.frame(adata_result$obs)
  expect_true(all(paste0("cell_", 1:5) %in% rownames(obs_result)))
  
  # Clean up
  exp$close()
})
