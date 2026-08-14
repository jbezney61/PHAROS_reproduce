#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(monocle3)
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(Matrix)
})

input_file <- "GSM7056151_sciPlex_4_preprocessed_cds.list.rds"
output_dir <- "GSM7056151_sciPlex_4_export"

dir.create(
    output_dir,
    showWarnings = FALSE,
    recursive = TRUE
)

cat("Reading:", input_file, "\n")

cds_list <- readRDS(input_file)

if (!is.list(cds_list) || length(cds_list) == 0) {
    stop("The RDS does not contain a non-empty list.")
}

if (is.null(names(cds_list)) || any(names(cds_list) == "")) {
    names(cds_list) <- paste0("cds_", seq_along(cds_list))
}

cat("Objects found:\n")

for (object_name in names(cds_list)) {
    cds <- cds_list[[object_name]]

    cat(
        sprintf(
            "  %s: %s cells x %s genes\n",
            object_name,
            format(ncol(cds), big.mark = ","),
            format(nrow(cds), big.mark = ",")
        )
    )
}


# ---------------------------------------------------------------------
# Validate objects
# ---------------------------------------------------------------------
valid_objects <- vapply(
    cds_list,
    function(x) inherits(x, "cell_data_set"),
    logical(1)
)

if (!all(valid_objects)) {
    stop(
        "These elements are not Monocle cell_data_set objects: ",
        paste(names(cds_list)[!valid_objects], collapse = ", ")
    )
}

has_counts <- vapply(
    cds_list,
    function(x) "counts" %in% assayNames(x),
    logical(1)
)

if (!all(has_counts)) {
    stop(
        "These objects do not contain a counts assay: ",
        paste(names(cds_list)[!has_counts], collapse = ", ")
    )
}


# ---------------------------------------------------------------------
# Validate gene identifiers and ordering
# ---------------------------------------------------------------------
reference_gene_ids <- rownames(cds_list[[1]])

if (is.null(reference_gene_ids)) {
    stop("The first CDS object has no gene row names.")
}

same_gene_order <- vapply(
    cds_list,
    function(x) identical(rownames(x), reference_gene_ids),
    logical(1)
)

if (!all(same_gene_order)) {
    stop(
        "Gene identifiers or gene order differ among objects: ",
        paste(names(cds_list)[!same_gene_order], collapse = ", ")
    )
}

if (anyDuplicated(reference_gene_ids)) {
    duplicate_ids <- unique(
        reference_gene_ids[duplicated(reference_gene_ids)]
    )

    stop(
        "Duplicate gene identifiers detected. Examples: ",
        paste(head(duplicate_ids, 10), collapse = ", ")
    )
}

cat("All objects have identical genes and gene order.\n")


# ---------------------------------------------------------------------
# Validate globally unique cell identifiers
# ---------------------------------------------------------------------
all_cell_ids <- unlist(
    lapply(cds_list, colnames),
    use.names = FALSE
)

if (anyNA(all_cell_ids)) {
    stop("Missing cell identifiers were detected.")
}

if (anyDuplicated(all_cell_ids)) {
    duplicate_ids <- unique(
        all_cell_ids[duplicated(all_cell_ids)]
    )

    stop(
        "Duplicate cell identifiers detected across objects. Examples: ",
        paste(head(duplicate_ids, 10), collapse = ", ")
    )
}

cat(
    "Validated",
    format(length(all_cell_ids), big.mark = ","),
    "unique cells.\n"
)


# ---------------------------------------------------------------------
# Helper for exporting metadata
# ---------------------------------------------------------------------
clean_metadata <- function(metadata, metadata_type) {

    metadata <- as.data.frame(
        metadata,
        stringsAsFactors = FALSE
    )

    supported <- vapply(
        metadata,
        function(x) {
            is.atomic(x) ||
                is.factor(x) ||
                inherits(x, "Date") ||
                inherits(x, "POSIXct")
        },
        logical(1)
    )

    if (any(!supported)) {
        warning(
            "Dropping unsupported ",
            metadata_type,
            " columns: ",
            paste(names(metadata)[!supported], collapse = ", ")
        )
    }

    metadata <- metadata[, supported, drop = FALSE]

    metadata[] <- lapply(
        metadata,
        function(x) {
            if (is.factor(x)) {
                return(as.character(x))
            }

            if (inherits(x, "Date") ||
                inherits(x, "POSIXct")) {
                return(as.character(x))
            }

            x
        }
    )

    metadata
}


# ---------------------------------------------------------------------
# Export gene metadata once
# ---------------------------------------------------------------------
gene_metadata <- clean_metadata(
    rowData(cds_list[[1]]),
    "gene metadata"
)

if ("gene_id" %in% colnames(gene_metadata)) {
    colnames(gene_metadata)[
        colnames(gene_metadata) == "gene_id"
    ] <- "gene_id_metadata"
}

gene_metadata <- cbind(
    gene_id = reference_gene_ids,
    gene_metadata
)

gene_file <- file.path(
    output_dir,
    "genes.tsv"
)

write.table(
    gene_metadata,
    file = gene_file,
    sep = "\t",
    quote = TRUE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
)

cat("Wrote gene metadata:", gene_file, "\n")


# ---------------------------------------------------------------------
# Export each object separately
# ---------------------------------------------------------------------
manifest_rows <- list()

for (i in seq_along(cds_list)) {

    object_name <- names(cds_list)[i]
    cds <- cds_list[[i]]

    safe_name <- gsub(
        "[^A-Za-z0-9._-]",
        "_",
        object_name
    )

    count_filename <- paste0(
        sprintf("%02d", i),
        "_",
        safe_name,
        "_counts.mtx"
    )

    obs_filename <- paste0(
        sprintf("%02d", i),
        "_",
        safe_name,
        "_obs.tsv"
    )

    count_path <- file.path(
        output_dir,
        count_filename
    )

    obs_path <- file.path(
        output_dir,
        obs_filename
    )

    cat("Exporting:", object_name, "\n")

    # Counts are genes x cells in the Monocle object.
    count_matrix <- assay(
        cds,
        "counts"
    )

    if (!inherits(count_matrix, "sparseMatrix")) {
        count_matrix <- as(
            count_matrix,
            "dgCMatrix"
        )
    } else {
        count_matrix <- as(
            count_matrix,
            "dgCMatrix"
        )
    }

    if (!identical(dim(count_matrix), dim(cds))) {
        stop(
            "Count matrix dimensions do not match the CDS dimensions for ",
            object_name
        )
    }

    Matrix::writeMM(
        count_matrix,
        file = count_path
    )

    cell_metadata <- clean_metadata(
        colData(cds),
        "cell metadata"
    )

    if ("cell_id" %in% colnames(cell_metadata)) {
        colnames(cell_metadata)[
            colnames(cell_metadata) == "cell_id"
        ] <- "cell_id_metadata"
    }

    if ("source_cds" %in% colnames(cell_metadata)) {
        colnames(cell_metadata)[
            colnames(cell_metadata) == "source_cds"
        ] <- "source_cds_metadata"
    }

    cell_metadata <- cbind(
        cell_id = colnames(cds),
        source_cds = object_name,
        cell_metadata
    )

    write.table(
        cell_metadata,
        file = obs_path,
        sep = "\t",
        quote = TRUE,
        row.names = FALSE,
        col.names = TRUE,
        na = ""
    )

    manifest_rows[[i]] <- data.frame(
        source_cds = object_name,
        counts_file = count_filename,
        obs_file = obs_filename,
        n_cells = ncol(cds),
        n_genes = nrow(cds),
        stringsAsFactors = FALSE
    )

    cat(
        "  Counts:",
        count_path,
        "\n"
    )

    cat(
        "  Metadata:",
        obs_path,
        "\n"
    )

    rm(count_matrix)
    gc()
}


# ---------------------------------------------------------------------
# Write manifest
# ---------------------------------------------------------------------
manifest <- do.call(
    rbind,
    manifest_rows
)

manifest_file <- file.path(
    output_dir,
    "manifest.tsv"
)

write.table(
    manifest,
    file = manifest_file,
    sep = "\t",
    quote = TRUE,
    row.names = FALSE,
    col.names = TRUE
)

cat("\nExport complete.\n")
cat("Manifest:", manifest_file, "\n")
cat("Export directory:", normalizePath(output_dir), "\n")