#!/usr/bin/env Rscript

# Export raw single-cell counts and metadata from a Seurat, Monocle3,
# or SingleCellExperiment RDS file into files that can be assembled as h5ad.
#
# Usage:
#   Rscript export_raw_counts_for_h5ad.R \
#       Integrated_Dataset_normalized_seurobj.rds \
#       breast_cancer_export \
#       RNA
#
# Arguments:
#   1. Input RDS file
#   2. Output directory
#   3. Preferred Seurat assay (default: RNA)
#
# Important:
#   - Only a real raw-count assay/layer is exported.
#   - The script never substitutes normalized data, scale.data, or logcounts.
#   - Exported Matrix Market matrices are genes x cells.

suppressPackageStartupMessages({
    library(Matrix)
    library(SingleCellExperiment)
    library(SummarizedExperiment)
})

# Load optional class/method packages when available.
if (requireNamespace("SeuratObject", quietly = TRUE)) {
    suppressPackageStartupMessages(library(SeuratObject))
}
if (requireNamespace("monocle3", quietly = TRUE)) {
    suppressPackageStartupMessages(library(monocle3))
}

args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1) {
    args[[1]]
} else {
    "Integrated_Dataset_normalized_seurobj.rds"
}

output_dir <- if (length(args) >= 2) {
    args[[2]]
} else {
    "breast_cancer_export"
}

preferred_seurat_assay <- if (length(args) >= 3) {
    args[[3]]
} else {
    "RNA"
}

if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Reading:", input_file, "\n")
loaded_object <- readRDS(input_file)


# -----------------------------------------------------------------------------
# Object helpers
# -----------------------------------------------------------------------------
is_supported_object <- function(x) {
    inherits(x, "Seurat") ||
        inherits(x, "SingleCellExperiment") ||
        inherits(x, "cell_data_set")
}

object_type <- function(x) {
    if (inherits(x, "Seurat")) {
        return("Seurat")
    }
    if (inherits(x, "cell_data_set")) {
        return("Monocle3_cell_data_set")
    }
    if (inherits(x, "SingleCellExperiment")) {
        return("SingleCellExperiment")
    }
    paste(class(x), collapse = ";")
}

# Accept either one object or a named list of objects.
if (is_supported_object(loaded_object)) {
    object_list <- list(dataset = loaded_object)
} else if (is.list(loaded_object) && length(loaded_object) > 0) {
    object_list <- loaded_object
} else {
    stop(
        "The RDS must contain a Seurat, Monocle3 cell_data_set, ",
        "SingleCellExperiment, or a non-empty list of these objects."
    )
}

valid_objects <- vapply(object_list, is_supported_object, logical(1))
if (!all(valid_objects)) {
    stop(
        "Unsupported elements in the RDS list: ",
        paste(names(object_list)[!valid_objects], collapse = ", ")
    )
}

if (is.null(names(object_list))) {
    names(object_list) <- rep("", length(object_list))
}
empty_names <- is.na(names(object_list)) | names(object_list) == ""
names(object_list)[empty_names] <- paste0(
    "object_",
    which(empty_names)
)

if (anyDuplicated(names(object_list))) {
    names(object_list) <- make.unique(names(object_list), sep = "_")
}

safe_object_names <- make.unique(
    gsub("[^A-Za-z0-9._-]", "_", names(object_list)),
    sep = "_"
)

cat("Objects found:\n")
for (i in seq_along(object_list)) {
    x <- object_list[[i]]
    cat(
        sprintf(
            "  %s [%s]: %s cells x %s features\n",
            names(object_list)[i],
            object_type(x),
            format(ncol(x), big.mark = ","),
            format(nrow(x), big.mark = ",")
        )
    )
}


# -----------------------------------------------------------------------------
# Metadata helpers
# -----------------------------------------------------------------------------
clean_metadata <- function(metadata, metadata_type) {
    metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

    if (ncol(metadata) == 0) {
        return(metadata)
    }

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
            "Dropping unsupported ", metadata_type, " columns: ",
            paste(names(metadata)[!supported], collapse = ", ")
        )
    }

    metadata <- metadata[, supported, drop = FALSE]

    metadata[] <- lapply(
        metadata,
        function(x) {
            if (is.factor(x) || inherits(x, "Date") || inherits(x, "POSIXct")) {
                return(as.character(x))
            }
            if (inherits(x, "integer64")) {
                return(as.character(x))
            }
            x
        }
    )

    colnames(metadata) <- make.unique(colnames(metadata), sep = "_")
    metadata
}

rename_reserved_columns <- function(metadata, reserved_names, suffix = "_metadata") {
    for (nm in reserved_names) {
        if (nm %in% colnames(metadata)) {
            colnames(metadata)[colnames(metadata) == nm] <- paste0(nm, suffix)
        }
    }
    colnames(metadata) <- make.unique(colnames(metadata), sep = "_")
    metadata
}

reorder_metadata_rows <- function(metadata, expected_ids, metadata_type) {
    metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

    metadata_ids <- rownames(metadata)
    if (!is.null(metadata_ids) && all(expected_ids %in% metadata_ids)) {
        metadata <- metadata[expected_ids, , drop = FALSE]
        rownames(metadata) <- expected_ids
        return(metadata)
    }

    if (nrow(metadata) == 0 && length(expected_ids) > 0) {
        return(data.frame(row.names = expected_ids))
    }

    if (nrow(metadata) != length(expected_ids)) {
        stop(
            metadata_type,
            " has ", nrow(metadata),
            " rows but expected ", length(expected_ids), "."
        )
    }

    rownames(metadata) <- expected_ids
    metadata
}


# -----------------------------------------------------------------------------
# Seurat helpers
# -----------------------------------------------------------------------------
require_seurat_object <- function() {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
        stop(
            "A Seurat object was detected, but the SeuratObject package is not installed. ",
            "Install it with install.packages('SeuratObject')."
        )
    }
}

seurat_assay_names <- function(x) {
    require_seurat_object()
    as.character(SeuratObject::Assays(x))
}

seurat_layer_names <- function(x, assay_name) {
    require_seurat_object()

    layers <- tryCatch(
        SeuratObject::Layers(x, assay = assay_name, search = NA),
        error = function(e) character(0)
    )

    if (length(layers) == 0) {
        assay_object <- x[[assay_name]]
        available_slots <- methods::slotNames(assay_object)
        layers <- intersect(c("counts", "data", "scale.data"), available_slots)
    }

    as.character(layers)
}

raw_count_layer_names <- function(x, assay_name) {
    layers <- seurat_layer_names(x, assay_name)
    layers[grepl("^counts($|[._])", layers)]
}

select_seurat_assay <- function(x, preferred_assay = "RNA") {
    require_seurat_object()

    assays <- seurat_assay_names(x)
    default_assay <- tryCatch(
        SeuratObject::DefaultAssay(x),
        error = function(e) NA_character_
    )

    candidates <- unique(c(preferred_assay, "RNA", default_assay, assays))
    candidates <- candidates[!is.na(candidates) & candidates %in% assays]

    for (candidate in candidates) {
        if (length(raw_count_layer_names(x, candidate)) > 0) {
            if (!identical(candidate, preferred_assay)) {
                warning(
                    "Preferred Seurat assay '", preferred_assay,
                    "' did not contain raw counts. Using assay '", candidate, "'."
                )
            }
            return(candidate)
        }
    }

    details <- vapply(
        assays,
        function(a) {
            paste0(a, "={", paste(seurat_layer_names(x, a), collapse = ","), "}")
        },
        character(1)
    )

    stop(
        "No Seurat assay contains a raw counts layer. Available assays/layers: ",
        paste(details, collapse = "; "),
        ". Raw counts cannot be reconstructed from normalized data."
    )
}

get_seurat_layer <- function(x, assay_name, layer_name) {
    require_seurat_object()
    exports <- getNamespaceExports("SeuratObject")

    if ("LayerData" %in% exports) {
        return(
            SeuratObject::LayerData(
                object = x,
                assay = assay_name,
                layer = layer_name
            )
        )
    }

    # Compatibility fallback for SeuratObject versions before LayerData.
    if (identical(layer_name, "counts") && "GetAssayData" %in% exports) {
        return(
            SeuratObject::GetAssayData(
                object = x,
                assay = assay_name,
                slot = "counts"
            )
        )
    }

    stop(
        "The installed SeuratObject version cannot read layer '",
        layer_name,
        "'. Update SeuratObject to a version that provides LayerData()."
    )
}

extract_seurat_counts <- function(x, assay_name) {
    require_seurat_object()

    count_layers <- raw_count_layer_names(x, assay_name)
    if (length(count_layers) == 0) {
        stop("No counts layer found in Seurat assay: ", assay_name)
    }

    # Prefer an exact counts layer when present.
    if ("counts" %in% count_layers) {
        counts <- get_seurat_layer(
            x = x,
            assay_name = assay_name,
            layer_name = "counts"
        )
        source_layer <- "counts"
    } else if (length(count_layers) == 1) {
        counts <- get_seurat_layer(
            x = x,
            assay_name = assay_name,
            layer_name = count_layers[[1]]
        )
        source_layer <- count_layers[[1]]
    } else {
        # Merged/split Seurat v5 objects may store counts in counts.1, counts.2, ...
        # Join only those raw-count layers; normalized data layers are untouched.
        if (!"JoinLayers" %in% getNamespaceExports("SeuratObject")) {
            stop(
                "Multiple split raw-count layers were found (",
                paste(count_layers, collapse = ", "),
                "), but the installed SeuratObject version does not provide JoinLayers()."
            )
        }
        joined_object <- SeuratObject::JoinLayers(
            object = x,
            assay = assay_name,
            layers = count_layers,
            new = "counts_for_export"
        )
        counts <- get_seurat_layer(
            x = joined_object,
            assay_name = assay_name,
            layer_name = "counts_for_export"
        )
        source_layer <- paste(count_layers, collapse = "+")
        rm(joined_object)
        gc()
    }

    list(
        counts = counts,
        assay = assay_name,
        source_layer = source_layer
    )
}

extract_seurat_cell_metadata <- function(x, cell_ids) {
    metadata <- x[[]]
    reorder_metadata_rows(metadata, cell_ids, "Seurat cell metadata")
}

extract_seurat_gene_metadata <- function(x, assay_name, feature_ids) {
    assay_object <- x[[assay_name]]
    available_slots <- methods::slotNames(assay_object)

    if ("meta.data" %in% available_slots) {
        metadata <- methods::slot(assay_object, "meta.data")
    } else if ("meta.features" %in% available_slots) {
        metadata <- methods::slot(assay_object, "meta.features")
    } else {
        metadata <- data.frame(row.names = feature_ids)
    }

    reorder_metadata_rows(metadata, feature_ids, "Seurat feature metadata")
}


# -----------------------------------------------------------------------------
# Generic extraction helpers
# -----------------------------------------------------------------------------
extract_object_contents <- function(x, preferred_assay = "RNA") {
    if (inherits(x, "Seurat")) {
        assay_name <- select_seurat_assay(x, preferred_assay)
        count_result <- extract_seurat_counts(x, assay_name)

        count_cells <- colnames(count_result$counts)
        count_genes <- rownames(count_result$counts)

        return(list(
            counts = count_result$counts,
            cells = count_cells,
            genes = count_genes,
            cell_metadata = extract_seurat_cell_metadata(x, count_cells),
            gene_metadata = extract_seurat_gene_metadata(x, assay_name, count_genes),
            assay = count_result$assay,
            source_layer = count_result$source_layer,
            object_type = object_type(x)
        ))
    }

    assay_names <- SummarizedExperiment::assayNames(x)
    if (!"counts" %in% assay_names) {
        stop(
            object_type(x),
            " does not contain an assay named 'counts'. Available assays: ",
            paste(assay_names, collapse = ", "),
            ". Raw counts cannot be reconstructed from normalized data."
        )
    }

    counts <- SummarizedExperiment::assay(x, "counts")

    list(
        counts = counts,
        cells = colnames(x),
        genes = rownames(x),
        cell_metadata = reorder_metadata_rows(
            SummarizedExperiment::colData(x),
            colnames(x),
            paste0(object_type(x), " cell metadata")
        ),
        gene_metadata = reorder_metadata_rows(
            SummarizedExperiment::rowData(x),
            rownames(x),
            paste0(object_type(x), " feature metadata")
        ),
        assay = "counts",
        source_layer = "counts",
        object_type = object_type(x)
    )
}

as_sparse_counts <- function(count_matrix, object_name) {
    if (is.null(rownames(count_matrix)) || is.null(colnames(count_matrix))) {
        stop("The raw-count matrix lacks row or column names for: ", object_name)
    }

    if (!inherits(count_matrix, "sparseMatrix")) {
        count_matrix <- methods::as(count_matrix, "dgCMatrix")
    } else {
        count_matrix <- methods::as(count_matrix, "dgCMatrix")
    }

    values <- count_matrix@x

    if (length(values) > 0) {
        if (anyNA(values) || any(!is.finite(values))) {
            stop("Raw counts contain missing or non-finite values for: ", object_name)
        }
        if (any(values < 0)) {
            stop("Raw counts contain negative values for: ", object_name)
        }

        non_integer <- abs(values - round(values)) > 1e-8
        if (any(non_integer)) {
            examples <- head(values[non_integer], 10)
            stop(
                "The selected counts matrix for '", object_name,
                "' contains non-integer values (examples: ",
                paste(signif(examples, 6), collapse = ", "),
                "). This likely represents normalized expression rather than raw UMI/read counts. ",
                "The script will not export it as raw counts."
            )
        }

        # Preserve sparse structure while ensuring exact whole-number values.
        count_matrix@x <- round(values)
    }

    count_matrix
}


# -----------------------------------------------------------------------------
# Establish the reference gene set from the first object
# -----------------------------------------------------------------------------
first_object_name <- names(object_list)[1]
first_extracted <- extract_object_contents(
    object_list[[1]],
    preferred_assay = preferred_seurat_assay
)
first_count_matrix <- as_sparse_counts(
    first_extracted$counts,
    first_object_name
)
reference_gene_ids <- rownames(first_count_matrix)
reference_gene_metadata <- first_extracted$gene_metadata

if (anyNA(reference_gene_ids) || any(reference_gene_ids == "")) {
    stop("Missing gene identifiers detected in: ", first_object_name)
}
if (anyDuplicated(reference_gene_ids)) {
    duplicate_ids <- unique(reference_gene_ids[duplicated(reference_gene_ids)])
    stop(
        "Duplicate gene identifiers detected in '", first_object_name,
        "'. Examples: ", paste(head(duplicate_ids, 10), collapse = ", ")
    )
}

cat(
    "Reference gene set:",
    format(length(reference_gene_ids), big.mark = ","),
    "unique genes.\n"
)


# -----------------------------------------------------------------------------
# Export gene metadata once
# -----------------------------------------------------------------------------
reference_gene_metadata <- reorder_metadata_rows(
    reference_gene_metadata,
    reference_gene_ids,
    "reference gene metadata"
)
reference_gene_metadata <- clean_metadata(
    reference_gene_metadata,
    "gene metadata"
)
reference_gene_metadata <- rename_reserved_columns(
    reference_gene_metadata,
    c("gene_id")
)

gene_metadata <- cbind(
    data.frame(
        gene_id = reference_gene_ids,
        stringsAsFactors = FALSE,
        check.names = FALSE
    ),
    reference_gene_metadata
)

gene_file <- file.path(output_dir, "genes.tsv")
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


# -----------------------------------------------------------------------------
# Export each object separately
# -----------------------------------------------------------------------------
manifest_rows <- vector("list", length(object_list))
all_exported_cell_ids <- character(0)

for (i in seq_along(object_list)) {
    object_name <- names(object_list)[i]
    safe_name <- safe_object_names[i]
    x <- object_list[[i]]

    if (i == 1) {
        extracted <- first_extracted
        count_matrix <- first_count_matrix
    } else {
        extracted <- extract_object_contents(
            x,
            preferred_assay = preferred_seurat_assay
        )
        count_matrix <- as_sparse_counts(extracted$counts, object_name)
    }

    gene_ids <- rownames(count_matrix)
    if (anyNA(gene_ids) || any(gene_ids == "")) {
        stop("Missing gene identifiers detected in: ", object_name)
    }
    if (anyDuplicated(gene_ids)) {
        duplicate_ids <- unique(gene_ids[duplicated(gene_ids)])
        stop(
            "Duplicate gene identifiers detected in '", object_name,
            "'. Examples: ", paste(head(duplicate_ids, 10), collapse = ", ")
        )
    }
    if (!setequal(gene_ids, reference_gene_ids)) {
        missing_from_object <- setdiff(reference_gene_ids, gene_ids)
        extra_in_object <- setdiff(gene_ids, reference_gene_ids)
        stop(
            "Gene sets differ for '", object_name, "'. Missing reference genes: ",
            paste(head(missing_from_object, 10), collapse = ", "),
            "; extra genes: ", paste(head(extra_in_object, 10), collapse = ", ")
        )
    }

    # Reorder genes to match genes.tsv when feature order differs.
    if (!identical(gene_ids, reference_gene_ids)) {
        count_matrix <- count_matrix[reference_gene_ids, , drop = FALSE]
    }

    original_cell_ids <- colnames(count_matrix)
    if (anyNA(original_cell_ids) || any(original_cell_ids == "")) {
        stop("Missing cell identifiers detected in: ", object_name)
    }
    if (anyDuplicated(original_cell_ids)) {
        duplicate_ids <- unique(original_cell_ids[duplicated(original_cell_ids)])
        stop(
            "Duplicate cell identifiers within '", object_name,
            "'. Examples: ", paste(head(duplicate_ids, 10), collapse = ", ")
        )
    }

    # Prefix cells when exporting multiple objects. This prevents repeated 10x
    # barcodes from colliding while retaining the original barcode separately.
    if (length(object_list) > 1) {
        exported_cell_ids <- paste0(
            sprintf("%02d", i), "_", safe_name, "__", original_cell_ids
        )
    } else {
        exported_cell_ids <- original_cell_ids
    }

    if (any(exported_cell_ids %in% all_exported_cell_ids)) {
        stop("Exported cell identifiers are not globally unique for: ", object_name)
    }
    all_exported_cell_ids <- c(all_exported_cell_ids, exported_cell_ids)

    count_filename <- paste0(
        sprintf("%02d", i), "_", safe_name, "_raw_counts.mtx"
    )
    obs_filename <- paste0(
        sprintf("%02d", i), "_", safe_name, "_obs.tsv"
    )

    count_path <- file.path(output_dir, count_filename)
    obs_path <- file.path(output_dir, obs_filename)

    cat("Exporting:", object_name, "\n")
    cat("  Raw source:", extracted$assay, "/", extracted$source_layer, "\n")

    Matrix::writeMM(count_matrix, file = count_path)

    cell_metadata <- reorder_metadata_rows(
        extracted$cell_metadata,
        original_cell_ids,
        paste0(object_name, " cell metadata")
    )
    cell_metadata <- clean_metadata(cell_metadata, "cell metadata")
    cell_metadata <- rename_reserved_columns(
        cell_metadata,
        c("cell_id", "cell_id_original", "source_object", "source_object_type")
    )

    cell_metadata <- cbind(
        data.frame(
            cell_id = exported_cell_ids,
            cell_id_original = original_cell_ids,
            source_object = rep(object_name, length(exported_cell_ids)),
            source_object_type = rep(extracted$object_type, length(exported_cell_ids)),
            stringsAsFactors = FALSE,
            check.names = FALSE
        ),
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
        source_object = object_name,
        source_object_type = extracted$object_type,
        assay = extracted$assay,
        raw_count_layer = extracted$source_layer,
        counts_file = count_filename,
        obs_file = obs_filename,
        n_cells = ncol(count_matrix),
        n_genes = nrow(count_matrix),
        n_nonzero = length(count_matrix@x),
        total_counts = sum(count_matrix),
        matrix_rows = "genes",
        matrix_columns = "cells",
        matrix_values = "raw_counts",
        stringsAsFactors = FALSE
    )

    cat("  Counts:", count_path, "\n")
    cat("  Metadata:", obs_path, "\n")

    rm(extracted, count_matrix, cell_metadata)
    if (i == 1) {
        first_extracted <- NULL
        first_count_matrix <- NULL
    }
    gc()
}


# -----------------------------------------------------------------------------
# Write manifest and instructions
# -----------------------------------------------------------------------------
manifest <- do.call(rbind, manifest_rows)
manifest_file <- file.path(output_dir, "manifest.tsv")
write.table(
    manifest,
    file = manifest_file,
    sep = "\t",
    quote = TRUE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
)

readme_file <- file.path(output_dir, "README_h5ad_export.txt")
writeLines(
    c(
        "This directory contains raw-count data ready for h5ad assembly.",
        "",
        "Files:",
        "  genes.tsv      Feature metadata. Row order matches every Matrix Market file.",
        "  manifest.tsv   Object names, count files, metadata files, dimensions, and provenance.",
        "  *_raw_counts.mtx  Sparse raw-count matrices stored as genes x cells.",
        "  *_obs.tsv         Cell metadata stored in the same cell order as the matching matrix.",
        "",
        "Important:",
        "  The h5ad builder transposes each matrix to AnnData's cells x genes orientation.",
        "  The final adata.X contains raw, nonnegative, integer-valued counts.",
        "  Normalized/log-transformed values are not exported by this script.",
        "",
        "Build h5ad:",
        "  python build_h5ad_from_export.py --export-dir <this_directory> --output output_raw_counts.h5ad"
    ),
    con = readme_file
)

cat("\nExport complete.\n")
cat("Manifest:", manifest_file, "\n")
cat("Total cells:", format(length(all_exported_cell_ids), big.mark = ","), "\n")
cat("Export directory:", normalizePath(output_dir), "\n")
