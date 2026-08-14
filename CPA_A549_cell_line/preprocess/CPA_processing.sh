#CPA A549 2-drug combinatorial perturb-seq

#Data: (combiantorial indexing from CPA paper)
# https://pmc.ncbi.nlm.nih.gov/articles/PMC10258562/#msb202211517-sec-0011
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206741
# 1 set of 2 drugs = Panobinostat and Crizotinib (both are present in our dataset)
# 13 drugs but some of them have drugs included with comparable MOA
# in A549 cells (is included in our dataset)

#norm and log1p the data - and remove bad cells
python prepare_CPA_positive_controls.py

SE_DIR=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_embedding/SE-600M
SE_CKPT=$SE_DIR/se600m_epoch16.ckpt

#run the embedding 
state emb transform \
  --model-folder "$SE_DIR" \
  --checkpoint "$SE_CKPT" \
  --input positive_controls/GSE206741_qc_mad_scrublet_log1p.h5ad \
  --output positive_controls/GSE206741_qc_mad_scrublet_log1p.SE600M.h5ad \
  --embed-key X_state \
  --batch-size 32

#outputs for this positive control
# 16731 genes mapped to embedding file (out of 31854)

