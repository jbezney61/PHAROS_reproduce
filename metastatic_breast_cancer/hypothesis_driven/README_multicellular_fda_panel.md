# Multicellular FDA-pair figure

`make_multicellular_fda_panel_figure.py` combines the eight cohort-level breast-cancer
outputs from `pharos hypothesis-driven panel` into publication-ready Panels B, C, and D.

Run it from the cluster project directory after `PC_conv0` through `PC_conv7` finish:

```bash
python metastatic_breast_cancer/hypothesis_driven/make_multicellular_fda_panel_figure.py \
  --runs-root breast_cancer_immune_runs \
  --pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_immune_runs/multicellular_fda_figure \
  --overwrite
```

The `--runs-root` directory must contain:

```text
PC_conv0/  PC_conv1/  PC_conv2/  PC_conv3/
PC_conv4/  PC_conv5/  PC_conv6/  PC_conv7/
```

Alternatively, pass each run explicitly. This is useful when the directories do not
share a parent:

```bash
python metastatic_breast_cancer/hypothesis_driven/make_multicellular_fda_panel_figure.py \
  --run-dir conv0=/path/to/PC_conv0 \
  --run-dir conv1=/path/to/PC_conv1 \
  --run-dir conv2=/path/to/PC_conv2 \
  --run-dir conv3=/path/to/PC_conv3 \
  --run-dir conv4=/path/to/PC_conv4 \
  --run-dir conv5=/path/to/PC_conv5 \
  --run-dir conv6=/path/to/PC_conv6 \
  --run-dir conv7=/path/to/PC_conv7 \
  --pairs-file /path/to/FDA_drug_pairs.csv \
  --output-dir /path/to/multicellular_fda_figure
```

## Evaluation behavior

The upstream panel workflow searches both orders and all available concentration
combinations on batch 0. The figure command therefore uses batches 1 onward by default
so that order/concentration selection and reported performance are separated. Use
`--evaluation-batches all` only for the corresponding all-batch sensitivity analysis.

For each conversion, random-control rows are averaged within random drug pair before
the FDA-pair percentile is calculated. Thus, batches are not treated as independent
random pairs.

## Outputs

- **Panel B:** clustered heatmap of the percentage of random pairs with pair-level mean
  Sinkhorn OT greater than or equal to each FDA pair. Higher values are better.
- **Panel C:** clustered heatmap of the mean percentage of baseline source-to-target
  distance closed. Positive values move toward target; negative values move away.
- **Panel D:** malignant-state reversal versus median immune rescue across the three
  T-cell and two macrophage conversions. It identifies response-profile clusters,
  NK/B compatibility liabilities, and Pareto-optimal combinations.

Every panel is written as PNG, PDF, and SVG by default. The output also contains the
long-format summary, both heatmap matrices, the clustered drug order, the Panel-D
prioritization table, and a machine-readable figure configuration.

The random-control percentiles are empirical descriptive rankings. With 100 random
pairs they have approximately one-percentage-point resolution and should not be
described as multiple-testing-adjusted P values.
