# Random clustering baseline (R)

This module implements Question 4 of STA426 Exercise 06. It reuses the
benchmark's `r_osta` environment. Run it through Terminal/OmniBenchmark,
not by pasting `run.R` into the RStudio Console: it expects command-line inputs.

Each cell receives an independent, uniformly sampled label from 1 to k using
`set.seed(42)`. The value k comes from `data.clusters_truth_num`; individual
reference labels and expression values are not used to choose the labels.

The existing metrics script requires `leiden` and `louvain` columns, so both
contain the SAME random assignments. Neither Leiden nor Louvain is run here.
Consequently, the reported `agreement.ari.leiden` is the RANDOM baseline ARI.

For the other metrics, the module computes real PCA coordinates from
library-size-normalized, log-transformed counts using up to 1,500 variable
genes and 30 components. These coordinates do not influence random labels.
Cells with zero or invalid library sizes are excluded; no other QC is applied.
This may use a different cell subset from Scanpy's manual filtering, so do
not describe the pipelines as having identical preprocessing. ARI is expected
to be near zero for random assignments, but its actual value must be measured.

## 1. Save and test the module

Unzip the download and place the `omni-random-r` folder in
`/Users/meretaltorfer/STA426/`, next to `omni-pipelines-plan`.
Do not put it in the benchmark's generated `out/.modules` directory.

Run the following in Terminal. This direct test uses your existing R environment
and writes four files to `test-output` without changing the benchmark plan.

```bash
cd "$HOME/STA426/omni-random-r"
chmod +x run.R
env -u RHOME -u R_HOME PATH="$HOME/conda-envs/sta426-ex06/e4db74c7467820a10bcb6a7add30e1a4_/bin:/usr/bin:/bin:/usr/sbin:/sbin" "$HOME/conda-envs/sta426-ex06/e4db74c7467820a10bcb6a7add30e1a4_/bin/Rscript" run.R --output_dir test-output --name datasets-r --data.h5ad "$HOME/STA426/omni-pipelines-plan/out/datasets/datasets-r/.6fbd9e19/datasets-r.h5ad" --data.clusters_truth_num "$HOME/STA426/omni-pipelines-plan/out/datasets/datasets-r/.6fbd9e19/datasets-r.clusters_truth_num.txt" --seed 42 --n_comp 30
```

The output should end with `Wrote random baseline outputs ...`. The script
uses the same Python H5AD reader as your successful existing R pipeline.
The code was checked against the supplied interfaces, but has not been run
on your data by the assistant. If this test fails, keep its full error message.

## 2. Put the module on GitHub

Create an empty GitHub repository named `omni-random-r` in your account.
Do not initialize it with another README, license, or gitignore. Follow your
course's visibility rules. For a private repository, your Mac's Git must be
authenticated with permission to clone it before OmniBenchmark can use it.

In Terminal, inside `~/STA426/omni-random-r`, run:

```bash
git init
git add run.R omnibenchmark.yaml CITATION.cff LICENSE README.md benchmark-snippet.yaml .gitignore
git commit -m "Add reproducible random clustering baseline"
git branch -M main
```

Replace `YOUR_GITHUB_USERNAME` below with your actual GitHub username:

```bash
git remote add origin https://github.com/YOUR_GITHUB_USERNAME/omni-random-r.git
git push -u origin main
git rev-parse HEAD
```

Copy the full commit hash printed by the last command.

## 3. Add the module to benchmark.yaml

Open `~/STA426/omni-pipelines-plan/benchmark.yaml` in your editor.
Copy the stanza in `benchmark-snippet.yaml` into the `methods` stage's
`modules` list, alongside `osta-r` and `scanpy`, before the `metrics` stage.
Replace the GitHub username and commit placeholder. Keep your Question 3
Scanpy `n_comp: [30, 50]` change. Do not change the stage outputs.

`omnibenchmark.yaml` in this module folder identifies `run.R` as its entrypoint.
It is different from the main benchmark plan, `benchmark.yaml`.

## 4. Run the extended benchmark

In Terminal:

```bash
conda activate omni
cd "$HOME/STA426/omni-pipelines-plan"
ob validate plan benchmark.yaml
env -u RHOME -u R_HOME PATH="$CONDA_PREFIX/bin:$HOME/anaconda3/condabin:/usr/bin:/bin:/usr/sbin:/sbin" "$CONDA_PREFIX/bin/ob" run benchmark.yaml -c 2 -- --conda-prefix "$HOME/conda-envs/sta426-ex06"
```

Completed existing outputs can be reused. The new random module and its
metrics should run; a full rerun with `--forceall` is not needed.

## 5. Retrieve the actual ARI in your QMD

Put this R code inside an R chunk in your Question 4 answer. `all.files = TRUE`
includes OmniBenchmark's parameter directories whose names start with a dot.

```r
out_dir <- path.expand("~/STA426/omni-pipelines-plan/out")
metric_files <- list.files(
  out_dir, pattern = "\\.metrics\\.json$",
  recursive = TRUE, full.names = TRUE, all.files = TRUE
)
random_files <- metric_files[
  grepl("/methods/random-r/", metric_files, fixed = TRUE)
]
stopifnot(length(random_files) > 0L)

random_results <- do.call(rbind, lapply(random_files, function(path) {
  result <- jsonlite::fromJSON(path)
  data.frame(
    method = "Random baseline",
    seed = result$clustering_info$seed,
    ARI = result$agreement$ari$leiden,
    metrics_file = path
  )
}))
random_results
```

Explain in your answer that the existing metrics interface calls the column
`leiden`, but this module fills it with random labels. Report the measured
ARI and compare it with the Scanpy result from Questions 2/3. An ARI near zero
means agreement close to chance; negative values are possible.

Submit your updated benchmark.yaml, QMD, and rendered HTML to your private
exercise repository as requested. Keep the module code in its own GitHub
repository referenced by the benchmark plan.
