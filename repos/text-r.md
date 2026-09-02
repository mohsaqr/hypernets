# text — CRAN (R, Python wrapper)

- **What**: R interface to Hugging Face transformers via reticulate
  (Kjell et al.) — `textEmbed()` for contextual embeddings by HF model id
  (SimCSE/E5/BGE-class included), plus downstream analysis/prediction
  helpers.
- **Verified 2026-08-24**: on CRAN (presence checked); capability description
  from package documentation background — exact function surface not audited
  this session.
- **Role for us**: an optional external embedding route. sbert is the native
  honets-facing frontend; either package can supply a matrix to
  `knn_hypergraph()`. PLM training remains outside honets.
- **Links**: https://cran.r-project.org/package=text · https://r-text.org
