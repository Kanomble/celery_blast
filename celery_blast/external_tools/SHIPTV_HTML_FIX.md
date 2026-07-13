# shiptv Standalone HTML Compatibility Fix

shiptv 0.4.1 emits standalone phylogeny HTML that bundles a modern AG Grid build but initializes it through the removed `new agGrid.Grid(...)` constructor. Generated pages now pass through a shared source-controlled postprocessor that uses the public `agGrid.createGrid(...)` API.

The generated application code no longer accesses AG Grid private internals such as `.context.beans` or `.beanInstance`. Tree-to-grid selection uses row-node public methods, and subtree filtering updates row data with `gridApi.setGridOption('rowData', ...)`.

Metadata column names are promoted to shared generated-code scope so subtree filtering can reuse them after initial grid setup. All-zero Newick branch-length trees are normalized only at render time, preserving the original Newick string while giving Phylocanvas equal positive branch lengths for display.

The postprocessor also fixes invalid header markup, removes the invalid `overflow: none` style, and keeps the HTML standalone by avoiding CDN injection.

Tests added:
- generator-level assertions against real shiptv output;
- all-zero and non-zero Newick normalization coverage;
- a Playwright browser regression test that runs when browser tooling is available in the environment.
